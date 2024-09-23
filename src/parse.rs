use super::cli::Parse;
use anyhow::{bail, Result};
use bio_rascal::fasta::FastaReader;
use bio_rascal::gff::{Annotation, GffReader};
use bio_rascal::samtools::read_depth_file;
use bio_rascal::sequence::SequenceRecord;
use bio_rascal::snps::PnPs;
use console::style;
use indicatif::ProgressBar;
use log::{error, info};
use serde_json::to_writer;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use uuid::Uuid;

pub type SampleInfo = HashMap<String, (String, String)>;
pub type SamplePnPs = HashMap<String, HashMap<Uuid, PnPs>>;

fn read_config_file<P: AsRef<Path>>(file_name: P) -> Result<SampleInfo> {
    info!("Reading Config file: {}", file_name.as_ref().display());
    let reader = BufReader::new(File::open(file_name)?);

    let mut sample_info = SampleInfo::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<_> = line.trim().split('\t').collect();
        if fields.len() < 3 {
            bail!(
                "Cannot parse sample information, expect 3 columns, got {}",
                fields.len()
            );
        }
        sample_info.insert(fields[1].into(), (fields[0].into(), fields[2].into()));
    }

    Ok(sample_info)
}

fn read_gff_file<P: AsRef<Path>>(file_name: &P) -> Result<HashMap<Uuid, Annotation>> {
    let mut annotations: HashMap<Uuid, Annotation> = HashMap::new();

    for annotation in GffReader::new(file_name).unwrap() {
        match annotation.feature_type.as_str() {
            "CDS" => annotations.insert(annotation.uid, annotation),
            _ => continue,
        };
    }
    Ok(annotations)
}

fn read_fasta_file<P: AsRef<Path>>(file_name: P) -> Result<HashMap<String, SequenceRecord>> {
    let reader = FastaReader::new(file_name)?;
    let mut hm: HashMap<String, SequenceRecord> = HashMap::new();

    for record in reader {
        _ = hm.insert(record.id.clone(), record);
    }

    Ok(hm)
}

/// For one sample
fn prepare_annotations(
    annotations: &HashMap<Uuid, Annotation>,
    fasta_records: &HashMap<String, SequenceRecord>,
) -> Result<Vec<PnPs>> {
    let mut pnps_list: Vec<PnPs> = vec![];
    info!("Preparing pN/pS data");
    let pb = ProgressBar::new(annotations.len() as u64);

    for annotation in annotations.values() {
        let record = match fasta_records.get(&annotation.seq_id) {
            None => bail!("Cannot find sequence for {}", &annotation.seq_id),
            Some(value) => value,
        };
        let (exp_syn, exp_nonsyn) = annotation.get_exp_syn(&record.seq);
        let pnps = PnPs {
            uid: annotation.uid,
            exp_syn,
            exp_nonsyn,
            ..Default::default()
        };
        pnps_list.push(pnps);
        pb.inc(1);
    }

    Ok(pnps_list)
}

fn prepare_sample_pnps<P: AsRef<Path>>(
    pnps_base: &[PnPs],
    depth_file: P,
    annotations: &HashMap<Uuid, Annotation>,
    min_cov: u32,
) -> Result<HashMap<Uuid, PnPs>> {
    let dm = read_depth_file(depth_file)?;

    let mut sp: HashMap<Uuid, PnPs> = HashMap::with_capacity(pnps_base.len());

    for pnps in pnps_base {
        let a = match annotations.get(&pnps.uid) {
            None => bail!("Cannot find annotation {}", &pnps.uid),
            Some(value) => value,
        };
        let mut p = PnPs {
            uid: a.uid,
            exp_nonsyn: pnps.exp_nonsyn,
            exp_syn: pnps.exp_syn,
            ..Default::default()
        };
        p.coverage = match dm.get(&a.seq_id) {
            None => 0u32,
            Some(depth) => depth.coverage_at(&a.start, &a.end),
        };
        if p.coverage >= min_cov {
            sp.insert(p.uid, p);
        }
    }
    Ok(sp)
}

fn add_depth_sample_data(
    sample_info: &SampleInfo,
    pnps_list: &[PnPs],
    annotations: &HashMap<Uuid, Annotation>,
    min_cov: u32,
) -> Result<SamplePnPs> {
    let pb = indicatif::ProgressBar::new(sample_info.len() as u64);

    let mut pnps_map: HashMap<String, HashMap<Uuid, PnPs>> =
        HashMap::with_capacity(sample_info.len());

    for (_v, (sample_id, d)) in sample_info.iter() {
        pb.println(format!(
            "[+] Reading depth information for sample {} from file {}",
            style(sample_id).blue(),
            style(d).blue()
        ));
        let sp = prepare_sample_pnps(pnps_list, d, annotations, min_cov)?;

        let max_pnps = sp.values().max_by_key(|c| c.coverage).unwrap();
        pb.println(format!(
            " | -> Max sample coverage {} in annotation: {}",
            style(max_pnps.coverage).yellow(),
            style(max_pnps.uid).yellow()
        ));
        pnps_map.insert(sample_id.clone(), sp);
        pb.inc(1);
    }

    Ok(pnps_map)
}

fn parse_vcf_file<P: AsRef<Path>>(
    file_name: P,
    pnps_map: &mut SamplePnPs,
    fasta_records: &HashMap<String, SequenceRecord>,
    annotations: &HashMap<Uuid, Annotation>,
    sample_info: &SampleInfo,
    min_qual: f64,
    min_depth: u32
) -> Result<()> {
    info!("Preparing annotations");
    let mut ann_seq: HashMap<&String, Vec<&Annotation>> = HashMap::new();
    for annotation in annotations.values() {
        ann_seq
            .entry(&annotation.seq_id)
            .and_modify(|e| e.push(annotation))
            .or_insert(vec![annotation]);
    }

    let vcf_reader = bio_rascal::snps::VcfReader::new(file_name)?;
    info!("Number of VCF samples: {}", vcf_reader.sample_names.len());

    let pb = indicatif::ProgressBar::new_spinner().with_message("VCF Reading");
    let mut count = 0u32;
    let mut skipped_dp = 0u32;
    let mut skipped_indel = 0u32;
    let mut skipped_qual = 0u32;

    for record in vcf_reader {
        pb.inc(1);
        count += 1;
        if record.info.dp < min_depth {
            skipped_dp += 1;
            continue;
        } else if record.qual < min_qual {
            skipped_qual += 1;
            continue;
        } else if record.info.indel || record.ref_c.len() > 1 {
            skipped_indel += 1;
            continue;
        }

        let ann = match ann_seq.get(&record.chrom) {
            None => {
                //error!("{}", record.chrom);
                continue;
            }
            Some(value) => value.iter().filter(|a| a.contains(record.pos)), //.collect()
        };
        for a in ann {
            if let Some(seqr) = fasta_records.get(&a.seq_id) {
                for (sample_id, alt) in record.get_sample_snps() {
                    let sample_id = match sample_info.get(&sample_id) {
                        None => {
                            error!("Cannot find the sample {sample_id}");
                            continue;
                        }
                        Some(value) => &value.0,
                    };
                    if let Some(sample_pnps_map) = pnps_map.get_mut(sample_id) {
                        if let Some(mut sample_pnps) = sample_pnps_map.get_mut(&a.uid) {
                            match a.is_syn(&seqr.seq, record.pos, &alt) {
                                Ok(is_syn) => {
                                    if is_syn {
                                        sample_pnps.syn += 1;
                                    } else {
                                        sample_pnps.nonsyn += 1;
                                    }
                                },
                                Err(err) => pb.println(style(err).red().to_string()),
                            }
                        }
                    }
                    //info!("{} -> {}", sample_id, alt);
                }
            }
        }
        //info!("{} -> {}, {}", record.chrom, ann.len(), record.info.ac.len());
    }

    info!(
        "VCF records {count}, Skipped INDEL: {skipped_indel}, Skipped for low QUAL: {skipped_qual}, Skipped for low DP (depth) {:.2}%",
        skipped_dp as f64 / count as f64 * 100f64
    );

    Ok(())
}

pub fn parse_command(options: Parse) -> Result<()> {
    let output_file = match options.output_file {
        None => File::create("pnps.json")?,
        Some(value) => File::create(value)?,
    };
    
    info!("Minimum Depth {}, Qual {}, Coverage {}", options.min_depth, options.min_qual, options.min_coverage);

    // starts reading the GFF file
    let annotations = read_gff_file(&options.gff_file)?;
    info!("Number of Annotations: {}", annotations.len());
    let sample_info = read_config_file(&options.config_file)?;
    info!("Number of Samples in Config file: {}", sample_info.len());
    let fasta_records = read_fasta_file(&options.fasta_file)?;
    info!("Number of Fasta records: {}", fasta_records.len());
    let pnps_list = prepare_annotations(&annotations, &fasta_records)?;

    let mut pnps_map =
        add_depth_sample_data(&sample_info, &pnps_list, &annotations, options.min_coverage)?;
    parse_vcf_file(
        options.vcf_file,
        &mut pnps_map,
        &fasta_records,
        &annotations,
        &sample_info,
        options.min_qual,
        options.min_depth
    )?;

    to_writer(output_file, &pnps_map)?;

    Ok(())
}
