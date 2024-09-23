use super::parse::SamplePnPs;
use anyhow::{Context, Result};
use bio_rascal::snps::{CalculatePnPs, GroupPnPs};
use bio_rascal::taxon::Taxonomy;
use log::{info, warn};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;
use std::str::FromStr;
use uuid::Uuid;

// TODO: account for situation where -g is passed but not `-l` or a taxonomy is

#[allow(non_camel_case_types)]
enum ResultType {
    pNpS,
    pN,
    pS,
}

type GeneMap = HashMap<Uuid, Vec<String>>;
type TaxonMap = HashMap<Uuid, u32>;
type LineageMap = HashMap<Uuid, String>;
type SampleGroupPnPs<'a> = HashMap<String, HashMap<(String, u32, String), GroupPnPs<'a>>>;

fn group_pnps<'a>(
    pnps_map: &'a SamplePnPs,
    gene_map: &GeneMap,
    taxon_map: &TaxonMap,
    lineage_map: &LineageMap,
) -> SampleGroupPnPs<'a> {
    let mut grouped_pnps = SampleGroupPnPs::with_capacity(pnps_map.len());
    for (sample_id, pnps_values) in pnps_map.iter() {
        let mut sample_map: HashMap<(String, u32, String), GroupPnPs<'a>> =
            HashMap::with_capacity(pnps_values.len());
        for (uid, value) in pnps_values.iter() {
            // make Keys
            let gene_ids = match gene_map.get(uid) {
                None => vec![uid.to_string()],
                Some(values) => values.clone(),
            };
            let taxon_id = match taxon_map.get(uid) {
                None => 0u32,
                Some(taxon_id) => *taxon_id,
            };
            let taxon_lineage = match lineage_map.get(uid) {
                None => "".to_string(),
                Some(taxon_lineage) => taxon_lineage.clone(),
            };
            // Use Entry to insert/modify
            for gene_id in gene_ids {
                let key = (gene_id.clone(), taxon_id, taxon_lineage.clone());
                sample_map
                    .entry(key)
                    .and_modify(|v| v.pnps.push(value))
                    .or_insert(GroupPnPs {
                        gene_id,
                        taxon_id,
                        pnps: vec![value],
                        taxon_lineage: taxon_lineage.clone(),
                    });
            }
        }
        grouped_pnps.insert(sample_id.clone(), sample_map);
    }

    grouped_pnps
}

fn write_grouped_output<P: AsRef<Path>>(
    file_name: P,
    pnps_map: &SampleGroupPnPs,
    result_type: &ResultType,
    taxonomy: &Taxonomy,
) -> Result<()> {
    info!("Writing results to file {}", file_name.as_ref().display());

    let mut non_null_index: HashSet<&(String, u32, String)> = HashSet::new();

    // prepares the index to be able to write records
    for map in pnps_map.values() {
        for key in map.keys() {
            non_null_index.insert(key);
        }
    }

    let mut writer = csv::Writer::from_path(file_name).context("Problem opening file")?;
    let mut record = Vec::with_capacity(pnps_map.len() + 2);

    record.push("gene_id".to_string());
    record.push("taxon".to_string());
    record.push("lineage".to_string());
    for sample_id in pnps_map.keys() {
        record.push(sample_id.clone());
    }
    writer
        .write_record(&record)
        .context("Problem writing Header")?;

    for (gene_id, taxon_id, lineage) in non_null_index {
        record.clear();
        record.push(gene_id.clone());
        record.push(taxon_id.to_string());
        if lineage.is_empty() {
            record.push(
                taxonomy
                    .get_taxon_lineage_string(taxon_id)
                    .context("Cannot build lineage string")?,
            );
        } else {
            record.push(lineage.clone());
        }
        let mut values: Vec<f64> = Vec::with_capacity(pnps_map.len());

        for pmap in pnps_map.values() {
            let p = match pmap.get(&(gene_id.clone(), *taxon_id, lineage.clone())) {
                None => f64::NAN,
                Some(value) => match result_type {
                    ResultType::pNpS => value.get_pnps(),
                    ResultType::pN => value.get_pn(),
                    ResultType::pS => value.get_ps(),
                },
            };
            // push value first
            values.push(p);
        }
        // counts how many "normal", so no infinite, NAN are inside and
        // only convert and write if at least one is there.
        if values.iter().filter(|e| e.is_normal()).count() > 0 {
            record.extend(values.iter().map(|e| e.to_string()));
            writer
                .write_record(&record)
                .context("Problem writing Record")?;
        }
    }

    writer.flush().context("Problem flushing to disk")?;

    Ok(())
}

fn write_output<P: AsRef<Path>>(
    file_name: P,
    pnps_map: &SamplePnPs,
    result_type: &ResultType,
) -> Result<()> {
    info!("Writing results to file {}", file_name.as_ref().display());

    let mut non_null_index: HashSet<Uuid> = HashSet::new();

    for uids in pnps_map.values() {
        for uid in uids.keys() {
            non_null_index.insert(*uid);
        }
    }

    let mut writer = csv::Writer::from_path(file_name).context("Problem opening file")?;
    let mut record = Vec::with_capacity(pnps_map.len() + 1);
    record.push("uid".to_string());

    for sample_id in pnps_map.keys() {
        record.push(sample_id.clone());
    }

    writer
        .write_record(&record)
        .context("Problem writing Header")?;

    for uid in non_null_index {
        record.clear();
        record.push(uid.to_string());
        let mut values: Vec<f64> = Vec::with_capacity(pnps_map.len());

        for pmap in pnps_map.values() {
            let p = match pmap.get(&uid) {
                None => f64::NAN,
                Some(value) => match result_type {
                    ResultType::pNpS => value.get_pnps(),
                    ResultType::pN => value.get_pn(),
                    ResultType::pS => value.get_ps(),
                },
            };
            // push value first
            values.push(p);
        }
        // counts how many "normal", so no infinite, NAN are inside and
        // only convert and write if at least one is there.
        if values.iter().filter(|e| e.is_normal()).count() > 0 {
            record.extend(values.iter().map(|e| e.to_string()));
            writer
                .write_record(&record)
                .context("Problem writing Record")?;
        }
    }

    writer.flush().context("Problem flushing to disk")?;

    Ok(())
}

fn read_gene_map_file<P: AsRef<Path>>(file_name: P) -> Result<GeneMap> {
    info!("Reading Gene map file: {}", &file_name.as_ref().display());
    let file_handle = bio_rascal::io::open_file(file_name).context("Cannot open file")?;

    let mut gene_map = GeneMap::new();

    for line in file_handle.lines() {
        let line = line.context("Cannot parse line in map file")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some((uid, values)) = line.trim().split_once('\t') {
            let uid = Uuid::from_str(uid).context("Cannot parse UID in line")?;
            let values: Vec<String> = values.split(',').map(|s| s.to_string()).collect();
            gene_map.insert(uid, values);
        }
    }
    Ok(gene_map)
}

fn read_taxon_map_file<P: AsRef<Path>>(file_name: P) -> Result<TaxonMap> {
    info!("Reading Taxon map file: {}", &file_name.as_ref().display());
    let file_handle = bio_rascal::io::open_file(file_name).context("Cannot open file")?;
    let mut taxon_map = TaxonMap::new();

    for line in file_handle.lines() {
        let line = line.context("Problem parsing line")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some((uid, taxon_id)) = line.trim().split_once('\t') {
            let uid = Uuid::from_str(uid).context("Cannot parse UID in line")?;
            let taxon_id = u32::from_str(taxon_id).context("Cannot convert taxon ID")?;
            taxon_map.insert(uid, taxon_id);
        }
    }

    Ok(taxon_map)
}

fn read_lineage_map_file<P: AsRef<Path>>(file_name: P) -> Result<LineageMap> {
    info!("Reading Taxon map file: {}", &file_name.as_ref().display());
    let file_handle = bio_rascal::io::open_file(file_name).context("Cannot open file")?;
    let mut lineage_map = LineageMap::new();

    for line in file_handle.lines() {
        let line = line.context("Problem parsing line")?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some((uid, lineage)) = line.trim().split_once('\t') {
            let uid = Uuid::from_str(uid).context("Cannot parse UID in line")?;
            lineage_map.insert(uid, lineage.to_string());
        }
    }

    Ok(lineage_map)
}

pub fn calc_command(options: super::cli::Calc) -> Result<()> {
    let mut gene_map = GeneMap::new();
    if let Some(gene_map_file) = options.gene_map {
        gene_map = read_gene_map_file(gene_map_file)?;
    }

    let mut taxonomy = Taxonomy::default();
    if let Some(taxonomy_file) = options.taxonomy {
        taxonomy = bio_rascal::taxon::Taxonomy::read_from_file(taxonomy_file)
            .context("Problem reading taxonomy file")?;
    }

    let mut taxon_map = TaxonMap::new();
    if let Some(taxon_map_file) = options.taxon_map {
        taxon_map = read_taxon_map_file(taxon_map_file)?;
    }

    let mut lineage_map = LineageMap::new();
    if let Some(lineage_map_file) = options.lineage_map {
        lineage_map = read_lineage_map_file(lineage_map_file)?;
    }

    if let Some(taxon_rank) = options.taxon_rank {
        warn!("Using a rank is not implemented, passed: {}", taxon_rank);
    }

    info!(
        "Reading pN/pS data from file: {}",
        options.input_file.display()
    );
    let pnps_file = bio_rascal::io::open_file_base(&options.input_file)
        .context("Cannot open the input file")?;
    let pnps_map: SamplePnPs =
        serde_json::from_reader(pnps_file).context("Problem parsing the input file")?;

    let result_type = match (options.output_pn, options.output_ps) {
        (true, false) => {
            info!("Calculating pN");
            ResultType::pN
        }
        (false, true) => {
            info!("Calculating pS");
            ResultType::pS
        }
        _ => {
            info!("Calculating pN/pS");
            ResultType::pNpS
        }
    };

    if taxon_map.is_empty() && gene_map.is_empty() && lineage_map.is_empty() {
        write_output(&options.output_file, &pnps_map, &result_type)
            .context("Problem writing output file")?;
    } else {
        let grouped_pnps = group_pnps(&pnps_map, &gene_map, &taxon_map, &lineage_map);
        write_grouped_output(
            &options.output_file,
            &grouped_pnps,
            &result_type,
            &taxonomy,
        )
        .context("Problem writing output file")?;
    }

    Ok(())
}
