use super::cli::Config;
use super::utils::file_or_stdout;
use anyhow::{bail, Result};
use log::info;
use std::io::Write;

use std::path::Path;

static HEADER: &str =
    "#Rearrange to make files and columns correspond\n#SAMPLE_ID\tVCF_COLUMN\tDEPTH_FILE";

fn write_config_file<R: Write, P: AsRef<Path>>(
    file_handle: &mut R,
    sample_ids: &[String],
    vcf_samples: &[String],
    depth_files: &[P],
) -> Result<()> {
    writeln!(file_handle, "{}", &HEADER)?;
    for (index, sample_id) in sample_ids.iter().enumerate() {
        writeln!(
            file_handle,
            "{}\t{}\t{}",
            sample_id,
            vcf_samples[index],
            depth_files[index].as_ref().display()
        )?;
    }

    file_handle.flush()?;
    Ok(())
}

pub fn config_command(options: Config) -> Result<()> {
    let mut ouput_file = file_or_stdout(&options.output_file)?;

    let vcf_reader = bio_rascal::snps::VcfReader::new(&options.vcf_file)?;
    let vcf_samples = vcf_reader.sample_names.clone();
    drop(vcf_reader);

    let sample_ids: Vec<String> = vcf_samples
        .iter()
        .flat_map(|c| Path::new(c).file_stem())
        .flat_map(|e| e.to_str())
        .map(|e| e.to_string())
        .collect();

    if options.depth_files.len() != sample_ids.len() {
        bail!(
            "Length of samples ({}) in VCF file and number of Depth files ({}) is not the same",
            sample_ids.len(),
            options.depth_files.len()
        );
    }
    info!("Writing config");

    write_config_file(
        &mut ouput_file,
        &sample_ids,
        &vcf_samples,
        &options.depth_files,
    )?;

    Ok(())
}
