use clap::{Args, Command, Parser, Subcommand};
use clap_complete::{generate, Generator, Shell};
use std::path::PathBuf;

/// pN/pS Utilities
#[derive(Parser, Debug)]
#[command(author, version, about, arg_required_else_help(true))]
pub struct Cli {
    /// Generates Shell completion code
    ///
    /// It prints the code to the standard output and the way to
    /// use depends on the Shell. For Fish, redirect to a file
    /// with `.fish` extension in `~/.config/fish/completion`.
    #[arg(long)]
    pub complete: Option<Shell>,

    #[command(subcommand)]
    pub command: Option<Commands>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    Config(Config),
    Parse(Parse),
    Calc(Calc),
}

/// Generates the config file for command `parse`
///
/// The config file contains information to map
/// the BAM files used to generate the VCF file
/// and the Depth files. Another column is used
/// as Sample names, to be used in the CSV file
/// that is the final output.
#[derive(Args, Debug)]
pub struct Config {
    /// VCF file to get the BAM files used
    ///
    /// This is used as a column in the config file
    /// so it can be compared to the depth files
    /// and the sample names
    ///
    /// The file will contains three `tab` separated columns
    /// Sample ID, VCF Column, Depth file
    #[arg(short, long, required = true)]
    pub vcf_file: PathBuf,
    /// Output to file, instead of stdout
    #[arg(short, long)]
    pub output_file: Option<PathBuf>,
    /// Depth files
    ///
    /// File names will be checked against the BAM
    /// files from the VCF file.
    pub depth_files: Vec<PathBuf>,
}

/// Parse a VCF file and saves data to calculate pN/pS
///
/// The VCF file is expected to be created with
/// `samtools` and `bcftools`, only the `CDS` feature type
/// is used
#[derive(Args, Debug)]
pub struct Parse {
    /// This file can be created with the `config` command.
    /// Config file
    ///
    #[arg(short, long, required = true)]
    pub config_file: PathBuf,
    // The GFF with annotations
    #[arg(short, long, required = true)]
    pub gff_file: PathBuf,
    /// The Fasta file
    #[arg(short, long, required = true)]
    pub fasta_file: PathBuf,
    /// Minimum accepted coverage, corresponding to `DP` in VCF
    #[arg(short, long, default_value_t = 4, value_parser = clap::value_parser!(u32).range(1..=20))]
    pub min_depth: u32,
    /// Minimum read coverage from BAM file
    #[arg(short = 'a', long, default_value_t = 4, value_parser = clap::value_parser!(u32).range(1..=20))]
    pub min_coverage: u32,
    /// Minimum Quality `QUAL` in VCF file
    #[arg(short = 'q', long, default_value_t = 30.)]
    pub min_qual: f64,
    /// VCF file with SNPs
    pub vcf_file: PathBuf,
    /// file name for the output, defaults to `pnps.json.gz`
    pub output_file: Option<PathBuf>,
}

/// Uses the result of `parse` and map files to calculate
/// pN/pS
#[derive(Args, Debug)]
pub struct Calc {
    /// Gene map, mapping a UID to another ID
    #[arg(short, long)]
    pub gene_map: Option<PathBuf>,
    /// Taxonomy file, use `taxa-utils`
    #[arg(short, long, group = "taxa")]
    pub taxonomy: Option<PathBuf>,
    /// Taxon map, mapping a UID to a taxon ID
    /// 
    /// It will resolve to the lineage string using the taxonomy supplied
    #[arg(short = 'm', long, requires = "taxonomy", group = "lineage")]
    pub taxon_map: Option<PathBuf>,
    /// Alternative to `--taxon_map` and the map contains strings showing the full lineage
    #[arg(short = 'l', long, group = "lineage")]
    pub lineage_map: Option<PathBuf>,
    /// Taxon rank to map taxa from the map (not implemented)
    #[arg(short = 'r', long, requires = "taxon_map")]
    pub taxon_rank: Option<String>,
    /// Only save pS value, not pN/pS
    #[arg(short = 's', long, group = "split")]
    pub output_ps: bool,
    /// Only save pN value, not pN/pS
    #[arg(short = 'n', long, group = "split")]
    pub output_pn: bool,
    /// Taxonomy file, use `taxa-utils` `import` or `download` to create
    pub input_file: PathBuf,
    /// Output file
    pub output_file: PathBuf,
}

/// Generates the completion for the specified shell
///
/// Slightly modified from example
pub fn print_completions<G: Generator>(gen: G, cmd: &mut Command) {
    generate(gen, cmd, cmd.get_name().to_string(), &mut std::io::stdout());
}
