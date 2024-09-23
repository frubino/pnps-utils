mod calc;
mod cli;
mod config;
mod parse;
mod utils;

use anyhow::Result;
use calc::calc_command;
use clap::{CommandFactory, Parser};
use cli::print_completions;
use config::config_command;
use env_logger::Env;
use log::{error, info};
use parse::parse_command;

fn main() -> Result<()> {
    let cli = cli::Cli::parse();

    if let Some(generator) = cli.complete {
        let mut cmd = cli::Cli::command();
        print_completions(generator, &mut cmd);
    } else if let Some(command) = cli.command {
        // Starts the logging
        // possible to define log level with RUST_LOG
        env_logger::Builder::from_env(Env::default().default_filter_or("info"))
            .format_timestamp_millis()
            .init();

        let result = match command {
            cli::Commands::Config(options) => config_command(options),
            cli::Commands::Parse(options) => parse_command(options),
            cli::Commands::Calc(options) => calc_command(options),
            //_ => todo!(),
        };

        match result {
            Err(err) => error!("{}", err),
            Ok(_) => info!("Finished execution"),
        }
    }

    /* let base_dir = PathBuf::from(var("CARGO_MANIFEST_DIR").unwrap());

    */
    Ok(())
}
