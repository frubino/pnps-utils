use anyhow::{bail, Result};
use log::{error, info};
use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

/// Checks passed Option<Path> and returns an open file, if the Option is
/// `None` in which case the `stdout` is used.
pub fn file_or_stdout<P: AsRef<Path>>(output_file: &Option<P>) -> Result<Box<dyn Write>> {
    let result = match output_file {
        None => {
            info!("Opening stdout");
            Box::new(std::io::stdout()) as Box<dyn Write>
        }
        Some(value) => match File::create(value) {
            Err(err) => {
                error!("Cannot create file {}", value.as_ref().display());
                bail!("{}", err.to_string())
            }
            Ok(handle) => {
                info!("Opening file: {:?}", value.as_ref().display());
                Box::new(handle)
            }
        },
    };
    Ok(result)
}

/// Checks passed Option<Path> and returns an open file, if the Option is
/// `None` in which case the `stdin` is used.
pub fn file_or_stdin<P: AsRef<Path>>(input_file: &Option<P>) -> Result<Box<dyn Read>> {
    let result = match input_file {
        None => {
            info!("Opening stdin");
            Box::new(std::io::stdin().lock()) as Box<dyn Read>
        }
        Some(value) => bio_rascal::io::open_file_base(value)?,
    };
    Ok(result)
}
