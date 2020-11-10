use std::error::Error;
use clap::ArgMatches;

pub fn callback(_sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    Ok(())
}