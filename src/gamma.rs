use clap::ArgMatches;
use std::error::Error;

pub fn callback(_sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    Ok(())
}
