extern crate clap;
extern crate crossbeam;
extern crate csv;
extern crate indicatif;
extern crate pretty_env_logger;

#[macro_use]
extern crate log;

extern crate carina;
extern crate rand;
extern crate sce;

use clap::{App, Arg, SubCommand};
use std::error::Error;

mod spatial;

fn main() -> Result<(), Box<dyn Error>> {
    let matches = App::new("indus")
        .version("0.1.0")
        .author("Avi Srivastava")
        .about("Generate summary stats for multimodal data.")
        .subcommand(
            SubCommand::with_name("autocorr")
                .about("A subcommand to generate auto-correlation summary statistics.")
                .arg(
                    Arg::with_name("weights")
                        .long("weights")
                        .short("w")
                        .takes_value(true)
                        .required(true)
                        .help("path to the weight matrix."),
                )
                .arg(
                    Arg::with_name("values")
                        .long("values")
                        .short("v")
                        .takes_value(true)
                        .required(true)
                        .help("path to the value matrix."),
                )
                .arg(
                    Arg::with_name("method")
                        .long("method")
                        .short("m")
                        .takes_value(true)
                        .required(true)
                        .possible_values(&["Moransi", "Gearyc"]),
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output file."),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    if let Some(sub_m) = matches.subcommand_matches("autocorr") {
        spatial::callback(&sub_m)?
    }

    Ok(())
}
