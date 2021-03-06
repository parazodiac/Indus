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

mod configs;
mod gibbs;
mod links;
mod multimodal;
mod spatial;
mod unify;

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
        .subcommand(
            SubCommand::with_name("gamma")
                .about("A subcommand to generate gamma fields.")
                .arg(
                    Arg::with_name("ipaths")
                        .long("ipaths")
                        .short("i")
                        .takes_value(true)
                        .required(true)
                        .multiple(true)
                        .help("path to the parent folders of matrices."),
                )
                .arg(
                    Arg::with_name("links")
                        .long("links")
                        .short("l")
                        .takes_value(true)
                        .required(true)
                        .help("path to the file with feature links."),
                )
                .arg(
                    Arg::with_name("anchors")
                        .long("anchors")
                        .short("a")
                        .takes_value(true)
                        .help("path to the file with cellular barcode anchors."),
                )
                .arg(
                    Arg::with_name("microclusters")
                        .long("microclusters")
                        .short("m")
                        .takes_value(true)
                        .help("path to the file with microclusters of pivot assay."),
                )
                .arg(
                    Arg::with_name("output")
                        .long("output")
                        .short("o")
                        .takes_value(true)
                        .required(true)
                        .help("path to the output path file."),
                ),
        )
        .get_matches();
    pretty_env_logger::init_timed();

    if let Some(sub_m) = matches.subcommand_matches("gamma") {
        unify::callback(&sub_m)?
    }

    if let Some(sub_m) = matches.subcommand_matches("autocorr") {
        spatial::callback(&sub_m)?
    }

    Ok(())
}
