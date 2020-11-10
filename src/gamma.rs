use clap::ArgMatches;
use std::error::Error;

use crate::carina;
use crate::multimodal;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let ipaths = carina::file::files_path_from_clap(sub_m, "ipaths")?;
    assert!(ipaths.len() > 1, "indus expects at least two matrices");

    let mm_obj = multimodal::MultiModalExperiment::from_paths(ipaths);
    info!("Found total {} modalities", mm_obj.num_modalities());

    Ok(())
}
