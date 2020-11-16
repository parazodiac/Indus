use clap::ArgMatches;
use std::error::Error;

use crate::carina;
use crate::links;
use crate::multimodal;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let ipaths = carina::file::files_path_from_clap(sub_m, "ipaths")?;
    assert!(ipaths.len() > 1, "indus expects at least two matrices");

    let mm_obj = multimodal::MultiModalExperiment::from_paths(ipaths);
    println!("{:?}", mm_obj);

    let olap_path = carina::file::file_path_from_clap(sub_m, "links")?;
    let links_obj = links::Links::new(&mm_obj, olap_path);
    println!("{:?}", links_obj);

    Ok(())
}
