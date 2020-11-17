use clap::ArgMatches;
use std::error::Error;

use crate::carina;
use crate::gibbs;
use crate::links;
use crate::multimodal;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let ipaths = carina::file::files_path_from_clap(sub_m, "ipaths")?;
    assert!(ipaths.len() > 1, "indus expects at least two matrices");

    info!("Reading quant matrices");
    let mm_obj = multimodal::MultiModalExperiment::from_paths(ipaths);
    info!("{:?}", mm_obj);

    info!("Reading Overlap file");
    let olap_path = carina::file::file_path_from_clap(sub_m, "links")?;
    let links_obj = links::Links::new(&mm_obj, olap_path);
    info!("{:?}", links_obj);

    info!("Finding Independantly quantifiable regions");
    let regions = links_obj.extract_iqr()?;
    info!("Found total {:?} regions", regions.len());

    info!("Starting gibbs sampling");
    let ofile = carina::file::bufwriter_from_clap(sub_m, "output")?;
    gibbs::callback(&mm_obj, links_obj, regions, ofile)?;

    info!("All done");
    Ok(())
}