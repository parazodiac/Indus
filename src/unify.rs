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

    info!("Creating Link object");
    let olap_path = carina::file::file_path_from_clap(sub_m, "links")?;
    let links_obj = match carina::file::try_file_path_from_clap(sub_m, "microclusters") {
        Some(mpath) => match carina::file::try_file_path_from_clap(sub_m, "anchors") {
            Some(apath) => links::Links::new_with_microclusters_and_anchors(&mm_obj, olap_path, mpath, apath),
            None => links::Links::new_with_microclusters(&mm_obj, olap_path, mpath),
        },
        None => match carina::file::try_file_path_from_clap(sub_m, "anchors") {
            Some(apath) => links::Links::new_with_anchors(&mm_obj, olap_path, apath),
            None => links::Links::new(&mm_obj, olap_path),
        },
    };
    info!("{:?}", links_obj);

    info!("Finding Independantly quantifiable regions");
    let regions = links_obj.extract_iqr()?;
    info!("Found total {:?} regions", regions.len());

    info!("Starting gibbs sampling");
    match links_obj.has_microclusters() {
        false => {
            let ofile = carina::file::bufwriter_from_clap(sub_m, "output")?;
            gibbs::callback(&mm_obj, &links_obj, &regions, ofile, None)?;
        }
        true => {
            for (key, value) in links_obj.microcluster().unwrap() {
                info!("Working on microcluster {}", key);
                let ofile = carina::file::bufwriter_from_clap_with_suffix(sub_m, "output", key)?;
                gibbs::callback(&mm_obj, &links_obj, &regions, ofile, Some(value))?;
            }
        }
    }

    info!("All done");
    Ok(())
}
