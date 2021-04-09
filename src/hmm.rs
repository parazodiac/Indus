use std::error::Error;
use clap::ArgMatches;
use std::ops::Range;
use std::io::BufRead;

use crate::fragment::Fragment;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let mut cells: Vec<String> = Vec::with_capacity(5000);
    { // reading in cell names of common assay.
        let cells_file_path = carina::file::file_path_from_clap(sub_m, "common_cells")?;
        let cells_reader = carina::file::bufreader_from_filepath(cells_file_path)?;
        
        for line in cells_reader.lines() {
            cells.push(line.unwrap());
        }
        info!("Found {} cells in common assay, first cell={}", cells.len(), cells[0]);
    }

    let fragment_file_paths = carina::file::files_path_from_clap(sub_m, "fragments")?;
    let mut frags: Vec<Fragment> = fragment_file_paths.into_iter()
        .map(|x| Fragment::from_pathbuf(x))
        .collect();

    let range = Range {start: 0, end: 250_000_000};
    let tids: Vec<u64> = frags.iter().map(|x| x.tid("chr1")).collect();
    let data: Vec<usize> = frags.iter_mut().enumerate().map(|(i,x)| x.fetch(tids[i], &range).len()).collect();

    println!("{:?}", data);
    Ok(())
}