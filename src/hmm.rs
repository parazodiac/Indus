use crate::fragment::Fragment;
use crate::model;

use clap::ArgMatches;
use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;
use std::ops::Range;

fn get_cells(sub_m: &ArgMatches) -> Result<Vec<String>, Box<dyn Error>> {
    // reading in cell names of common assay.
    let cells_file_path = carina::file::file_path_from_clap(sub_m, "common_cells")?;
    let cells_reader = carina::file::bufreader_from_filepath(cells_file_path)?;

    let mut common_cells: Vec<String> = Vec::with_capacity(5000);
    for line in cells_reader.lines() {
        common_cells.push(line.unwrap());
    }

    Ok(common_cells)
}

fn get_anchors(
    sub_m: &ArgMatches,
    string_index_common_cells: &HashMap<String, u32>,
) -> Result<(Vec<Vec<(u32, u32, f32)>>, Vec<HashMap<String, u32>>), Box<dyn Error>> {
    // reading anchors
    let mut vec_anchor_triplets = Vec::with_capacity(5);
    let mut vec_string_index_assay_cells = Vec::with_capacity(5);

    let anchor_file_paths = carina::file::files_path_from_clap(sub_m, "anchors")?;
    for file_path in anchor_file_paths {
        let mut anchor_triplets: Vec<(u32, u32, f32)> = Vec::with_capacity(1_000_000);
        let mut string_index_assay_cells = HashMap::with_capacity(10_000);

        let reader = carina::file::bufreader_from_filepath(file_path)?;
        for line in reader.lines() {
            let record = line.unwrap();
            let toks: Vec<&str> = record.split_whitespace().collect();

            let common_cell_index: u32 = *string_index_common_cells
                .get(toks[0])
                .expect("can't find cell name in the common cell list");

            let num_assay_cells = string_index_assay_cells.len();
            let assay_cell_index: u32 = *string_index_assay_cells
                .entry(toks[1].to_owned())
                .or_insert(num_assay_cells as u32);

            let score = toks[2].parse::<f32>().unwrap();
            anchor_triplets.push((common_cell_index, assay_cell_index, score));
        }

        vec_anchor_triplets.push(anchor_triplets);
        vec_string_index_assay_cells.push(string_index_assay_cells);
    }

    Ok((vec_anchor_triplets, vec_string_index_assay_cells))
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let hmm = model::get_hmm_params(&sub_m)?;
    info!("Read HMM paramers: {:?}", hmm);

    let common_cells = get_cells(&sub_m)?;
    info!(
        "Found {} cells in common assay, first cell={}",
        common_cells.len(),
        common_cells[0]
    );

    let string_index_common_cells: HashMap<String, u32> = common_cells
        .iter()
        .enumerate()
        .map(|(i, x)| (x.clone(), i as u32))
        .collect();

    let (vec_anchor_triplets, vec_string_index_assay_cells) =
        get_anchors(&sub_m, &string_index_common_cells)?;
    let num_assays = vec_anchor_triplets.len();
    let assay_num_cells: Vec<usize> = vec_string_index_assay_cells
        .iter()
        .map(|x| x.len())
        .collect();
    let assay_num_anchors: Vec<usize> = vec_anchor_triplets.iter().map(|x| x.len()).collect();
    info!(
        "Found {} assays with {:?} cells and {:?} anchors",
        num_assays, assay_num_cells, assay_num_anchors
    );

    let fragment_file_paths = carina::file::files_path_from_clap(sub_m, "fragments")?;
    let mut frags: Vec<Fragment> = fragment_file_paths
        .into_iter()
        .map(|x| Fragment::from_pathbuf(x))
        .collect();

    let range = Range {
        start: 0,
        end: 250_000_000,
    };
    let tids: Vec<u64> = frags.iter().map(|x| x.tid("chr1")).collect();
    let data: Vec<usize> = frags
        .iter_mut()
        .enumerate()
        .map(|(i, x)| x.fetch(tids[i], &range).len())
        .collect();

    println!("{:?}", data);
    Ok(())
}
