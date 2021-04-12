use crate::config::ProbT;
use crate::fragment::Fragment;
use crate::model;
use crate::quantify;
use crate::record::{AssayRecords, Experiment};

use clap::ArgMatches;
use indicatif::{ProgressBar, ProgressStyle};

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
    common_cells: &Vec<String>,
) -> Result<Vec<HashMap<u64, HashMap<u32, ProbT>>>, Box<dyn Error>> {
    // reading anchors
    let string_index_common_cells: HashMap<String, u32> = common_cells
        .iter()
        .enumerate()
        .map(|(i, x)| (x.clone(), i as u32))
        .collect();

    let mut vec_anchor_triplets = Vec::with_capacity(5);

    let anchor_file_paths = carina::file::files_path_from_clap(sub_m, "anchors")?;
    for file_path in anchor_file_paths {
        let mut anchor_triplets: HashMap<u64, HashMap<u32, ProbT>> = HashMap::with_capacity(10_000);

        let reader = carina::file::bufreader_from_filepath(file_path)?;
        for line in reader.lines() {
            let record = line.unwrap();
            let toks: Vec<&str> = record.split_whitespace().collect();

            let common_cell_index: u32 = *string_index_common_cells
                .get(toks[0])
                .expect("can't find cell name in the common cell list");

            let cb_str: Vec<&str> = toks[1].split("-").collect();
            let cb_id = cb_str.get(1).unwrap().parse::<u8>().unwrap();

            let assay_cell_index: u64 =
                carina::barcode::cb_string_to_u64_with_id(cb_str[0].as_bytes(), 16, cb_id)?;

            let score = toks[2].parse::<ProbT>().unwrap();
            anchor_triplets
                .entry(assay_cell_index)
                .or_insert(HashMap::new())
                .insert(common_cell_index, score);
        }

        vec_anchor_triplets.push(anchor_triplets);
    }

    Ok(vec_anchor_triplets)
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let hmm = model::get_hmm_params(&sub_m)?;
    info!("Read HMM model paramers: {:?}", hmm);

    let common_cells = get_cells(&sub_m)?;
    let num_common_cells = common_cells.len();
    info!(
        "Found {} cells in common assay, first cell={}",
        common_cells.len(),
        common_cells[0]
    );

    let vec_anchor_triplets = get_anchors(&sub_m, &common_cells)?;
    let num_assays = vec_anchor_triplets.len();
    let assay_num_cells: Vec<usize> = vec_anchor_triplets.iter().map(|x| x.len()).collect();
    let assay_num_anchors: Vec<usize> = vec_anchor_triplets
        .iter()
        .map(|x| x.iter().map(|(_, z)| z.len()).sum())
        .collect();
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
    let assay_data: Vec<AssayRecords<ProbT>> = frags
        .iter_mut()
        .enumerate()
        .map(|(i, x)| {
            let cell_records = x.fetch(
                tids[i],
                &range,
                &vec_anchor_triplets.get(i).unwrap(),
                num_common_cells,
            );

            AssayRecords::new(cell_records)
        })
        .collect();
    let exp = Experiment::new(assay_data);

    let num_common_cells = 1;
    info!("Starting forward backward");
    let pbar = ProgressBar::new(num_common_cells as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    for cell_id in 0..num_common_cells {
        pbar.inc(1);
        let quants = quantify::run_fwd_bkw(exp.get_cell_data(cell_id), &hmm)?;
        println!("{:?}", quants);
    }

    pbar.finish();
    info!("All Done");

    Ok(())
}
