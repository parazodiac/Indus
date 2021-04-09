use std::error::Error;
use clap::ArgMatches;
use std::ops::Range;
use std::io::BufRead;
use std::collections::HashMap;
use crate::fragment::Fragment;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let mut common_cells: Vec<String> = Vec::with_capacity(5000);
    { // reading in cell names of common assay.
        let cells_file_path = carina::file::file_path_from_clap(sub_m, "common_cells")?;
        let cells_reader = carina::file::bufreader_from_filepath(cells_file_path)?;
        
        for line in cells_reader.lines() {
            common_cells.push(line.unwrap());
        }
        info!("Found {} cells in common assay, first cell={}", common_cells.len(), common_cells[0]);
    }
    let string_index_common_cells: HashMap<String, u32> = common_cells.iter()
        .enumerate()
        .map(|(i, x)| (x.clone(), i as u32))
        .collect();
    let num_common_cells = common_cells.len();

    let mut vec_anchor_triplets = Vec::with_capacity(5);
    let mut vec_string_index_assay_cells = Vec::with_capacity(5);
    { // reading anchors
        let anchor_file_paths = carina::file::files_path_from_clap(sub_m, "anchors")?;
        for file_path in anchor_file_paths {
            let mut anchor_triplets: Vec<(u32, u32, f32)> = Vec::with_capacity(num_common_cells);
            let mut string_index_assay_cells = HashMap::new();

            let reader = carina::file::bufreader_from_filepath(file_path)?;
            for line in reader.lines() {
                let record = line.unwrap();
                let toks: Vec<&str> = record.split_whitespace().collect();
                
                let common_cell_index: u32 = *string_index_common_cells.get(toks[0])
                    .expect("can't find cell name in the common cell list");
                
                let num_assay_cells = string_index_assay_cells.len();
                let assay_cell_index: u32 = *string_index_assay_cells.entry(toks[1].to_owned())
                    .or_insert(num_assay_cells as u32);

                let score = toks[2].parse::<f32>().unwrap();
                anchor_triplets.push((common_cell_index, assay_cell_index, score));
            }

            vec_anchor_triplets.push(anchor_triplets);
            vec_string_index_assay_cells.push(string_index_assay_cells);
        }

        let num_assays = vec_anchor_triplets.len();
        let assay_num_cells: Vec<usize> = vec_string_index_assay_cells.iter().map(|x| x.len()).collect();
        let assay_num_anchors: Vec<usize> = vec_anchor_triplets.iter().map(|x| x.len()).collect();
        info!("Found {} assays with {:?} cells and {:?} anchors", num_assays, assay_num_cells, assay_num_anchors);
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