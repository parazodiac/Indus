use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::error::Error;

use rand::Rng;
use rand::distributions::{Distribution, Uniform};

use crate::links;
use crate::multimodal;

pub struct Gamma {
    num_rows: usize,
    num_cols: usize,
    stats: Vec<u32>
}

impl Gamma {
    pub fn write(&self, ofile: &mut BufWriter<File>) -> Result<(), Box<dyn Error>> {
        let mut row_index = 0;
        let mut col_index = 0;
        for val in &self.stats {
            col_index += 1;
            if col_index == self.num_cols {
                col_index = 0;
                row_index += 1;
            }

            assert!(row_index < self.num_rows, "wrong increment while writing");
            write!(ofile, "{}\t{}\t{}", row_index, col_index, val)?;
        }
        Ok(())
    }
}

pub struct State {
    sec: usize,
    pivot: usize
}

impl State {
    pub fn row_major_index(&self, num_cols: usize) -> usize {
        (self.sec * num_cols) + self.pivot
    }
}

pub fn process_region(
    sec_feats: &Vec<usize>,
    pivot_feats: &Vec<usize>,
    links_obj: &links::Links<f32>,
    mm_obj: &multimodal::MultiModalExperiment<f32>,
) -> Result<Gamma, Box<dyn Error>> {
    //let sec_matrix = mm_obj.get_dense_submatrix(0, &sec_feats);
    //let pivot_matrix = mm_obj.get_dense_submatrix(1, &pivot_feats);
    let num_cells = mm_obj.num_cells();
    let num_sec_feats = sec_feats.len();
    let num_pivot_feats = pivot_feats.len();
    
    let pivot_dist = Uniform::from(0..num_pivot_feats);
    let sec_dist = Uniform::from(0..num_sec_feats);
    let cells_dist = Uniform::from(0..num_cells);

    let mut rng = rand::thread_rng();
    let mut state = State {
        sec: sec_dist.sample(&mut rng), 
        pivot: pivot_dist.sample(&mut rng)
    };

    let mut stats = vec![0_u32; num_sec_feats * num_pivot_feats];
    for _ in 0..crate::configs::NUM_SAMPLES {
        // sample a cell
        let cell_id = cells_dist.sample(&mut rng);

        // sample from sec
        let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);
        let sec_hits = links_obj.entry_from_pivot(state.pivot);
        state.sec = mm_obj.choose_feature(sec_hits, coin_toss_value, cell_id, false)?;

        // sample from pivot
        let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);
        let pivot_hits = links_obj.entry_to_pivot(state.sec);
        state.pivot = mm_obj.choose_feature(pivot_hits, coin_toss_value, cell_id, true)?;

        stats[state.row_major_index(num_pivot_feats)] += 1;
    }

    Ok(Gamma {
        num_rows: num_sec_feats, 
        num_cols: num_pivot_feats, 
        stats
    })
}

pub fn callback(
    mm_obj: &multimodal::MultiModalExperiment<f32>,
    links_obj: links::Links<f32>,
    regions: links::IQRegions,
    mut ofile: BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    for pivot_feats in regions.groups() {
        let sec_feats = links_obj.get_from_pivot_hits(&pivot_feats);
        let gamma = process_region(&sec_feats, &pivot_feats, &links_obj, &mm_obj)?;
        gamma.write(&mut ofile)?;
    }

    Ok(())
}
