use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use indicatif::{ProgressBar, ProgressStyle};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;

use crate::links;
use crate::multimodal;

#[derive(Debug)]
pub struct Gamma {
    num_rows: usize,
    num_cols: usize,
    stats: Vec<u32>,
}

impl Gamma {
    pub fn write(
        &self,
        ofile: &mut BufWriter<File>,
        mm_obj: &multimodal::MultiModalExperiment<f32>,
        sec_feats: &Vec<usize>,
        pivot_feats: &Vec<usize>,
    ) -> Result<(), Box<dyn Error>> {
        for (mat_index, val) in self.stats.iter().enumerate() {
            if *val == 0 {
                continue;
            }
            let state = State::new(mat_index, self.num_rows, self.num_cols);

            write!(
                ofile,
                "{}\t{}\t{}\n",
                mm_obj.get_feature_string(false, sec_feats[state.sec]),
                mm_obj.get_feature_string(true, pivot_feats[state.pivot]),
                val
            )?;
        }
        Ok(())
    }
}

#[derive(PartialEq, Debug)]
pub struct State {
    sec: usize,
    pivot: usize,
}

impl State {
    pub fn new(index: usize, num_secs: usize, num_pivots: usize) -> State {
        if num_secs == 1 {
            return State {
                sec: 0,
                pivot: index,
            };
        }

        if num_pivots == 1 {
            return State {
                sec: index,
                pivot: 0,
            };
        }

        let pindex = index % num_pivots;
        let sindex = index / num_pivots;
        State {
            sec: sindex,
            pivot: pindex,
        }
    }

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
        pivot: pivot_dist.sample(&mut rng),
    };

    let pivot_mat = mm_obj.get_dense_submatrix(None, pivot_feats, true);
    let sec_mat = mm_obj.get_dense_submatrix(None, sec_feats, false);

    let mut stats = vec![0_u32; num_sec_feats * num_pivot_feats];
    for _ in 0..crate::configs::NUM_SAMPLES {
        // sample a cell
        let cell_id = cells_dist.sample(&mut rng);

        {
            // sample from sec
            let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);

            let pivot_feat = pivot_feats[state.pivot];
            let sec_hits = links_obj.entry_from_pivot(pivot_feat);
            let sec_indices: Vec<usize> = sec_hits.into_iter()
                .map(|x| sec_feats.iter().position(|y| y == x).unwrap())
                .collect();

            state.sec = mm_obj.choose_feature(&sec_mat, &sec_indices, coin_toss_value, cell_id, false)?;
        }

        {
            // sample from pivot
            let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);

            let sec_feat = sec_feats[state.sec];
            let pivot_hits = links_obj.entry_to_pivot(sec_feat);
            let pivot_indices: Vec<usize> = pivot_hits.into_iter()
                .map(|x| pivot_feats.iter().position(|y| y == x).unwrap())
                .collect();

            state.pivot = mm_obj.choose_feature(&pivot_mat, &pivot_indices, coin_toss_value, cell_id, true)?;
        }

        stats[state.row_major_index(num_pivot_feats)] += 1;
    }

    Ok(Gamma {
        num_rows: num_sec_feats,
        num_cols: num_pivot_feats,
        stats,
    })
}

pub fn callback(
    mm_obj: &multimodal::MultiModalExperiment<f32>,
    links_obj: links::Links<f32>,
    regions: links::IQRegions,
    mut ofile: BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    let num_regions = regions.len();
    let pbar = ProgressBar::new(num_regions as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    for pivot_feats in regions.groups() {
        // progress bar increment
        pbar.inc(1);

        let sec_feats = links_obj.get_from_pivot_hits(&pivot_feats);
        let gamma = process_region(&sec_feats, &pivot_feats, &links_obj, &mm_obj)?;
        gamma.write(&mut ofile, mm_obj, &sec_feats, pivot_feats)?;
    }

    pbar.finish();
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::Path;

    use crate::gibbs;
    use crate::links::Links;
    use crate::multimodal::MultiModalExperiment;

    #[test]
    fn test_state() {
        let state = gibbs::State::new(45, 9, 20);
        assert_eq!(state, gibbs::State { sec: 2, pivot: 5 });
    }

    #[test]
    fn test_gibbs() {
        let ppath = Path::new("test/pivot");
        let spath = Path::new("test/sec");
        let mm_obj =
            MultiModalExperiment::from_paths(vec![spath.to_path_buf(), ppath.to_path_buf()]);

        let opath = Path::new("test/olaps.tsv");
        let links_obj = Links::new(&mm_obj, opath.to_path_buf());

        let pivot_feats = vec![0, 1, 3];
        let sec_feats = links_obj.get_from_pivot_hits(&pivot_feats);
        assert_eq!(sec_feats, vec![0, 6, 7]);

        let gamma = gibbs::process_region(&sec_feats, &pivot_feats, &links_obj, &mm_obj).unwrap();
        println!("gamma: {:?}", gamma);

        assert_eq!(0, 1);
    }
}
