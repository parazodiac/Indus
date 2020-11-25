use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

use indicatif::{ProgressBar, ProgressStyle};
use rand::distributions::{Distribution, Uniform};
use rand::Rng;

use crossbeam::queue::ArrayQueue;
use std::sync::{mpsc, Arc};

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
        region_id: usize,
    ) -> Result<(), Box<dyn Error>> {
        let _norm: u32 = self.stats.iter().sum();
        for (mat_index, val) in self.stats.iter().enumerate() {
            if *val == 0 {
                continue;
            }
            let state = State::new(mat_index, self.num_rows, self.num_cols);

            write!(
                ofile,
                "{}\t{}\t{}\t{}\n",
                mm_obj.get_feature_string(false, sec_feats[state.sec]),
                mm_obj.get_feature_string(true, pivot_feats[state.pivot]),
                val,
                region_id,
                //*val as f32 / norm as f32
            )?;
        }
        Ok(())
    }

    pub fn _stats(&self) -> &Vec<u32> {
        &self.stats
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
    num_samples: usize,
    links_obj: &links::Links<f32>,
    mm_obj: &multimodal::MultiModalExperiment<f32>,
    cells: Option<&Vec<usize>>,
) -> Result<Gamma, Box<dyn Error>> {
    let num_cells = match cells {
        Some(cells) => cells.len(),
        None => mm_obj.num_cells(),
    };

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

    let pivot_mat = mm_obj.get_dense_submatrix(cells, pivot_feats, true);
    let sec_mat = mm_obj.get_dense_submatrix(cells, sec_feats, false);

    let mut stats = vec![0_u32; num_sec_feats * num_pivot_feats];
    for _ in 0..num_samples {
        // sample a cell
        // looking for cell id in the submatrix
        let cell_id_sec = cells_dist.sample(&mut rng);
        // match cells.is_some() {
        //     true => cells.unwrap()[cells_dist.sample(&mut rng)],
        //     false => cells_dist.sample(&mut rng),
        // };

        {
            // sample from sec
            let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);

            let pivot_feat = pivot_feats[state.pivot];
            let sec_hits = links_obj.entry_from_pivot(pivot_feat);
            let sec_indices: Vec<usize> = sec_hits
                .into_iter()
                .map(|x| sec_feats.iter().position(|y| y == x).unwrap())
                .collect();

            state.sec =
                mm_obj.choose_feature(&sec_mat, &sec_indices, coin_toss_value, cell_id_sec)?;
        }

        {
            // sample from the anchors
            let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);
            let pivot_cell = links_obj.jump_cell_id(cell_id_sec, coin_toss_value);
            let pivot_cell_sub_matrix = match links_obj.has_anchors() {
                true => cells.unwrap().iter().position(|&x| x == pivot_cell).unwrap(),
                false => pivot_cell,
            };

            // sample from pivot
            let coin_toss_value: f32 = rng.gen_range(0.0, 1.0);
            let sec_feat = sec_feats[state.sec];
            let pivot_hits = links_obj.entry_to_pivot(sec_feat);
            let pivot_indices: Vec<usize> = pivot_hits
                .into_iter()
                .map(|x| pivot_feats.iter().position(|y| y == x).unwrap())
                .collect();

            state.pivot = mm_obj.choose_feature(
                &pivot_mat,
                &pivot_indices,
                coin_toss_value,
                pivot_cell_sub_matrix,
            )?;
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
    links_obj: &links::Links<f32>,
    regions: &links::IQRegions,
    mut ofile: BufWriter<File>,
    cells: Option<&Vec<usize>>,
) -> Result<(), Box<dyn Error>> {
    let num_regions = regions.len();
    let pbar = ProgressBar::new(num_regions as u64);
    pbar.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
            )
            .progress_chars("╢▌▌░╟"),
    );

    let num_threads = 10;
    let q = Arc::new(ArrayQueue::<usize>::new(num_regions));
    (0..num_regions).for_each(|x| q.push(x).unwrap());

    let (tx, rx) = mpsc::sync_channel(num_threads);
    crossbeam::scope(|scope| {
        for _worker in 0..num_threads {
            let tx = tx.clone();
            let reader = Arc::clone(&q);

            scope.spawn(move |_| loop {
                match reader.pop() {
                    Some(index) => {
                        let pivot_feats = regions.get(index);
                        let sec_feats = links_obj.get_from_pivot_hits(&pivot_feats);
                        let gamma = process_region(
                            &sec_feats,
                            &pivot_feats,
                            crate::configs::NUM_SAMPLES,
                            &links_obj,
                            &mm_obj,
                            cells,
                        )
                        .expect("can't process gamma region");
                        tx.send(Some((gamma, sec_feats, pivot_feats)))
                            .expect("Could not send mid data!");
                    }
                    None => {
                        tx.send(None).expect("Could not send end data!");
                        break;
                    }
                }
            });
        }

        let mut num_regions = 0;
        let mut dead_thread_count = 0;
        for out_data in rx.iter() {
            match out_data {
                Some((gamma, sec_feats, pivot_feats)) => {
                    pbar.inc(1);
                    num_regions += 1;
                    gamma
                        .write(&mut ofile, mm_obj, &sec_feats, &pivot_feats, num_regions)
                        .expect("can't write gamma");
                } // end-Some
                None => {
                    dead_thread_count += 1;
                    if dead_thread_count == num_threads {
                        drop(tx);

                        // consume what's remaining
                        for out_data in rx.iter() {
                            pbar.inc(1);
                            num_regions += 1;
                            out_data.map_or((), |(gamma, sec_feats, pivot_feats)| {
                                gamma
                                    .write(
                                        &mut ofile,
                                        mm_obj,
                                        &sec_feats,
                                        &pivot_feats,
                                        num_regions,
                                    )
                                    .expect("can't write gamma");
                            });
                        }

                        break;
                    } // end if
                } // end-None
            } // end-match
        } // end-for
    })
    .unwrap(); //end crossbeam

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

        let pivot_feats = vec![2];
        let sec_feats = links_obj.get_from_pivot_hits(&pivot_feats);
        assert_eq!(sec_feats, vec![1, 2, 3, 4, 5]);

        let gamma =
            gibbs::process_region(&sec_feats, &pivot_feats, 100_000, &links_obj, &mm_obj, None)
                .unwrap();
        let norm: u32 = gamma._stats().clone().iter().sum();

        let exp_gamma = vec![0.223, 0.211, 0.548, 0.018, 0.000];
        let is_reasonable: bool = gamma
            ._stats()
            .into_iter()
            .enumerate()
            .any(|(index, &x)| ((x as f32 / norm as f32) - exp_gamma[index]).abs() > 1e-2);

        assert!(!is_reasonable);
    }
}
