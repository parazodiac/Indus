use std::error::Error;
use std::io::Read;
use std::io::Write;
use std::convert::TryInto;
use std::sync::{mpsc, Arc};

use clap::ArgMatches;
use crossbeam::queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};

use crate::config::CHR_LENS;
use crate::hmm;

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let common_cells = hmm::get_cells(&sub_m)?;
    let num_common_cells = common_cells.len();
    info!(
        "Found {} cells in common assay, first cell={}",
        num_common_cells, common_cells[0]
    );

    let num_states = 12;
    let num_chrs = CHR_LENS.len();
    info!("Found total {} chromosomes", num_chrs);

    let in_path = carina::file::file_path_from_clap(&sub_m, "in_directory").unwrap();
    info!("Found input directory path: {:?}", in_path);

    let out_path = carina::file::file_path_from_clap(&sub_m, "out_directory").unwrap();
    info!("Found output directory path: {:?}", out_path);

    info!("Starting to read");
    (0..num_chrs).rev().for_each(|chr_id| {
        let chr_name = format!("chr{}", chr_id+1);
        info!("Working on {}", chr_name);

        let pbar = ProgressBar::new(num_common_cells as u64);
        pbar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
                )
                .progress_chars("╢▌▌░╟"),
        );

        let q = Arc::new(ArrayQueue::<usize>::new(num_common_cells));
        //(0..num_common_cells).filter(|&x| x == 2840).for_each(|x| q.push(x).unwrap());
        (0..num_common_cells).for_each(|x| q.push(x).unwrap());

        let num_threads = 5;
        let (tx, rx) = mpsc::sync_channel(num_threads);

        let chr_path = in_path.join(&chr_name);
        std::fs::create_dir_all(&chr_path).unwrap();

        let arc_in_path = Arc::new(&chr_path);
        let arc_common_cells = Arc::new(&common_cells);

        crossbeam::scope(|scope| {
            for _ in 0..num_threads {
                let tx = tx.clone();
                let reader = Arc::clone(&q);
                let arc_in_path = Arc::clone(&arc_in_path);
                let arc_common_cells = Arc::clone(&arc_common_cells);

                scope.spawn(move |_| loop {
                    match reader.pop() {
                        Some(cell_id) => {
                            let cell_file = arc_in_path.join(format!("{}.bin", arc_common_cells[cell_id]));
                            let mut file_handle = carina::file::bufreader_from_filepath(cell_file).unwrap();

                            let mut mat_u8_sizes = vec![0 as u8; 12];
                            file_handle.read_exact(&mut mat_u8_sizes).expect("can't read sizes");
                            let nnz: u32 = u32::from_le_bytes(mat_u8_sizes[0..4].try_into().unwrap());
                            let _nrows: u32 = u32::from_le_bytes(mat_u8_sizes[4..8].try_into().unwrap());
                            let _ncols: u32 = u32::from_le_bytes(mat_u8_sizes[8..12].try_into().unwrap());

                            let mut indices = vec![0 as u8; nnz as usize * 4];
                            let mut states = vec![0 as u8; nnz as usize];
                            let mut probs = vec![0 as u8; nnz as usize];

                            file_handle.read_exact(&mut probs).expect("can't read probs");
                            file_handle.read_exact(&mut states).expect("can't read states");
                            file_handle.read_exact(&mut indices).expect("can't read indices");

                            let mut state_indices = vec![Vec::<u8>::new(); num_states];
                            let mut state_probs = vec![Vec::new(); num_states];
                            for (idx, state) in states.into_iter().enumerate() {
                                let state: usize = u8::from_le(state) as usize;
                                state_probs[state].push(probs[idx]);
                                state_indices[state].extend(&indices[idx*4..(idx+1)*4]);
                            }

                            tx.send(Some((state_indices, state_probs, cell_id)))
                                .expect("Could not send mid data!");
                        }
                        None => {
                            tx.send(None).expect("Could not send end data!");
                            break;
                        }
                    }
                });
            }

            let chr_path = out_path.join(&chr_name);
            std::fs::create_dir_all(&chr_path).unwrap();

            let mut file_handles: Vec<std::io::BufWriter<std::fs::File>> = (0..num_states).map(|x| {
                let file_path = chr_path.join(&format!("{}.bin", x+1));
                std::io::BufWriter::new(std::fs::File::create(file_path).unwrap())
            }).collect();
            let mut cell_id_handle = std::io::BufWriter::new(std::fs::File::create(chr_path.join(&"cells.txt")).unwrap());

            let mut dead_thread_count = 0;
            for out_data in rx.iter() {
                match out_data {
                    Some((state_indices, state_probs, cell_id)) => {
                        for i in 0..num_states {
                            file_handles[i].write_all(&state_probs[i]).unwrap();
                            file_handles[i].write_all(&state_indices[i]).unwrap();
                        }
                        writeln!(cell_id_handle, "{}", common_cells[cell_id]).unwrap();
                        pbar.inc(1);
                    } // end-Some
                    None => {
                        dead_thread_count += 1;
                        if dead_thread_count == num_threads {
                            drop(tx);
                            break;
                        }
                    } // end-None
                } // end-match
            } // end-for
        })
        .unwrap(); //end crossbeam
        pbar.finish();
    }); // end for loop over chromosomes

    info!("All Done");
    Ok(())
}