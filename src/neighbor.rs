use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::PathBuf;

use clap::ArgMatches;
use granne::Builder;
use crossbeam::queue::ArrayQueue;
use indicatif::{ProgressBar, ProgressStyle};
use std::sync::{mpsc, Arc};


fn parse_line(line: &str) -> std::io::Result<granne::angular::Vector<'static>> {
    let line_iter = line.split(",");
    let data: granne::angular::Vector = line_iter.map(|d| d.parse::<f32>().unwrap()).collect();

    Ok(data)
}

pub fn generate_neighbors(
    matrix_file_path: PathBuf,
    ofile: BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    let file = BufReader::new(File::open(matrix_file_path)?);

    // reading the input data
    let mut num_rows = 0;
    let mut elements = granne::angular::Vectors::new();
    for line in file.lines() {
        let vector = parse_line(&line?)?;
        elements.push(&vector);
        num_rows += 1;
    }

    println!("Building the index");
    let build_config = granne::BuildConfig::default()
        .show_progress(true);

    let mut builder = granne::GranneBuilder::new(build_config, elements);
    builder.build();

    println!("Querying");
    let index = builder.get_index();
    {
        let pbar = ProgressBar::new(num_rows as u64);
        pbar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "{spinner:.red} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos:>7}/{len:7} {msg}",
                )
                .progress_chars("╢▌▌░╟"),
        );
  
        let num_threads = 10;
        let q = Arc::new(ArrayQueue::<usize>::new(num_rows));
        let arc_index = Arc::new(index);
        (0..num_rows).for_each(|x| q.push(x).unwrap());
    
        let (tx, rx) = mpsc::sync_channel(num_threads);
        crossbeam::scope(|scope| {
            for _worker in 0..num_threads {
                let tx = tx.clone();
                let reader = Arc::clone(&q);
                let index = Arc::clone(&arc_index);
    
                scope.spawn(move |_| loop {
                    match reader.pop() {
                        Some(i) => {
                            let res = index.search(&index.get_element(i), 200, 20);
                            let res: Vec<_> = res.into_iter().map(|(idx, _dist)| idx.to_string()).collect();

                            tx.send(Some(res))
                                .expect("Could not send mid data!");
                        }
                        None => {
                            tx.send(None).expect("Could not send end data!");
                            break;
                        }
                    }
                });
            }

            let mut wtr = csv::WriterBuilder::new()
                .has_headers(false)
                .delimiter(b'\t')
                .from_writer(ofile);
    
            let mut dead_thread_count = 0;
            for out_data in rx.iter() {
                match out_data {
                    Some(res) => {
                        pbar.inc(1);
                        wtr.write_record(&res).unwrap();
                    } // end-Some
                    None => {
                        dead_thread_count += 1;
                        if dead_thread_count == num_threads {
                            drop(tx);
    
                            // consume what's remaining
                            for out_data in rx.iter() {
                                pbar.inc(1);
                                out_data.map_or((), |res| {
                                    wtr.write_record(&res).unwrap();
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
    }

    Ok(())
}

pub fn callback(sub_m: &ArgMatches) -> Result<(), Box<dyn Error>> {
    let matrix_file_path = carina::file::file_path_from_clap(sub_m, "matrix")?;
    let ofile = carina::file::bufwriter_from_clap(sub_m, "output")?;

    generate_neighbors(matrix_file_path, ofile)?;

    info!("All done");
    Ok(())
}
