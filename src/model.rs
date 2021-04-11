use crate::config::THRESHOLDS;
use clap::ArgMatches;

use std::error::Error;
use std::fmt;
use std::io::BufRead;

use crate::config::ProbT;
pub struct Hmm {
    init: Vec<ProbT>,
    emission: Vec<Vec<ProbT>>,
    transition: Vec<Vec<ProbT>>,
}

impl fmt::Debug for Hmm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Found ")
            .field("#states", &self.num_states())
            .field("#assays", &self.num_assays())
            //.field("initialization proabilities", &self.transition)
            .finish()
    }
}

impl Hmm {
    pub fn get_init_prob(&self, state: usize) -> ProbT {
        self.init[state]
    }

    pub fn get_transition_prob(&self, pstate: usize, state: usize) -> ProbT {
        self.transition[pstate][state]
    }

    pub fn get_emission_prob(&self, state: usize, observations: &Vec<ProbT>) -> ProbT {
        let mut prob = 0.0;
        for assay_id in 0..observations.len() {
            let emission_prob = self.emission[state][assay_id];
            match observations[assay_id] >= THRESHOLDS[assay_id] {
                true => prob *= emission_prob,
                false => prob *= 1.0 - emission_prob,
            }
        }

        prob
    }

    pub fn num_states(&self) -> usize {
        self.init.len()
    }

    pub fn num_assays(&self) -> usize {
        self.emission.get(0).expect("emission not filled").len()
    }

    pub fn new(mut reader: std::io::BufReader<std::fs::File>) -> Hmm {
        let mut first_line = String::new();
        reader
            .read_line(&mut first_line)
            .expect("Unable to read line");

        let toks: Vec<&str> = first_line.split_whitespace().collect();
        let num_states = toks[0].parse::<usize>().unwrap();
        let num_assays = toks[1].parse::<usize>().unwrap();

        let mut init = vec![0.0; num_states];
        let mut emission = vec![vec![0.0; num_assays]; num_states];
        let mut transition = vec![vec![0.0; num_states]; num_states];

        let (mut pcounter, mut tcounter, mut ecounter) = (0, 0, 0);
        for line in reader.lines() {
            let record = line.unwrap();
            let toks: Vec<&str> = record.split_whitespace().collect();

            match toks[0] {
                "probinit" => {
                    assert!(toks.len() == 3);
                    let state = toks[1].parse::<usize>().unwrap() - 1;
                    let probability = toks[2].parse::<ProbT>().unwrap();
                    init[state] = probability + 1.15358E-31;
                    pcounter += 1;
                }
                "transitionprobs" => {
                    assert!(toks.len() == 4);
                    let fstate = toks[1].parse::<usize>().unwrap() - 1;
                    let sstate = toks[2].parse::<usize>().unwrap() - 1;
                    let probability = toks[3].parse::<ProbT>().unwrap();
                    transition[fstate][sstate] = probability;
                    tcounter += 1;
                }
                "emissionprobs" => {
                    assert!(toks.len() == 6);
                    let is_presence = toks[4].parse::<u8>().unwrap();
                    if is_presence != 1 {
                        continue;
                    }

                    let state = toks[1].parse::<usize>().unwrap() - 1;
                    let assay = toks[2].parse::<usize>().unwrap();
                    let probability = toks[5].parse::<ProbT>().unwrap();
                    emission[state][assay] = probability;
                    ecounter += 1;
                }
                _ => unreachable!(),
            }
        } // end-for
        assert_eq!(pcounter, num_states);
        assert_eq!(tcounter, num_states * num_states);
        assert_eq!(ecounter, num_states * num_assays);

        Hmm {
            init,
            emission,
            transition,
        }
    }
}

pub fn get_hmm_params(sub_m: &ArgMatches) -> Result<Hmm, Box<dyn Error>> {
    let hmm_file_path = carina::file::file_path_from_clap(sub_m, "hmm")?;
    let file_reader = carina::file::bufreader_from_filepath(hmm_file_path)?;
    let hmm = Hmm::new(file_reader);

    Ok(hmm)
}
