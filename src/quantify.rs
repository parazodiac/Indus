use bio::data_structures::interval_tree::IntervalTree;
use sprs::CsMat;

use crate::config::{ProbT, WINDOW_SIZE};
use crate::model::Hmm;
use crate::record::CellRecords;

use std::error::Error;

fn forward(
    observations: &Vec<Vec<ProbT>>,
    hmm: &Hmm,
    num_states: usize,
    num_observations: usize,
) -> Result<(ProbT, Vec<Vec<ProbT>>), Box<dyn Error>> {
    let mut f_curr = vec![0.0; num_states];
    let mut f_prev = vec![0.0; num_states];
    let mut fprob = vec![vec![0.0; num_states]; num_observations];

    for i in 0..num_observations {
        for state in 0..num_states {
            let mut prev_f_sum = 0.0;
            if i == 0 {
                prev_f_sum = hmm.get_init_prob(state);
            } else {
                for prev_state in 0..num_states {
                    prev_f_sum += f_prev[prev_state] * hmm.get_transition_prob(prev_state, state);
                }
            }

            f_curr[state] = hmm.get_emission_prob(state, &observations[i]) * prev_f_sum;
        }
        let prob_norm: ProbT = f_curr.iter().sum();
        f_curr.iter_mut().for_each(|x| *x /= prob_norm);

        fprob[i] = f_curr.clone();
        f_prev = f_curr.clone();
    }

    let mut norm = f_curr.into_iter().map(|x| x * 0.1).sum();
    if norm == 0.0 {
        norm = 1.0;
    }

    Ok((norm, fprob))
}

fn backward(
    observations: Vec<Vec<ProbT>>,
    hmm: &Hmm,
    norm: ProbT,
    num_states: usize,
    num_observations: usize,
    fprob: Vec<Vec<ProbT>>,
    posterior: &mut sprs::TriMat<ProbT>,
) -> Result<(), Box<dyn Error>> {
    let mut b_curr = vec![0.1; num_states];
    let mut b_prev = vec![0.1; num_states];

    let update_triplet = |i: usize, posterior: &mut sprs::TriMat<ProbT>, b_curr: &Vec<ProbT>| {
        let probs: Vec<ProbT> = (0..num_states)
            .map(|state| fprob[i][state] * b_curr[state] / norm)
            .collect();
        let state_norm: ProbT = probs.iter().sum();
        probs.into_iter().enumerate().for_each(|(state, prob)| {
            let prob = prob / state_norm;
            if prob > 1e-4 {
                posterior.add_triplet(i, state, prob)
            }
        });
    };
    update_triplet(num_observations - 1, posterior, &b_curr);

    for i in (1..num_observations).rev() {
        let obv_emissions: Vec<ProbT> = (0..num_states)
            .into_iter()
            .map(|state| hmm.get_emission_prob(state, &observations[i]))
            .collect();

        b_curr.iter_mut().for_each(|i| *i = 0.0);
        for state in 0..num_states {
            for next_state in 0..num_states {
                b_curr[state] += hmm.get_transition_prob(state, next_state)
                    * obv_emissions[next_state]
                    * b_prev[next_state];
            }
        }

        let prob_norm: ProbT = b_curr.iter().sum();
        b_curr.iter_mut().for_each(|x| *x /= prob_norm);

        b_prev = b_curr.clone();
        update_triplet(i - 1, posterior, &b_curr);
    }

    Ok(())
}

fn get_posterior(observations: Vec<Vec<ProbT>>, hmm: &Hmm) -> Result<CsMat<ProbT>, Box<dyn Error>> {
    let num_states = hmm.num_states();
    let _num_assays = hmm.num_assays();
    let num_observations = observations.len();

    //assert!(num_assays == observations[0].len());
    let (norm, fprob) = forward(&observations, &hmm, num_states, num_observations)?;
    let mut posterior = sprs::TriMat::new((num_observations, num_states));
    backward(
        observations,
        &hmm,
        norm,
        num_states,
        num_observations,
        fprob,
        &mut posterior,
    )?;

    Ok(posterior.to_csr())
}

pub fn run_fwd_bkw(
    cell_records: Vec<&CellRecords<ProbT>>,
    hmm: &Hmm,
) -> Result<CsMat<ProbT>, Box<dyn Error>> {
    let itrees: Vec<IntervalTree<u32, ProbT>> = cell_records
        .into_iter()
        .map(|cell_records| {
            let mut tree = IntervalTree::new();
            for record in cell_records.records() {
                tree.insert(record.range(), record.id());
            }
            tree
        })
        .collect();

    let observation_list: Vec<Vec<ProbT>> = (0..250_000_000)
        .step_by(WINDOW_SIZE)
        .map(|qstart| {
            let qrange = qstart as u32..(qstart + WINDOW_SIZE) as u32;
            let cts: Vec<ProbT> = itrees
                .iter()
                .map(|tree| {
                    let vals: Vec<ProbT> = tree.find(&qrange).map(|x| *x.data()).collect();
                    vals.iter().sum()
                })
                .collect();
            cts
        })
        .collect();

    let posterior = get_posterior(observation_list, hmm)?;
    Ok(posterior)
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_fwd_bkw() {
        let path =
            std::path::PathBuf::from("/home/srivastavaa/parazodiac/Indus/data/model_test.txt");
        let file_reader = carina::file::bufreader_from_filepath(path).unwrap();
        let hmm = crate::model::Hmm::new(file_reader);

        let mat =
            crate::quantify::get_posterior(vec![vec![0.0], vec![1.0], vec![2.0]], &hmm).unwrap();
        let probs: Vec<String> = mat
            .data()
            .into_iter()
            .map(|&x| format!("{:.4}", x))
            .collect();
        assert_eq!(
            probs,
            ["0.5616", "0.4384", "0.8942", "0.1058", "0.9050", "0.0950"]
        );
    }
}
