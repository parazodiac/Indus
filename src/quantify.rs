use sprs::CsMat;
use bio::data_structures::interval_tree::IntervalTree;

use crate::config::ProbT;
use crate::record::CellRecords;

use std::error::Error;

pub fn get_posterior(cell_records: Vec<&CellRecords<ProbT>>) -> Result<CsMat<ProbT>, Box<dyn Error>> {
    let itrees: Vec<IntervalTree<u32, f32>> = cell_records.into_iter().map(|cell_records| {
        let mut tree = IntervalTree::new();
        for record in cell_records.records() {
            tree.insert(record.range(), record.id());
        }
        tree
    }).collect();

    let qrange = 1..1_000_000;
    let cts: Vec<Vec<f32>> = itrees.iter().map(|tree| {
        let vals: Vec<f32> = tree.find(&qrange).map(|x| *x.data()).collect();
        vals
    }).collect();
    println!("{:?}", cts);

    Ok(CsMat::eye(3))
}
