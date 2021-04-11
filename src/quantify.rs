use sprs::CsMat;

use crate::config::ProbT;
use crate::record::CellRecords;

use std::error::Error;

pub fn get_posterior(cell_records: &CellRecords<ProbT>) -> Result<CsMat<ProbT>, Box<dyn Error>> {
    info!("{:?}", cell_records.len());

    Ok(CsMat::eye(3))
}
