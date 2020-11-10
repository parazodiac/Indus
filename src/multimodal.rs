use sce::SingleCellExperiment;
use std::path::PathBuf;

pub struct MultiModalExperiment<T> {
    matrices: Vec<SingleCellExperiment<T>>,
    pivot: usize,
}

impl<T> MultiModalExperiment<T> {
    pub fn get_experiment(&self, index: usize) -> Option<&SingleCellExperiment<T>> {
        self.matrices.get(index)
    }

    pub fn num_modalities(&self) -> usize {
        self.matrices.len()
    }

    pub fn matrices(&self) -> &Vec<SingleCellExperiment<T>> {
        &self.matrices
    }

    pub fn pivot(&self) -> usize {
        self.pivot
    }
}

impl MultiModalExperiment<f32> {
    pub fn from_paths(paths: Vec<PathBuf>) -> MultiModalExperiment<f32> {
        let mut matrices: Vec<sce::SingleCellExperiment<f32>> = Vec::new();
        for path in paths {
            let experiment = sce::SingleCellExperiment::from_alevin(
                path
            ).expect("error reading the input matrix");
    
            info!("{:?}", experiment);
            matrices.push(experiment);
        }

        MultiModalExperiment {
            matrices: matrices,
            pivot: 0,
        }
    }
}