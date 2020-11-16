use sce::SingleCellExperiment;
use std::fmt;
use std::path::PathBuf;

pub struct MultiModalExperiment<T> {
    assays: Vec<SingleCellExperiment<T>>,
    _pivot: usize,
}

impl<T> fmt::Debug for MultiModalExperiment<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MultiModalExperiment: {} modalitles\n", self.len())?;
        for (index, exp) in self.assays().into_iter().enumerate() {
            write!(f, "Modality {} Shape: {:?}\n", index, exp.shape())?;
        }

        Ok(())
    }
}

impl<T> MultiModalExperiment<T> {
    pub fn _get_experiment(&self, index: usize) -> Option<&SingleCellExperiment<T>> {
        self.assays.get(index)
    }

    pub fn len(&self) -> usize {
        self.assays.len()
    }

    pub fn assays(&self) -> &Vec<SingleCellExperiment<T>> {
        &self.assays
    }

    pub fn _pivot(&self) -> usize {
        self._pivot
    }

    pub fn features(&self) -> Vec<&Vec<String>> {
        let mut features = Vec::new();
        for assay in &self.assays {
            features.push(assay.col_names());
        }

        features
    }
}

impl MultiModalExperiment<f32> {
    pub fn from_paths(paths: Vec<PathBuf>) -> MultiModalExperiment<f32> {
        let mut assays: Vec<sce::SingleCellExperiment<f32>> = Vec::new();
        for path in paths {
            let experiment = sce::SingleCellExperiment::from_tenx_v2(path)
                .expect("error reading the input matrix");

            info!("{:?}", experiment);
            assays.push(experiment);
        }

        MultiModalExperiment {
            assays: assays,
            _pivot: 0,
        }
    }
}
