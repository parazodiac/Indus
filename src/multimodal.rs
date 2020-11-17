use rand::Rng;
use sce::SingleCellExperiment;
use std::error::Error;
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
    pub fn get_experiment(&self, index: usize) -> Option<&SingleCellExperiment<T>> {
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

    pub fn num_cells(&self) -> usize {
        self.get_experiment(0).unwrap().rows()
    }

    pub fn features(&self) -> Vec<&Vec<String>> {
        let mut features = Vec::new();
        for assay in &self.assays {
            features.push(assay.col_names());
        }

        features
    }

    pub fn get_feature_string(&self, is_pivot: bool, index: usize) -> &str {
        let assay = match is_pivot {
            true => &self.assays()[1],
            false => &self.assays()[0],
        };
        &assay.col_names()[index]
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

    pub fn choose_feature(
        &self,
        features: &Vec<usize>,
        coin_val: f32,
        cell_id: usize,
        is_pivot_assay: bool,
    ) -> Result<usize, Box<dyn Error>> {
        if features.len() == 1 {
            return Ok(features[0]);
        }

        assert!(coin_val < 1.0 && coin_val >= 0.0, "wrong coin toss value");
        let mat = match is_pivot_assay {
            true => self.get_experiment(1).unwrap().counts(),
            false => self.get_experiment(0).unwrap().counts(),
        };

        let mut stats: Vec<f32> = features
            .iter()
            .map(|&feature| *mat.get(cell_id, feature).unwrap_or(&0.0))
            .collect();

        let norm: f32 = stats.iter().sum();
        if norm == 0.0 {
            let mut rng = rand::thread_rng();
            return Ok(features[rng.gen_range(0, features.len())]);
        }

        let mut cum_sum_iter = stats.iter_mut().scan(0.0_f32, |cusum, x| {
            *cusum = *cusum + (*x / norm);
            Some(*cusum)
        });

        let chosen_index = cum_sum_iter.position(|x| x > coin_val).unwrap();
        Ok(features[chosen_index])
    }
}

#[cfg(test)]
mod tests {
    use crate::multimodal::MultiModalExperiment;
    use std::path::Path;

    #[test]
    fn test_mmexp() {
        let ppath = Path::new("test/pivot");
        let spath = Path::new("test/sec");
        let mm_obj =
            MultiModalExperiment::from_paths(vec![spath.to_path_buf(), ppath.to_path_buf()]);

        assert_eq!(mm_obj.get_feature_string(true, 0), "FAM138A");
        assert_eq!(mm_obj.num_cells(), 5);
        assert_eq!(mm_obj.len(), 2);
        assert_eq!(mm_obj._pivot(), 0);

        assert_eq!(
            *mm_obj.features()[0],
            vec![
                "chr1-10126-10439",
                "chr1-180794-181148",
                "chr1-181394-181705",
                "chr1-191396-192016",
                "chr1-267971-268188",
                "chr1-280584-280784",
                "chr1-629918-630126",
                "chr1-633995-634215"
            ]
        );
    }
}
