use crate::multimodal;
use std::fmt;
use std::path::PathBuf;

pub struct Links<'a, T> {
    mm_obj: &'a multimodal::MultiModalExperiment<T>,
    maps: Vec<(usize, usize)>,
}

impl<'a, T> fmt::Debug for Links<'a, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Found total {} Links", self.len())?;

        Ok(())
    }
}

impl<'a, T> Links<'a, T> {
    pub fn new(mm_obj: &multimodal::MultiModalExperiment<T>, file_path: PathBuf) -> Links<T> {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(file_path)
            .expect("can't read the overlap file");

        let features = mm_obj.features();
        let get_index = |names: &Vec<String>, query: &str| -> usize {
            names.iter().position(|r| r == query).unwrap()
        };

        let mut maps = Vec::new();
        for line in rdr.records() {
            let record = line.unwrap();
            let values: Vec<String> = record.into_iter().flat_map(str::parse::<String>).collect();
            assert_eq!(values.len(), 2);

            let map = (
                get_index(features[0], &values[0]),
                get_index(features[1], &values[1]),
            );
            maps.push(map);
        }

        Links { mm_obj, maps }
    }

    pub fn len(&self) -> usize {
        self.maps.len()
    }
}
