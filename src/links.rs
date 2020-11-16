use crate::multimodal;

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::iter::FromIterator;
use std::path::PathBuf;

pub struct IQRegions {
    groups: Vec<Vec<usize>>,
}

impl IQRegions {
    pub fn len(&self) -> usize {
        self.groups.len()
    }
}

impl fmt::Debug for IQRegions {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Total {} groups in the region", self.len())?;

        Ok(())
    }
}

pub struct Links<'a, T> {
    _mm_obj: &'a multimodal::MultiModalExperiment<T>,
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

        let mut maps = Vec::new();
        {
            let all_features = mm_obj.features();
            let mut feature_string_to_index = Vec::<HashMap<String, usize>>::new();
            for features in all_features {
                let mut indexing = HashMap::<String, usize>::new();
                for (index, feature) in features.into_iter().enumerate() {
                    indexing.insert(feature.to_owned(), index);
                }
                feature_string_to_index.push(indexing);
            }

            for line in rdr.records() {
                let record = line.unwrap();
                let values: Vec<String> =
                    record.into_iter().flat_map(str::parse::<String>).collect();
                assert_eq!(values.len(), 2);

                let map = (
                    *feature_string_to_index[0].get(&values[0]).unwrap(),
                    *feature_string_to_index[1].get(&values[1]).unwrap(),
                );
                maps.push(map);
            }
        } // end populating maps

        Links {
            _mm_obj: mm_obj,
            maps,
        }
    }

    pub fn len(&self) -> usize {
        self.maps.len()
    }

    // TODO: Currently it's very ugly looking code, needs more work.
    pub fn extract_iqr(&self) -> Result<IQRegions, Box<dyn Error>> {
        let mut fwd = HashMap::<usize, Vec<usize>>::new();
        let mut rev = HashMap::<usize, Vec<usize>>::new();

        for map in &self.maps {
            let m1 = map.0;
            let m2 = map.1;

            let val = fwd.entry(m1).or_insert(Vec::new());
            val.push(m2);

            let val = rev.entry(m2).or_insert(Vec::new());
            val.push(m1);
        }

        let mut groups = Vec::new();
        let mut pivot_features: HashSet<usize> = rev.keys().into_iter().map(|x| *x).collect();
        while pivot_features.len() > 0 {
            let group = self.extract_region(&fwd, &rev, *pivot_features.iter().nth(0).unwrap());
            let group_set = HashSet::from_iter(group.iter().cloned());

            groups.push(group);
            pivot_features = pivot_features.difference(&group_set).map(|x| *x).collect();
        }

        Ok(IQRegions { groups })
    }

    fn extract_region(
        &self,
        fwd: &HashMap<usize, Vec<usize>>,
        rev: &HashMap<usize, Vec<usize>>,
        seed: usize,
    ) -> Vec<usize> {
        let mut rhits = vec![seed];
        let mut num_rhits = 1;
        let mut num_fhits = 0;

        let mut rchange = true;
        let mut fchange = true;

        let all_hits = |map: &HashMap<usize, Vec<usize>>, query: &Vec<usize>| -> Vec<usize> {
            let mut hits = HashSet::new();
            for elem in query {
                let hit: HashSet<usize> = map.get(elem).unwrap().iter().map(|x| *x).collect();
                hits = hits.union(&hit).into_iter().map(|x| *x).collect();
            }

            let hits: Vec<usize> = hits.into_iter().collect();
            hits
        };

        while rchange || fchange {
            let fhits = all_hits(rev, &rhits);
            rhits = all_hits(fwd, &fhits);

            if rhits.len() != num_rhits {
                num_rhits = rhits.len();
                rchange = true;
            } else {
                rchange = false;
            }

            if fhits.len() != num_fhits {
                num_fhits = fhits.len();
                fchange = true;
            } else {
                fchange = false;
            }
        }

        rhits
    }
}
