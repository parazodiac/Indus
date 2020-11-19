use crate::multimodal;

use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fmt;
use std::iter::FromIterator;
use std::path::PathBuf;

#[derive(PartialEq)]
pub struct IQRegions {
    groups: Vec<Vec<usize>>,
}

impl IQRegions {
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    pub fn groups(&self) -> &Vec<Vec<usize>> {
        &self.groups
    }

    pub fn _sort(&mut self) {
        self.groups.sort();
    }

    pub fn get(&self, index: usize) -> Vec<usize> {
        self.groups()[index].clone()
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
    to_pivot: HashMap<usize, Vec<usize>>,
    from_pivot: HashMap<usize, Vec<usize>>,
    microclusters: Option<HashMap<String, Vec<usize>>>,
    anchors: Option<HashMap<usize, Vec<usize>>>,
}

impl<'a, T> fmt::Debug for Links<'a, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Found total {} Linked pivot features", self.len())?;
        if let Some(microclusters) = &self.microclusters {
            write!(f, "\nFound total {} micorclusters", microclusters.len())?;
        }

        if let Some(anchors) = &self.anchors {
            write!(f, "\nFound total {} pivot anchors", anchors.len())?;
        }

        Ok(())
    }
}

impl<'a, T> Links<'a, T> {
    pub fn new(mm_obj: &multimodal::MultiModalExperiment<T>, links_file_path: PathBuf) -> Links<T> {
        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(links_file_path)
            .expect("can't read the links file");

        let (to_pivot, from_pivot) = Links::get_links(&mut rdr, &mm_obj);
        Links {
            _mm_obj: mm_obj,
            from_pivot,
            to_pivot,
            microclusters: None,
            anchors: None,
        }
    }

    pub fn new_with_microclusters(
        mm_obj: &multimodal::MultiModalExperiment<T>,
        links_file_path: PathBuf,
        microclusters_file_path: PathBuf,
    ) -> Links<T> {
        let mut links = Links::new(&mm_obj, links_file_path);

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(microclusters_file_path)
            .expect("can't read the microcluster file");

        let clusters = Links::get_microclusters(&mut rdr, &mm_obj);
        links.set_microclusters(clusters);

        links
    }

    pub fn new_with_anchors(
        mm_obj: &multimodal::MultiModalExperiment<T>,
        links_file_path: PathBuf,
        anchors_file_path: PathBuf,
    ) -> Links<T> {
        let mut links = Links::new(&mm_obj, links_file_path);

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .from_path(anchors_file_path)
            .expect("can't read the microcluster file");

        let clusters = Links::get_anchors(&mut rdr, &mm_obj);
        links.set_anchors(clusters);

        links
    }

    fn get_anchors(
        rdr: &mut csv::Reader<std::fs::File>,
        mm_obj: &multimodal::MultiModalExperiment<T>,
    ) -> HashMap<usize, Vec<usize>> {
        let mut anchors = HashMap::<usize, Vec<usize>>::new();
        {
            let all_cells = mm_obj.cells();
            let mut cells_string_to_index = HashMap::<String, usize>::new();
            for (index, cell_string) in all_cells.into_iter().enumerate() {
                cells_string_to_index.insert(cell_string.to_owned(), index);
            }

            for line in rdr.records() {
                let record = line.unwrap();
                let values: Vec<String> =
                    record.into_iter().flat_map(str::parse::<String>).collect();
                assert_eq!(values.len(), 2);

                // todo handle cases when value not in the input matrix
                let sec_cb_index = *cells_string_to_index.get(&values[0]).unwrap();
                let pivot_cb_index = *cells_string_to_index.get(&values[1]).unwrap();

                anchors
                    .entry(pivot_cb_index)
                    .or_insert(Vec::new())
                    .push(sec_cb_index);
            }
        } // end populating maps

        anchors
    }

    fn get_microclusters(
        rdr: &mut csv::Reader<std::fs::File>,
        mm_obj: &multimodal::MultiModalExperiment<T>,
    ) -> HashMap<String, Vec<usize>> {
        let mut clusters = HashMap::<String, Vec<usize>>::new();
        {
            let all_cells = mm_obj.cells();
            let mut cells_string_to_index = HashMap::<String, usize>::new();
            for (index, cell_string) in all_cells.into_iter().enumerate() {
                cells_string_to_index.insert(cell_string.to_owned(), index);
            }

            for line in rdr.records() {
                let record = line.unwrap();
                let values: Vec<String> =
                    record.into_iter().flat_map(str::parse::<String>).collect();
                assert_eq!(values.len(), 2);

                // todo handle cases when value not in the input matrix
                let cb_index = *cells_string_to_index.get(&values[0]).unwrap();
                let cl_id = values[1].clone();

                clusters.entry(cl_id).or_insert(Vec::new()).push(cb_index);
            }
        } // end populating maps

        clusters
    }

    fn get_links(
        rdr: &mut csv::Reader<std::fs::File>,
        mm_obj: &multimodal::MultiModalExperiment<T>,
    ) -> (HashMap<usize, Vec<usize>>, HashMap<usize, Vec<usize>>) {
        let mut to_pivot = HashMap::<usize, Vec<usize>>::new();
        let mut from_pivot = HashMap::<usize, Vec<usize>>::new();
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

                // todo handle cases when value not in the input matrix
                let sec_index = *feature_string_to_index[0].get(&values[0]).unwrap();
                let pivot_index = *feature_string_to_index[1].get(&values[1]).unwrap();

                to_pivot
                    .entry(sec_index)
                    .or_insert(Vec::new())
                    .push(pivot_index);
                from_pivot
                    .entry(pivot_index)
                    .or_insert(Vec::new())
                    .push(sec_index);
            }
        } // end populating maps

        (to_pivot, from_pivot)
    }

    pub fn set_microclusters(&mut self, clusters: HashMap<String, Vec<usize>>) {
        self.microclusters = Some(clusters);
    }

    pub fn has_microclusters(&self) -> bool {
        self.microclusters.is_some()
    }

    pub fn microcluster(&self) -> Option<&HashMap<std::string::String, std::vec::Vec<usize>>> {
        self.microclusters.as_ref()
    }

    pub fn set_anchors(&mut self, anchors: HashMap<usize, Vec<usize>>) {
        self.anchors = Some(anchors);
    }

    pub fn len(&self) -> usize {
        self.from_pivot.len()
    }

    pub fn get_pivot_features(&self) -> HashSet<usize> {
        self.from_pivot.keys().into_iter().map(|x| *x).collect()
    }

    pub fn _get_sec_features(&self) -> HashSet<usize> {
        self.to_pivot.keys().into_iter().map(|x| *x).collect()
    }

    pub fn get_to_pivot_hits(&self, query: &Vec<usize>) -> Vec<usize> {
        let mut all_hits = Vec::new();
        for elem in query {
            let mut elem_hits = self.to_pivot.get(&elem).unwrap().clone();
            all_hits.append(&mut elem_hits);
        }

        all_hits.sort();
        all_hits.dedup();
        all_hits
    }

    pub fn get_from_pivot_hits(&self, query: &Vec<usize>) -> Vec<usize> {
        let mut all_hits = Vec::new();
        for elem in query {
            let mut elem_hits = self.from_pivot.get(&elem).unwrap().clone();
            all_hits.append(&mut elem_hits);
        }

        all_hits.sort();
        all_hits.dedup();
        all_hits
    }

    pub fn entry_to_pivot(&self, query: usize) -> &Vec<usize> {
        self.to_pivot.get(&query).unwrap()
    }

    pub fn entry_from_pivot(&self, query: usize) -> &Vec<usize> {
        self.from_pivot.get(&query).unwrap()
    }

    // TODO: Currently it's very ugly looking code, needs more work.
    pub fn extract_iqr(&self) -> Result<IQRegions, Box<dyn Error>> {
        let mut groups = Vec::new();
        let mut pivot_features = self.get_pivot_features();
        while pivot_features.len() > 0 {
            let elem = *pivot_features.iter().nth(0).unwrap();
            let group = self.extract_region(elem);
            let group_set = HashSet::from_iter(group.iter().cloned());

            groups.push(group);
            pivot_features = pivot_features.difference(&group_set).map(|x| *x).collect();
        }

        Ok(IQRegions { groups })
    }

    fn extract_region(&self, seed: usize) -> Vec<usize> {
        let mut pivot_hits = vec![seed];
        let mut num_pivot_hits = 1;
        let mut num_sec_hits = 0;

        let mut pivot_change = true;
        let mut sec_change = true;

        while pivot_change || sec_change {
            let sec_hits = self.get_from_pivot_hits(&pivot_hits);
            pivot_hits = self.get_to_pivot_hits(&sec_hits);

            if pivot_hits.len() != num_pivot_hits {
                num_pivot_hits = pivot_hits.len();
                pivot_change = true;
            } else {
                pivot_change = false;
            }

            if sec_hits.len() != num_sec_hits {
                num_sec_hits = sec_hits.len();
                sec_change = true;
            } else {
                sec_change = false;
            }
        }

        pivot_hits
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::iter::FromIterator;
    use std::path::Path;

    use crate::links::IQRegions;
    use crate::links::Links;
    use crate::multimodal::MultiModalExperiment;

    #[test]
    fn test_links() {
        let ppath = Path::new("test/pivot");
        let spath = Path::new("test/sec");
        let mm_obj =
            MultiModalExperiment::from_paths(vec![spath.to_path_buf(), ppath.to_path_buf()]);

        let opath = Path::new("test/olaps.tsv");
        let links_obj = Links::new(&mm_obj, opath.to_path_buf());

        assert_eq!(links_obj.len(), 4);
        assert_eq!(
            links_obj.get_pivot_features(),
            HashSet::from_iter(vec![0, 1, 2, 3])
        );
        assert_eq!(
            links_obj._get_sec_features(),
            HashSet::from_iter(vec![6, 2, 1, 7, 3, 5, 4, 0])
        );
        assert_eq!(links_obj.entry_to_pivot(7), &vec![3, 1]);
        assert_eq!(links_obj.entry_from_pivot(1), &vec![0, 7]);
        assert_eq!(links_obj.get_to_pivot_hits(&vec![0, 3]), vec![0, 1, 2]);
        assert_eq!(links_obj.get_from_pivot_hits(&vec![0, 3]), vec![0, 6, 7]);
        assert_eq!(links_obj.extract_region(0), vec![0, 1, 3]);
        assert_eq!(
            links_obj.extract_iqr().unwrap()._sort(),
            IQRegions {
                groups: vec![vec![0, 1, 3], vec![2]]
            }
            ._sort()
        );
    }
}
