use crate::records::Records;
use std::ops::Range;

use rust_htslib::tbx::{self, Read};
use std::path::PathBuf;

pub struct Fragment {
    _filepath: PathBuf,
    reader: tbx::Reader,
}

impl Fragment {
    pub fn from_pathbuf(filepath: PathBuf) -> Fragment {
        let tbx_reader =
            tbx::Reader::from_path(&filepath).expect(&format!("Could not open {:?}", filepath));

        Fragment {
            _filepath: filepath,
            reader: tbx_reader,
        }
    }

    pub fn _filepath(&self) -> &str {
        self._filepath.to_str().unwrap()
    }

    pub fn tid(&self, seqname: &str) -> u64 {
        match self.reader.tid(seqname) {
            Ok(tid) => tid,
            Err(_) => panic!("Could not resolve to contig ID"),
        }
    }

    pub fn fetch(&mut self, tid: u64, region: &Range<u32>) -> Vec<Records<u32>> {
        // Set region to fetch.
        self.reader
            .fetch(tid, region.start as u64, region.end as u64)
            .expect("Could not seek to fetch region");

        let data: Vec<Records<u32>> = self
            .reader
            .records()
            .map(|x| String::from_utf8(x.unwrap()).expect("UTF8 conversion error"))
            .map(|x| Records::from_string(x))
            .collect();

        return data;
    }
}
