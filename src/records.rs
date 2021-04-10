use std::collections::HashSet;
use std::ops::Range;

#[derive(Debug)]
pub struct Records<T> {
    range: Range<T>,
    cb: u64,
}

impl Records<u32> {
    pub fn from_string(data: String, assay_cells: &HashSet<u64>) -> Option<Records<u32>> {
        let (mut start, mut end) = (0, 0);
        let mut cb: &str = "chr0";
        let mut id: u8 = 0;
        for (index, text) in data.split_whitespace().enumerate() {
            match index {
                1 => start = text.parse::<u32>().unwrap(),
                2 => end = text.parse::<u32>().unwrap(),
                3 => {
                    let toks: Vec<&str> = text.split("-").collect();
                    cb = toks[0];
                    id = toks
                        .get(1)
                        .expect(&format!("{:?}", toks.get(0)))
                        .parse::<u8>()
                        .unwrap();
                }
                4 => break,
                _ => (),
            }
        }

        let cb = match id {
            1 => carina::barcode::cb_string_to_u64_with_id(cb.as_bytes(), 16, 1).unwrap(),
            2 => carina::barcode::cb_string_to_u64_with_id(cb.as_bytes(), 16, 2).unwrap(),
            _ => unreachable!(),
        };

        match assay_cells.contains(&cb) {
            true => {
                let range = Range { start, end };
                Some(Records { range, cb })
            }
            false => None,
        }
    }
}
