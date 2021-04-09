use std::ops::Range;

#[derive(Debug)]
pub struct Records<T> {
    range: Range<T>,
    cb: u64
}

impl Records<u32> {
    pub fn from_string(data: String) -> Records<u32> {
        let (mut start, mut end) = (0, 0);
        let mut cb: &str = "chr0";
        for (index, text) in data.split_whitespace().enumerate() {
            match index {
                1 => start = text.parse::<u32>().unwrap(),
                2 => end = text.parse::<u32>().unwrap(),
                3 => cb = text.split("-").nth(0).unwrap(),
                4 => break,
                _ => (),
            }
        }

        let range = Range{start, end};
        let cb = carina::barcode::cb_string_to_u64(cb.as_bytes()).unwrap();
        Records {
            range, cb
        }
    }
}