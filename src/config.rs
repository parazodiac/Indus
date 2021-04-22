pub type RangeT = u32;
pub type ProbT = f32;

pub const MIN_PROB: f32 = 1e-2;
pub const WINDOW_SIZE: usize = 200;
pub static THRESHOLDS: &[ProbT] = &[0.000, 0.005, 0.005, 0.0023, 0.0038, 0.000];
