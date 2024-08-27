use rand::{SeedableRng, rngs::StdRng};

pub struct ThreadData {
    pub best_log_probabilities: Vec<Vec<f32>>,
    pub best_total_log_probability: f32,
    pub rng: StdRng,
    pub solves_per_thread: usize,
    pub thread_num: usize,
}

impl ThreadData {
    pub fn from_seed(seed: [u8; 32], solves_per_thread: usize, thread_num: usize) -> ThreadData {
        ThreadData {
            best_log_probabilities: Vec::new(),
            best_total_log_probability: f32::NEG_INFINITY,
            rng: SeedableRng::from_seed(seed),
            solves_per_thread,
            thread_num,
        }
    }
}