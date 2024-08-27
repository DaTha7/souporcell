use rand::{rngs::StdRng, Rng};

pub fn new_seed(rng: &mut StdRng) -> [u8; 32] {
    let mut seed = [0; 32];
    for i in 0..32 {
        seed[i] = rng.gen::<u8>();
    }
    seed
}