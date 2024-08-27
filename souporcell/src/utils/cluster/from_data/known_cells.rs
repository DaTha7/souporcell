use crate::vars::data_cell::*;
use crate::vars::config_params::*;

use rand::rngs::StdRng;

pub fn init_cluster_centers_known_cells(loci: usize,
                                        cell_data: &Vec<CellData>,
                                        params: &Params,
                                        rng: &mut StdRng) -> Vec<Vec<f32>> {
    assert!(false, "known cell assignments not yet implemented");
    Vec::new()
}