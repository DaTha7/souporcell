use crate::vars::data_cell::*;
use crate::vars::config_params::*;

use rand::rngs::StdRng;

pub fn init_cluster_centers_middle_variance(loci: usize,
                                            cell_data: &Vec<CellData>,
                                            params: &Params,
                                            rng: &mut StdRng) -> Vec<Vec<f32>> {
    assert!(false, "middle variance not yet implemented");
    Vec::new()
}