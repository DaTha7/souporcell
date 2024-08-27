use crate::vars::data_cell::*;
use crate::vars::config_params::*;

use rand::{rngs::StdRng, Rng};

pub fn init_cluster_centers_random_assignment(loci: usize,
                                              cell_data: &Vec<CellData>,
                                              params: &Params,
                                              rng: &mut StdRng) -> Vec<Vec<f32>> {
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.num_clusters {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _ in 0..loci {
            sums[cluster].push(rng.gen::<f32>()*0.01);
            denoms[cluster].push(0.01);
        }
    }
    for cell in cell_data {
        let cluster = rng.gen_range(0,params.num_clusters);
        for locus in 0..cell.loci.len() {
            let alt_c = cell.alt_counts[locus] as f32;
            let total = alt_c + (cell.ref_counts[locus] as f32);
            let locus_index = cell.loci[locus];
            sums[cluster][locus_index] += alt_c;
            denoms[cluster][locus_index] += total;
        }
    }
    for cluster in 0..params.num_clusters {
        for locus in 0..loci {
            sums[cluster][locus] = sums[cluster][locus]/denoms[cluster][locus] + (rng.gen::<f32>()/2.0 - 0.25);
            sums[cluster][locus] = sums[cluster][locus].min(0.9999).max(0.0001);
        }
    }
    let centers = sums;
    centers
}