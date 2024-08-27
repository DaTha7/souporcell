use crate::vars::config_params::*;

use rand::{rngs::StdRng, Rng};
use crate::vars::data_cell::CellData;

// pub fn init_cluster_centers_uniform(loci: usize,
//                                     params: &Params,
//                                     rng: &mut StdRng) -> Vec<Vec<f32>> {
//
//     let mut centers: Vec<Vec<f32>> = Vec::new();
//
//     for cluster in 0..params.num_clusters {
//         centers.push(Vec::new());
//         for _ in 0..loci {
//             centers[cluster].push(rng.gen::<f32>().min(0.9999).max(0.0001));
//         }
//     }
//
//     centers
// }

pub fn init_cluster_centers_uniform(loci: usize,
                                    cell_data: &Vec<CellData>,
                                    params: &Params,
                                    rng: &mut StdRng) -> Vec<Vec<(f32, f32)>> {

    let mut centers:Vec<Vec<(f32,f32)>> = Vec::new();

    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push((1.0, 1.0));
        }
    }

    for cluster in 0..params.num_clusters {
        let chosen_cell_idx= rng.gen_range(0, cell_data.len());
        let chosen_cell = &cell_data[chosen_cell_idx];

        for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
            centers[cluster][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
            centers[cluster][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
        }
    }

    centers
}