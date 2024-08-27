use crate::vars::data_cell::*;
use crate::vars::config_params::*;
use crate::vars::init_cluster::*;

use crate::utils::cluster::from_data::{known_genotypes::*, known_cells::*};
use crate::utils::cluster::from_strategy::{kmpp::*, midvar::*, randassign::*, uniform::*};

use hashbrown::HashMap;
use rand::rngs::StdRng;

// pub fn init_cluster_centers(loci_used: usize,
//                             cell_data: &Vec<CellData>,
//                             params: &Params,
//                             rng: &mut StdRng,
//                             locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<f32>> {
//
//     if let Some(known_genotypes) = &params.known_genotypes {
//         return init_cluster_centers_known_genotypes(loci_used, params, rng, locus_to_index);
//     }
//     else if let Some(assigned_cells) = &params.known_cell_assignments {
//         return init_cluster_centers_known_cells(loci_used, &cell_data, params, rng);
//     }
//     else {
//         match params.initialization_strategy {
//             ClusterInit::KmeansPP => Vec::new(),
//             ClusterInit::RandomUniform => init_cluster_centers_uniform(loci_used, params, rng),
//             ClusterInit::RandomAssignment => init_cluster_centers_random_assignment(loci_used, &cell_data, params, rng),
//             ClusterInit::MiddleVariance => init_cluster_centers_middle_variance(loci_used, &cell_data, params, rng),
//         }
//     }
//
// }


pub fn init_cluster_centers(loci_used: usize,
                            cell_data: &Vec<CellData>,
                            params: &Params,
                            rng: &mut StdRng,
                            locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<(f32, f32)>> {

    if let Some(known_genotypes) = &params.known_genotypes {
        return Vec::new();
    }
    else if let Some(assigned_cells) = &params.known_cell_assignments {
        return Vec::new();
    }
    else {
        match params.initialization_strategy {
            ClusterInit::KmeansPP => init_cluster_centers_kmeans_pp(loci_used, &cell_data, params, rng),
            ClusterInit::RandomUniform => init_cluster_centers_uniform(loci_used, &cell_data, params, rng),
            ClusterInit::RandomAssignment => Vec::new(),
            ClusterInit::MiddleVariance => Vec::new(),
        }
    }
}