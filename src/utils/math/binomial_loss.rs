use crate::vars::data_cell::*;
use crate::utils::math::beta_function::*;

use statrs::function::beta::*;
// pub fn binomial_loss(cell_data: &CellData,
//                      cluster_centers: &Vec<Vec<f32>>,
//                      log_prior: f32,
//                      cellnum: usize) -> Vec<f32> {
//
//     let mut log_probabilities: Vec<f32> = Vec::new();
//     let mut sum = 0.0;
//
//     for (cluster, center) in cluster_centers.iter().enumerate() {
//         log_probabilities.push(log_prior);
//
//         for (locus_index, locus) in cell_data.loci.iter().enumerate() {
//             log_probabilities[cluster] += cell_data.log_binomial_coefficient[locus_index] +
//                 (cell_data.alt_counts[locus_index] as f32) * center[*locus].ln() +
//                 (cell_data.ref_counts[locus_index] as f32) * (1.0 - center[*locus]).ln();
//         }
//
//         sum += log_probabilities[cluster];
//     }
//
//     log_probabilities
// }

pub fn beta_binomial_loss(cell_data: &CellData,
                          cluster_centers: &Vec<Vec<(f32, f32)>>,
                          log_prior: f32,
                          cellnum: usize) -> Vec<f32> {


    let mut log_probabilities: Vec<f32> = Vec::new();
    let mut sum = 0.0;

    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);

        for (locus_idx, locus) in cell_data.loci.iter().enumerate() {

            log_probabilities[cluster] += cell_data.log_binomial_coefficient[locus_idx] -
                     (ln_beta(center[*locus].0 as f64,
                              center[*locus].1 as f64) as f32) +
                     (ln_beta(center[*locus].0 as f64 + cell_data.alt_counts[locus_idx] as f64,
                              center[*locus].1 as f64 + cell_data.ref_counts[locus_idx] as f64) as f32);


            // log_probabilities[cluster] += cell_data.log_binomial_coefficient[locus_idx] -
            //     log_beta_function(center[*locus].0,
            //                       center[*locus].1) +
            //     log_beta_function(center[*locus].0 + cell_data.alt_counts[locus_idx] as f32,
            //                       center[*locus].1 + cell_data.ref_counts[locus_idx] as f32);

        }

        sum += log_probabilities[cluster];
    }

    log_probabilities
}