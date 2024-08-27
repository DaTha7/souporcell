use std::fs::File;
use std::io::{BufWriter, Write};

use crate::vars::data_cell::*;
use crate::vars::config_params::*;

use crate::utils::math::reset_sums_denoms::*;
use crate::utils::math::binomial_loss::*;
use crate::utils::math::log_sum_exp::*;
use crate::utils::math::normalize_in_log_wt_temp::*;

use crate::utils::cluster::update::*;


// pub fn EM (loci: usize,
//            mut cluster_centers: Vec<Vec<f32>>,
//            cell_data: &Vec<CellData>,
//            params: &Params,
//            epoch: usize,
//            thread_num: usize) -> (f32, Vec<Vec<f32>>) {
//
//     let mut sums: Vec<Vec<f32>> = Vec::new();
//     let mut denoms: Vec<Vec<f32>> = Vec::new();
//
//     for cluster in 0..params.num_clusters {
//         sums.push(Vec::new());
//         denoms.push(Vec::new());
//         for index in 0..loci {
//             sums[cluster].push(1.0);
//             denoms[cluster].push(2.0); // psuedocounts
//         }
//     }
//
//     let log_prior: f32 = (1.0/(params.num_clusters as f32)).ln();
//     let mut change = 1000.0;
//     let mut iterations = 0;
//
//     //let mut cell_probabilities: Vec<Vec<f32>> = Vec::new();
//     //for _cell in cell_data {
//     //    cell_probabilities.push(Vec::new());
//     //}
//
//     let mut total_log_loss = f32::NEG_INFINITY;
//     let mut total_log_loss_binom = f32::NEG_INFINITY;
//     let mut final_log_probabilities = Vec::new();
//
//     for _cell in 0..cell_data.len() {
//         final_log_probabilities.push(Vec::new());
//     }
//
//     let log_loss_change_limit = 0.01*(cell_data.len() as f32);
//     let temp_steps = 9;
//     let mut last_log_loss = f32::NEG_INFINITY;
//
//
//     for temp_step in 0..temp_steps {
//         let mut log_loss_change = 100000.0;
//
//         while (log_loss_change > log_loss_change_limit && iterations < 1000) {
//
//             let mut log_binom_loss = 0.0;
//
//             reset_sums_denoms(loci, &mut sums, &mut denoms, &cluster_centers, params.num_clusters);
//
//             for (celldex, cell) in cell_data.iter().enumerate() {
//
//                 let log_binoms = binomial_loss(cell, &cluster_centers, log_prior, celldex);
//
//                 log_binom_loss += log_sum_exp(&log_binoms);
//
//                 let mut temp = (cell.total_alleles/(20.0 * 2.0f32.powf((temp_step as f32)))).max(1.0);
//
//                 if temp_step == temp_steps - 1 {
//                     temp = 1.0;
//                 }
//
//                 let probabilities = normalize_in_log_with_temp(&log_binoms, temp);
//
//                 update_centers_average(&mut sums, &mut denoms, cell, &probabilities);
//                 final_log_probabilities[celldex] = log_binoms;
//             }
//
//             total_log_loss = log_binom_loss;
//             log_loss_change = log_binom_loss - last_log_loss;//log_loss - last_log_loss;
//             last_log_loss = log_binom_loss;//log_loss;
//
//             update_final(loci, &sums, &denoms, &mut cluster_centers);
//             iterations += 1;
//             eprintln!("binomial\t thread:{} \tepoch:{} \titer:{} \ttemp_step:{} \tlog_binom_loss:{} \tchange:{}", thread_num, epoch, iterations, temp_step, log_binom_loss, log_loss_change);//, cluster_cells_weighted);
//         }
//     }
//     //for (celldex, probabilities) in cell_probabilities.iter().enumerate() {
//     //    println!("cell {} with {} loci, cluster probabilities {:?}", celldex, cell_data[celldex].loci.len(), probabilities);
//     //}
//     //for center in 0..cluster_centers.len() {
//     //    for locus in 0..cluster_centers[0].len() {
//     //        println!("cluster {} locus {} {}", center, locus, cluster_centers[center][locus]);
//     //    }
//     //}
//     //println!("total log probability = {}",total_log_loss);
//
//     (total_log_loss, final_log_probabilities)
// }

pub fn EM (loci: usize,
           mut cluster_centers: Vec<Vec<(f32, f32)>>,
           cell_data: &Vec<CellData>,
           params: &Params,
           epoch: usize,
           thread_num: usize) -> (f32, Vec<Vec<f32>>) {

    let log_prior: f32 = 0.0;
    let mut iterations = 0;

    let mut total_log_loss = f32::NEG_INFINITY;
    let mut final_log_probabilities = Vec::new();

    for _cell in 0..cell_data.len() {
        final_log_probabilities.push(Vec::new());
    }

    let log_loss_change_limit = 0.01*(cell_data.len() as f32);
    let temp_steps = 9;
    let mut last_log_loss = f32::NEG_INFINITY;

    let temp_step = 1;

    let mut log_loss_change = 100000.0;

    while (log_loss_change > log_loss_change_limit && iterations < 1000) {
        let mut log_binom_loss = 0.0;

        let mut updated_cluster_centers = cluster_centers.clone();

        for (celldex, cell) in cell_data.iter().enumerate() {
            let log_binoms = beta_binomial_loss(cell, &cluster_centers, log_prior, celldex);

            log_binom_loss += log_sum_exp(&log_binoms);

            let mut temp = (cell.total_alleles / (20.0 * 2.0f32.powf((temp_step as f32)))).max(1.0);

            if temp_step == temp_steps - 1 {
                temp = 1.0;
            }

            let probabilities = normalize_in_log_with_temp(&log_binoms, temp);

            update_beta_variables(celldex, cell, &probabilities, &mut updated_cluster_centers);
            final_log_probabilities[celldex] = log_binoms;
        }

        total_log_loss = log_binom_loss;
        log_loss_change = log_binom_loss - last_log_loss;
        last_log_loss = log_binom_loss;

        cluster_centers = updated_cluster_centers.clone();
        iterations += 1;

        eprintln!("binomial \tthread:{} \tepoch:{} \titer:{} \ttemp_step:{} \tlog_binom_loss:{} \t\tchange:{}", thread_num, epoch, iterations, temp_step, log_binom_loss, log_loss_change); //, cluster_cells_weighted);
    }

    (total_log_loss, final_log_probabilities)
}