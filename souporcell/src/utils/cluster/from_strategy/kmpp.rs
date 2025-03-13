use rand::Rng;
use rand::rngs::StdRng;

use rand::prelude::*;
use rand::distributions::WeightedIndex;

use crate::vars::data_cell::*;
use crate::vars::config_params::*;
use crate::utils::math::log_distance::*;
use crate::utils::sys::search::*;

// pub fn init_cluster_centers_kmeans_pp(loci: usize,
//                                       cell_data: &Vec<CellData>,
//                                       params: &Params,
//                                       rng: &mut StdRng) -> Vec<Vec<(f32, f32)>> {
//
//     // Fetch all loci that have data for >10% of cells
//     let mut locus_population: Vec<u32> = Vec::new();
//     for _ in 0..loci {
//         locus_population.push(0);
//     }
//
//     for cell in cell_data.iter() {
//         for locus_idx in cell.loci.iter() {
//             locus_population[*locus_idx] += 1;
//         }
//     }
//
//     let mut new_population = locus_population.clone();
//
//     for (idx, population) in locus_population.iter().enumerate() {
//         let ratio = (*population as f32) / (cell_data.len() as f32);
//         if ratio < 0.05 {
//             new_population[idx] = 0;
//         }
//     }
//     // ********************************************
//
//
//     let mut centers:Vec<Vec<(f32, f32)>> = Vec::new();
//
//     for cluster in 0..params.num_clusters {
//         centers.push(Vec::new());
//         for _ in 0..loci {
//             centers[cluster].push((1.0, 1.0));
//         }
//     }
//
//     let mut chosen_cell_idx = rng.gen_range(0, cell_data.len());
//     let mut chosen_cell = &cell_data[chosen_cell_idx];
//
//     let mut list_of_chosen_cells: Vec<usize> = Vec::new();
//     list_of_chosen_cells.push(chosen_cell_idx);
//
//     for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
//         centers[0][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
//         centers[0][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
//     }
//
//     for cluster in 1..params.num_clusters {
//
//         let mut cell_log_likelihoods: Vec<f32> = Vec::new();
//         cell_log_likelihoods.push(0.0); // fake value to make the search work???
//
//         for cell in cell_data.iter() {
//             let mut min_distance = f32::INFINITY;
//
//             for cluster_index in 0..(cluster+1) {
//                 let distance = get_log_distance(cell, &new_population, &centers[cluster_index]);
//                 min_distance = min_distance.min(distance);
//             }
//
//             cell_log_likelihoods.push(min_distance);
//         }
//
//
//         // Debug
//         let saved = cell_log_likelihoods.clone();
//         let mut sum = 0.0;
//         let mut denom = 0.0;
//
//         for x in &saved {
//             sum += x;
//             denom += 1.0;
//         }
//
//         let mean = sum/denom;
//         // **********
//
//
//         cell_log_likelihoods.iter_mut().fold(0.0, |log_val, cumulative_sum| {*cumulative_sum += log_val; *cumulative_sum});
//
//         let min = cell_log_likelihoods.iter().fold(f32::INFINITY, |a, &b| a.min(b));
//         let max = cell_log_likelihoods.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
//
//         cell_log_likelihoods.iter_mut().for_each(|x| *x = (*x - min)/(max - min));
//
//         let mut recorded_sampling_value;
//         loop {
//             let sampling_value: f32 = rng.gen_range(0.0, 1.0);
//             recorded_sampling_value = sampling_value;
//             let sampling_index = relative_index_search(&cell_log_likelihoods, sampling_value).unwrap();
//
//             if list_of_chosen_cells.iter().find(|&&x| x == sampling_index - 1) == None {
//                 chosen_cell_idx = sampling_index;
//                 chosen_cell = &cell_data[chosen_cell_idx - 1];
//                 break
//             }
//         }
//
//
//         // Debug
//         let mut cumulative_tmp = Vec::new();
//         let mut distance_tmp = Vec::new();
//         for i in (chosen_cell_idx-3)..(chosen_cell_idx+3) {
//             cumulative_tmp.push(cell_log_likelihoods[i]);
//             distance_tmp.push(saved[i]);
//         }
//         println!("chose {} with distance metric {} vs mean distance metric {}", chosen_cell_idx, saved[chosen_cell_idx], mean);
//         println!("indices {:?}, sampling value {}, normalized cumulative values {:?}, distance metric {:?}", ((chosen_cell_idx-3)..(chosen_cell_idx+3)), recorded_sampling_value, cumulative_tmp, distance_tmp);
//         for cluster_index in 0..(cluster+1) {
//             println!("chosen cell {} with distance to cluster {} of {}", chosen_cell_idx, cluster_index, get_log_distance(chosen_cell, &new_population, &centers[cluster_index]));
//         }
//         // **********
//
//
//         list_of_chosen_cells.push(chosen_cell_idx);
//
//         for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
//             centers[cluster][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
//             centers[cluster][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
//         }
//
//     }
//
//     centers
//
// }


pub fn init_cluster_centers_kmeans_pp(loci: usize,
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

    let mut chosen_cell_idx = rng.gen_range(0, cell_data.len());
    let mut chosen_cell = &cell_data[chosen_cell_idx];

    let mut list_of_chosen_cells: Vec<usize> = Vec::new();
    list_of_chosen_cells.push(chosen_cell_idx);

    for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
        centers[0][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
        centers[0][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
    }

    for cluster in 1..params.num_clusters {

        let mut cell_log_likelihoods: Vec<f32> = Vec::new();

        for cell in cell_data.iter() {
            let mut min_distance = f32::INFINITY;

            for cluster_index in 0..(cluster+1) {
                let distance = get_log_distance(cell, &centers[cluster_index]);
                min_distance = min_distance.min(distance);
            }

            cell_log_likelihoods.push(min_distance);
        }

        let mut weights: Vec<f32> = Vec::new();
        let squared_sum: f32 = cell_log_likelihoods.iter().sum();
        for cell_log_likelihood in cell_log_likelihoods.iter() {
            weights.push(*cell_log_likelihood/squared_sum)
        }

        let choices: Vec<usize> = (0..cell_data.len()).collect();
        let weighted_distances = WeightedIndex::new(&weights).unwrap();

        loop {
            let sampling_value = weighted_distances.sample(rng);
            let sampling_index = choices[sampling_value];

            if list_of_chosen_cells.iter().find(|&&x| x == sampling_index) == None {
                chosen_cell_idx = sampling_index;
                chosen_cell = &cell_data[chosen_cell_idx - 1];
                break
            }
        }


        list_of_chosen_cells.push(chosen_cell_idx);

        for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
            centers[cluster][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
            centers[cluster][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
        }
    }
    centers
}
