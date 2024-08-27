use rand::Rng;
use rand::rngs::StdRng;

use rand::prelude::*;
use rand::distributions::WeightedIndex;

use crate::vars::data_cell::*;
use crate::vars::config_params::*;
use crate::utils::math::log_distance::*;
use crate::utils::sys::search::*;

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
        cell_log_likelihoods.push(0.0); // fake value to make the search work???

        for cell in cell_data.iter() {
            let mut min_distance = f32::INFINITY;

            for cluster_index in 0..(cluster+1) {
                let distance = get_log_distance(cell, &centers[cluster_index]);
                min_distance = min_distance.min(distance);
            }

            cell_log_likelihoods.push(min_distance);
        }
        
        cell_log_likelihoods.iter_mut().fold(0.0, |log_val, cummulative_sum| {*cummulative_sum += log_val; *cummulative_sum});
        
        let min = cell_log_likelihoods.iter().fold(f32::INFINITY, |a, &b| a.min(b));
        let max = cell_log_likelihoods.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        
        cell_log_likelihoods.iter_mut().for_each(|x| *x = (*x - min)/(max - min));
        
        loop {
            let sampling_value: f32 = rng.gen_range(0.0, 1.0);
            let sampling_index = relative_index_search(&cell_log_likelihoods, sampling_value).unwrap() - 1;
       
            if list_of_chosen_cells.iter().find(|&&x| x == sampling_index) == None {
                chosen_cell_idx = sampling_index;
                chosen_cell = &cell_data[chosen_cell_idx];
                break
            }
        }

        list_of_chosen_cells.push(chosen_cell_idx);

        for (locus_idx, locus) in chosen_cell.loci.iter().enumerate() {
            centers[cluster][*locus].0 += chosen_cell.alt_counts[locus_idx] as f32;
            centers[cluster][*locus].1 += chosen_cell.ref_counts[locus_idx] as f32;
        }
        
        /*
        for cell in cell_data.iter() {
            let distance = get_log_distance(cell, &centers[cluster]);
            cell_log_likelihoods.push(distance);
        }

        // This turns the % to 88.6, but it doesn't make sense
        // cell_log_likelihoods.iter_mut().fold(0.0, |log_val, cummulative_sum| {*cummulative_sum += log_val; *cummulative_sum});

        let min = cell_log_likelihoods.iter().fold(f32::INFINITY, |a, &b| a.min(b));
        let max = cell_log_likelihoods.iter().fold(f32::NEG_INFINITY, |a, &b| a.max(b));
        cell_log_likelihoods.iter_mut().for_each(|x| *x = (*x - min)/(max - min));

        let choices: Vec<usize> = (0..cell_data.len()).collect();
        let dist = WeightedIndex::new(&cell_log_likelihoods).unwrap();

        loop {
            let sampling_index = choices[dist.sample(rng)];
            if list_of_chosen_cells.iter().find(|&&x| x == sampling_index) == None {
                chosen_cell_idx = sampling_index;
                chosen_cell = &cell_data[chosen_cell_idx];
                break
            }
        }
        */
    }

    centers

}
