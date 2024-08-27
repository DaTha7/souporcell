use crate::vars::data_cell::CellData;
use statrs::function::beta::*;

pub fn get_log_distance(cell: &CellData,
                        cluster_center: &Vec<(f32,f32)>) -> f32 {

    let mut cell_log_likelihood: f32 = 0.0;
    let mut cell_implied_likelihood: f32 = 0.0;

    for (locus_idx, locus) in cell.loci.iter().enumerate() {

        cell_log_likelihood += cell.log_binomial_coefficient[locus_idx] -
            (ln_beta(cluster_center[*locus].0 as f64,
                     cluster_center[*locus].1 as f64) as f32) +
            (ln_beta(cluster_center[*locus].0 as f64 + cell.alt_counts[locus_idx] as f64,
                     cluster_center[*locus].1 as f64 + cell.ref_counts[locus_idx] as f64) as f32);

        cell_implied_likelihood += cell.log_binomial_coefficient[locus_idx] -
            (ln_beta(1.0 + cell.alt_counts[locus_idx] as f64,
                     1.0 + cell.ref_counts[locus_idx] as f64) as f32) +
            (ln_beta(1.0 + 2.0 * cell.alt_counts[locus_idx] as f64,
                     1.0 + 2.0 * cell.ref_counts[locus_idx] as f64) as f32);
    }

    cell_implied_likelihood - cell_log_likelihood
}