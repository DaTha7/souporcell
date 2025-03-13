use crate::vars::data_cell::CellData;

use statrs::function::beta::*;

pub fn get_log_distance(cell: &CellData,
                        // locus_population: &Vec<u32>,
                        cluster_center: &Vec<(f32,f32)>) -> f32 {

    let mut cell_log_likelihood: f32 = 0.0;
    let mut cell_implied_likelihood: f32 = 0.0;

    for (locus_idx, locus) in cell.loci.iter().enumerate() {

        // if locus_population[locus_idx] == 0 {
        //     continue
        // }

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

    //eprintln!("really wtf is happening");
    //return(-cell_log_likelihood/(cell.loci.len() as f32));
    //((cell_implied_likelihood - cell_log_likelihood).powf(2.0)/(cell.loci.len() as f32))
    //(cell_implied_likelihood - cell_log_likelihood).exp()
    (cell_implied_likelihood - cell_log_likelihood) // best so far 88%, still not as good as 96% for binomial
    //((cell_implied_likelihood - cell_log_likelihood)).exp()
    //1.0
}


// Mahalanobis distance

// pub fn get_log_distance(cell: &CellData,
//                         cluster_center: &Vec<(f32,f32)>) -> f32 {
//
//     let mut locus_mean: Vec<f32> = Vec::new();
//     let mut locus_likelihood: Vec<f32> = Vec::new();
//     let mut covariance: Vec<Vec<f32>> = Vec::new();
//
//     for locus_idx in 0..cell.loci.len() {
//         covariance.push(Vec::new());
//         for _ in 0..cell.loci.len() {
//             covariance[locus_idx].push(0.0);
//         }
//     }
//
//     for (locus_idx, locus) in cell.loci.iter().enumerate() {
//
//         locus_mean.push(((cell.alt_counts[locus_idx] + cell.ref_counts[locus_idx]) as f32 * cluster_center[*locus].0) /
//                         (cluster_center[*locus].0 + cluster_center[*locus].1));
//
//         locus_likelihood.push((cell.log_binomial_coefficient[locus_idx] -
//                         (ln_beta(cluster_center[*locus].0 as f64,
//                                  cluster_center[*locus].1 as f64) as f32) +
//                         (ln_beta(cluster_center[*locus].0 as f64 + cell.alt_counts[locus_idx] as f64,
//                                  cluster_center[*locus].1 as f64 + cell.ref_counts[locus_idx] as f64) as f32)).exp());
//
//         covariance[locus_idx][locus_idx] = (
//             (cluster_center[*locus].0 + cluster_center[*locus].1).powf(2.0) *
//             (cluster_center[*locus].0 + cluster_center[*locus].1 + 1.0)
//         ) /
//             (
//                 ((cell.alt_counts[locus_idx] + cell.ref_counts[locus_idx]) as f32) *
//                     cluster_center[*locus].0 *
//                     cluster_center[*locus].1 *
//                     (
//                         (cell.alt_counts[locus_idx] + cell.ref_counts[locus_idx]) as f32 +
//                             cluster_center[*locus].0 +
//                             cluster_center[*locus].1
//                     )
//             );
//     }
//
//     let mut mat_A: Vec<f32> = Vec::new();
//
//     for (a, b) in locus_likelihood.iter().zip(&locus_mean) {
//         mat_A.push(a - b);
//     }
//
//     let mut mat_B: Vec<f32> = Vec::new();
//     for locus_idx in 0..cell.loci.len() {
//         let cov_results: f32 = mat_A.iter().zip(covariance[locus_idx].iter()).map(|(x, y)| x * y).sum();
//         mat_B.push(cov_results);
//     }
//
//     let distance: f32 = mat_B.iter().zip(mat_A.iter()).map(|(x, y)| x * y).sum();
//
//     distance
// }