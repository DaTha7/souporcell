use crate::vars::data_cell::*;

pub fn update_centers_average(sums: &mut Vec<Vec<f32>>,
                          denoms: &mut Vec<Vec<f32>>,
                          cell: &CellData,
                          probabilities: &Vec<f32>) {

    for locus in 0..cell.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            sums[cluster][cell.loci[locus]] += probabilities[cluster] * (cell.alt_counts[locus] as f32);
            denoms[cluster][cell.loci[locus]] += probabilities[cluster] * ((cell.alt_counts[locus] + cell.ref_counts[locus]) as f32);
        }
    }
}

pub fn update_final(loci: usize, sums: &Vec<Vec<f32>>, denoms: &Vec<Vec<f32>>, cluster_centers: &mut Vec<Vec<f32>>) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            let update = sums[cluster][locus]/denoms[cluster][locus];
            cluster_centers[cluster][locus] = update.min(0.99).max(0.01);//max(0.0001, min(0.9999, update));
        }
    }
}

pub fn update_beta_variables(cellnum: usize,
                             cell: &CellData,
                             probabilities: &Vec<f32>,
                             cluster_centers: &mut Vec<Vec<(f32, f32)>>) {

    for (cluster, probability) in probabilities.iter().enumerate() {
        for (loci_idx, loci) in cell.loci.iter().enumerate() {
            cluster_centers[cluster][*loci].0 += probability * cell.alt_counts[loci_idx] as f32;
            cluster_centers[cluster][*loci].1 += probability * cell.ref_counts[loci_idx] as f32;
        }
    }
}