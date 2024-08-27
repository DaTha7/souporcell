pub fn reset_sums_denoms(loci: usize,
                         sums: &mut Vec<Vec<f32>>,
                         denoms: &mut Vec<Vec<f32>>,
                         cluster_centers: &Vec<Vec<f32>>,
                         num_clusters: usize) {

    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index] = 1.0;
            denoms[cluster][index] = 2.0;
        }
    }

}
