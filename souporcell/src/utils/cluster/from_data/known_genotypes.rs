use crate::vars::config_params::*;
use crate::utils::sys::read::*;

use rand::rngs::StdRng;
use hashbrown::HashMap;
use vcf::*;

pub fn init_cluster_centers_known_genotypes(loci: usize,
                                            params: &Params,
                                            rng: &mut StdRng,
                                            locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<f32>> {

    let mut centers: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(0.5);
        }
    }
    let mut vcf_reader = VCFReader::new(reader(params.known_genotypes.as_ref().unwrap())).unwrap();
    let mut locus_id: usize = 0;
    for record in vcf_reader {
        let record = record.unwrap();
        if let Some(loci_index) = locus_to_index.get(&locus_id) {
            if params.known_genotypes_sample_names.len() > 0 {
                for (sample_index, sample) in params.known_genotypes_sample_names.iter().enumerate() {
                    let gt = record.call[sample]["GT"][0].to_string();
                    // complicated way of getting the haplotype to numbers
                    let hap0 = gt.chars().nth(0).unwrap().to_string();
                    if hap0 == "." { continue; }
                    let hap0 = hap0.parse::<u32>().unwrap().min(1);
                    let hap1 = gt.chars().nth(2).unwrap().to_string().parse::<u32>().unwrap().min(1);
                    centers[sample_index][*loci_index] = (((hap0 + hap1) as f32)/2.0).min(0.99).max(0.01);
                }
            } else { assert!(false, "currently requiring known_genotypes_sample_names if known_genotypes set"); }
        }
        locus_id += 1;
    }
    centers
}