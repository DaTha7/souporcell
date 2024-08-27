#[derive(Debug)]
pub struct CellData {
    pub allele_fractions: Vec<f32>,
    pub log_binomial_coefficient: Vec<f32>,
    pub alt_counts: Vec<u32>,
    pub ref_counts: Vec<u32>,
    pub loci: Vec<usize>,
    pub total_alleles: f32,
}

impl CellData {
    pub(crate) fn new() -> CellData {
        CellData{
            allele_fractions: Vec::new(),
            log_binomial_coefficient: Vec::new(),
            alt_counts: Vec::new(),
            ref_counts: Vec::new(),
            loci: Vec::new(),
            total_alleles: 0.0,
        }
    }
}
