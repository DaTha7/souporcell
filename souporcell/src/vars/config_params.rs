use crate::vars::init_cluster::ClusterInit;

pub struct Params {
    pub ref_mtx: String,
    pub alt_mtx: String,
    pub barcodes: String,
    pub num_clusters: usize,
    pub min_alt: u32,
    pub min_ref: u32,
    pub min_alt_umis: u32,
    pub min_ref_umis: u32,
    pub restarts: u32,
    pub known_cell_assignments: Option<String>,
    pub known_genotypes: Option<String>,
    pub known_genotypes_sample_names: Vec<String>,
    pub initialization_strategy: ClusterInit,
    pub threads: usize,
    pub seed: u8,
}