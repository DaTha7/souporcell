use crate::vars::{config_params::*,
                  init_cluster::*,
                  data_cell::*};

use crate::utils::sys::read::*;

use std::{fs::File,
          io::{BufRead,
               BufReader}};

use clap::{App, load_yaml};

use hashbrown::{HashMap, HashSet};

use itertools::izip;

pub fn load_params() -> Params  {
    let yaml = load_yaml!("../../yaml/params.yml");
    let params = App::from_yaml(yaml).get_matches();

    let ref_mtx = params.value_of("ref_matrix").unwrap_or("/home/hvu/hvsp/data/A_SRR5398235/runs/default_run/ref.mtx");
    let alt_mtx = params.value_of("alt_matrix").unwrap_or("/home/hvu/hvsp/data/A_SRR5398235/runs/default_run/alt.mtx");
    let barcodes = params.value_of("barcodes").unwrap_or("/home/hvu/hvsp/data/A_SRR5398235/GSM2560245_barcodes.tsv");

    let num_clusters = params.value_of("num_clusters").unwrap_or("1");
    let num_clusters = num_clusters.to_string().parse::<usize>().unwrap();

    let min_alt = params.value_of("min_alt").unwrap_or("10");
    let min_alt = min_alt.to_string().parse::<u32>().unwrap();

    let min_ref = params.value_of("min_ref").unwrap_or("10");
    let min_ref = min_ref.to_string().parse::<u32>().unwrap();

    let restarts = params.value_of("restarts").unwrap_or("100");
    let restarts = restarts.to_string().parse::<u32>().unwrap();

    let known_cell_assignments = params.value_of("known_cell_assignments");
    let known_cell_assignments = match known_cell_assignments {
        Some(x) => Some(x.to_string()),
        None => None,
    };

    let known_genotypes = params.value_of("known_genotypes");
    let known_genotypes = match known_genotypes {
        Some(x) => {
            assert_eq!(known_cell_assignments, None, "Cannot set both known_genotypes and known_cell_assignments");
            Some(x.to_string())
        },
        None => None,
    };

    let known_genotypes_sample_names = params.values_of("known_genotypes_sample_names");
    let known_genotypes_sample_names: Vec<&str> = match known_genotypes_sample_names {
        Some(x) => x.collect(),
        None => Vec::new(),
    };

    let mut sample_names: Vec<String> = Vec::new();

    for name in known_genotypes_sample_names {
        sample_names.push(name.to_string());
    }

    let initialization_strategy = params.value_of("initialization_strategy").unwrap_or("random_uniform");
    let initialization_strategy = match initialization_strategy {
        "kmeans++" => ClusterInit::KmeansPP,
        "random_uniform" => ClusterInit::RandomUniform,
        "random_cell_assignment" => ClusterInit::RandomAssignment,
        "middle_variance" => ClusterInit::MiddleVariance,
        _ => {
            assert!(false, "initialization strategy must be one of kmeans++, random_uniform, random_cell_assignment, middle_variance");
            ClusterInit::RandomAssignment
        },
    };

    let threads = params.value_of("threads").unwrap_or("32");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap_or("4");
    let seed = seed.to_string().parse::<u8>().unwrap();

    let min_ref_umis = params.value_of("min_ref_umis").unwrap_or("0");
    let min_ref_umis = min_ref_umis.to_string().parse::<u32>().unwrap();

    let min_alt_umis = params.value_of("min_alt_umis").unwrap_or("0");
    let min_alt_umis = min_alt_umis.to_string().parse::<u32>().unwrap();

    Params {
        ref_mtx: ref_mtx.to_string(),
        alt_mtx: alt_mtx.to_string(),
        barcodes: barcodes.to_string(),
        num_clusters,
        min_alt,
        min_ref,
        restarts,
        known_cell_assignments,
        known_genotypes,
        known_genotypes_sample_names: sample_names,
        initialization_strategy,
        threads,
        seed,
        min_alt_umis,
        min_ref_umis,
    }
}

pub fn load_barcodes(params: &Params) -> Vec<String> {
    let reader = reader(&params.barcodes);
    let mut cell_barcodes: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        cell_barcodes.push(line.to_string());
    }
    cell_barcodes
}

pub fn load_cell_data(params: &Params) -> (usize, usize, Vec<CellData>, Vec<usize>, HashMap<usize, usize>) {
    let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");
    let alt_reader = BufReader::new(alt_reader);

    let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");
    let ref_reader = BufReader::new(ref_reader);

    let mut used_loci: HashSet<usize> = HashSet::new();
    let mut line_number = 0;
    let mut total_loci = 0;
    let mut total_cells = 0;

    let mut all_loci: HashSet<usize> = HashSet::new();
    let mut locus_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_umi_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_counts: HashMap<usize, HashMap<usize, [u32; 2]>> = HashMap::new();

    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let alt_line = alt_line.expect("cannot read alt mtx");
        let ref_line = ref_line.expect("cannot read ref mtx");

        if line_number > 2 {
            let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
            let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
            let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1;

            all_loci.insert(locus);

            let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1;
            let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
            let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();

            assert!(locus < total_loci);
            assert!(cell < total_cells);

            let cell_counts = locus_cell_counts.entry(locus).or_insert([0; 2]);
            let umi_counts = locus_umi_counts.entry(locus).or_insert([0; 2]);

            if ref_count > 0 {
                cell_counts[0] += 1;
                umi_counts[0] += ref_count;
            }

            if alt_count > 0 {
                cell_counts[1] += 1;
                umi_counts[1] += alt_count;
            }

            let cell_counts = locus_counts.entry(locus).or_insert(HashMap::new());

            cell_counts.insert(cell, [ref_count, alt_count]);
        }
        else if line_number == 2 {
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci = tokens[0].to_string().parse::<usize>().unwrap();
            total_cells = tokens[1].to_string().parse::<usize>().unwrap();
        }

        line_number += 1;
    }
    let mut all_loci2: Vec<usize> = Vec::new();

    for loci in all_loci {
        all_loci2.push(loci);
    }

    let mut all_loci = all_loci2;
    all_loci.sort();

    let mut index_to_locus: Vec<usize> = Vec::new();
    let mut locus_to_index: HashMap<usize, usize> = HashMap::new();

    let mut cell_data: Vec<CellData> = Vec::new();

    for _cell in 0..total_cells {
        cell_data.push(CellData::new());
    }

    let mut locus_index = 0;

    for locus in all_loci {
        let cell_counts = locus_cell_counts.get(&locus).unwrap();
        let umi_counts = locus_umi_counts.get(&locus).unwrap();

        if  cell_counts[0] >= params.min_ref        &&
            cell_counts[1] >= params.min_alt        &&
            umi_counts[0] >= params.min_ref_umis    &&
            umi_counts[1] >= params.min_alt_umis    {

            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);

            for (cell, counts) in locus_counts.get(&locus).unwrap() {

                if counts[0]+counts[1] == 0 {
                    continue;
                }

                cell_data[*cell].alt_counts.push(counts[1]);
                cell_data[*cell].ref_counts.push(counts[0]);
                cell_data[*cell].loci.push(locus_index);
                cell_data[*cell].allele_fractions.push((counts[1] as f32)/((counts[0] + counts[1]) as f32));
                cell_data[*cell].log_binomial_coefficient.push(
                    statrs::
                    function::
                    factorial::
                    ln_binomial((counts[1]+counts[0]) as u64, counts[1] as u64) as f32);

                cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;

                // println!("cell {} locus {} alt {} ref {} fraction {}",*cell, locus_index, counts[1], counts[0],
                // (counts[1] as f32)/((counts[0] + counts[1]) as f32));
            }

            locus_index += 1;
        }
    }

    eprintln!("total loci used {}",used_loci.len());

    (used_loci.len(), total_cells, cell_data, index_to_locus, locus_to_index)
}