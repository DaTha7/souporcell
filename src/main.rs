mod vars;
mod utils;

use std::fs::File;
use std::io::{BufWriter, Write};

use vars::{data_thread::*,
           config_params::*};

use utils::{sys::{load::*,
                  seed::*},
            cluster::init::*,
            cluster::em::*};

use crate::vars::data_cell;
use hashbrown::HashMap;

use rand::{Rng, rngs::StdRng, SeedableRng};
use utils::sys::seed::*;
use rayon::prelude::*;

use std::time::Instant;

fn main() {

    let now: Instant = Instant::now();

    let params = load_params();
    let cell_barcodes = load_barcodes(&params);
    let (loci_used,
        total_cells,
        cell_data,
        index_to_locus,
        locus_to_index) = load_cell_data(&params);

    souporcell_main(loci_used,
                    cell_data,
                    &params,
                    cell_barcodes,
                    locus_to_index);

    eprintln!("Elapsed: {:.4?}", now.elapsed());
}

fn souporcell_main(loci_used: usize,
                   cell_data: Vec<data_cell::CellData>,
                   params: &Params,
                   barcodes: Vec<String>,
                   locus_to_index: HashMap<usize, usize>) {


    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32)/(params.threads as f32)).ceil() as usize;

    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(new_seed(&mut rng), solves_per_thread, i));
    }

    threads.par_iter_mut().for_each(|thread_data| {
        for iteration in 0..thread_data.solves_per_thread {

            let cluster_centers = init_cluster_centers(loci_used, &cell_data, params, &mut thread_data.rng, &locus_to_index);

            let (log_loss, log_probabilities) = EM(loci_used, cluster_centers, &cell_data, params, iteration, thread_data.thread_num);

            if log_loss > thread_data.best_total_log_probability {
                thread_data.best_total_log_probability = log_loss;
                thread_data.best_log_probabilities = log_probabilities;
            }

            eprintln!("thread {} iteration {} done with {}, best so far {}",
                      thread_data.thread_num, iteration, log_loss, thread_data.best_total_log_probability);

        }
    });

    let mut best_log_probability = f32::NEG_INFINITY;
    let mut best_log_probabilities: Vec<Vec<f32>> = Vec::new();

    for thread_data in threads {
        if thread_data.best_total_log_probability > best_log_probability {
            best_log_probability = thread_data.best_total_log_probability;
            best_log_probabilities = thread_data.best_log_probabilities;
        }
    }
    eprintln!("best total log probability = {}", best_log_probability);

    for (bc, log_probs) in barcodes.iter().zip(best_log_probabilities.iter()) {
        let mut best = 0;
        let mut best_lp = f32::NEG_INFINITY;

        for index in 0..log_probs.len() {
            if log_probs[index] > best_lp {
                best = index;
                best_lp = log_probs[index];
            }
        }

        print!("{}\t{}\t",bc, best);

        for index in 0..log_probs.len() {
            print!("{}",log_probs[index]);

            if index < log_probs.len() - 1 {
                print!("\t");
            }
        }

        print!("\n");
    }
}