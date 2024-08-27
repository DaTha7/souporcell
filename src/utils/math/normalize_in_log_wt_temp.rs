use crate::utils::math::log_sum_exp::*;

pub fn normalize_in_log_with_temp(log_probs: &Vec<f32>,
                              temp: f32) -> Vec<f32> {

    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let mut new_log_probs: Vec<f32> = Vec::new();
    for log_prob in log_probs {
        new_log_probs.push(log_prob/temp);
    }

    let sum = log_sum_exp(&new_log_probs);
    for i in 0..log_probs.len() {
        normalized_probabilities.push((new_log_probs[i]-sum).exp());
    }
    normalized_probabilities
}
