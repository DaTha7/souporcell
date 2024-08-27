// use std::f64::consts::*;
//
// const LANCZOS_COEFFICIENTS: [&'static f64; 11] = [
//     &1.000001502363886407565018998866435140371322631835938,
//     &0.464895966191246401422176859341561794281005859375000,
//     &-0.04956300218599807294594938866794109344482421875000,
//     &0.387353748287360133417678298428654670715332031250000,
//     &-2.32544894194688822608441114425659179687500000000000,
//     &8.977536034497006767196580767631530761718750000000000,
//     &-21.8474347546784883888904005289077758789062500000000,
//     &33.28472248815523926168680191040039062500000000000000,
//     &-30.7328214590961579233407974243164062500000000000000,
//     &15.70077904063509777188301086425781250000000000000000,
//     &-3.40185381879564374685287475585937500000000000000000
// ];
//
// const G: f64 = 1.0;
//
// pub fn log_gamma(rez: f64) -> f64 {
//     if (rez < 1.0) {
//         PI.ln() - (PI * rez).sin().ln() + log_gamma(1.0 - rez)
//     }
//     else {
//         let z = rez - 1.0;
//         let mut series = *LANCZOS_COEFFICIENTS[0];
//
//         for i in 1..LANCZOS_COEFFICIENTS.len() {
//             series += LANCZOS_COEFFICIENTS[i] / (z + i as f64);
//         }
//
//         let t = z + G + 0.5;
//
//         (2.0 * PI).sqrt().ln() + (z + 0.5)*t.ln() + (-t) + series.ln()
//     }
// }
//
// pub fn log_beta_function(a: f32, b: f32) -> f32 {
//     (log_gamma(a as f64) + log_gamma(b as f64) - log_gamma((a + b) as f64)) as f32
// }


pub fn log_beta_binomial_pmf(alt_count: f64, ref_count: f64, alpha: f64, beta: f64, ln_coefficient: f64) -> f64 {
    let ln_numerator = log_beta_calc(alt_count + alpha, ref_count + beta);
    let ln_denominator = log_beta_calc(alpha, beta);
    let log_likelihood = ln_coefficient + ln_numerator - ln_denominator;

    log_likelihood
}

pub fn log_beta_calc(alpha: f64, beta: f64) -> f64 {
    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma(alpha + beta);

    log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta
}