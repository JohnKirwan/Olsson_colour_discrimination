// generated with brms 2.12.0
functions {
}
data {
  int<lower=1> N;  // number of observations
  int Y[N];  // response variable
  int<lower=1> K_base;  // number of population-level effects
  matrix[N, K_base] X_base;  // population-level design matrix
  int<lower=1> K_lapse;  // number of population-level effects
  matrix[N, K_lapse] X_lapse;  // population-level design matrix
  int<lower=1> K_threshold;  // number of population-level effects
  matrix[N, K_threshold] X_threshold;  // population-level design matrix
  int<lower=1> K_width;  // number of population-level effects
  matrix[N, K_width] X_width;  // population-level design matrix
  // covariate vectors for non-linear functions
  vector[N] C_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_lapse_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  int<lower=1> J_2[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_lapse_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  int<lower=1> J_3[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_threshold_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  int<lower=1> J_4[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_threshold_1;
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  int<lower=1> J_5[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_5_width_1;
  // data for group-level effects of ID 6
  int<lower=1> N_6;  // number of grouping levels
  int<lower=1> M_6;  // number of coefficients per level
  int<lower=1> J_6[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_6_width_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector<lower=0.25,upper=0.75>[K_base] b_base;  // population-level effects
  vector<upper=-1>[K_lapse] b_lapse;  // population-level effects
  vector[K_threshold] b_threshold;  // population-level effects
  vector[K_width] b_width;  // population-level effects
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // standardized group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector[N_3] z_3[M_3];  // standardized group-level effects
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // standardized group-level effects
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  vector[N_5] z_5[M_5];  // standardized group-level effects
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  vector[N_6] z_6[M_6];  // standardized group-level effects
}
transformed parameters {
  vector[N_1] r_1_lapse_1;  // actual group-level effects
  vector[N_2] r_2_lapse_1;  // actual group-level effects
  vector[N_3] r_3_threshold_1;  // actual group-level effects
  vector[N_4] r_4_threshold_1;  // actual group-level effects
  vector[N_5] r_5_width_1;  // actual group-level effects
  vector[N_6] r_6_width_1;  // actual group-level effects
  r_1_lapse_1 = (sd_1[1] * (z_1[1]));
  r_2_lapse_1 = (sd_2[1] * (z_2[1]));
  r_3_threshold_1 = (sd_3[1] * (z_3[1]));
  r_4_threshold_1 = (sd_4[1] * (z_4[1]));
  r_5_width_1 = (sd_5[1] * (z_5[1]));
  r_6_width_1 = (sd_6[1] * (z_6[1]));
}
model {
  // initialize linear predictor term
  vector[N] nlp_base = X_base * b_base;
  // initialize linear predictor term
  vector[N] nlp_lapse = X_lapse * b_lapse;
  // initialize linear predictor term
  vector[N] nlp_threshold = X_threshold * b_threshold;
  // initialize linear predictor term
  vector[N] nlp_width = X_width * b_width;
  // initialize non-linear predictor term
  vector[N] mu;
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_lapse[n] += r_1_lapse_1[J_1[n]] * Z_1_lapse_1[n] + r_2_lapse_1[J_2[n]] * Z_2_lapse_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_threshold[n] += r_3_threshold_1[J_3[n]] * Z_3_threshold_1[n] + r_4_threshold_1[J_4[n]] * Z_4_threshold_1[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    nlp_width[n] += r_5_width_1[J_5[n]] * Z_5_width_1[n] + r_6_width_1[J_6[n]] * Z_6_width_1[n];
  }
  for (n in 1:N) {
    // compute non-linear predictor values
    mu[n] = nlp_base[n] + (1 - inv_logit(nlp_lapse[n]) - nlp_base[n]) * inv_logit(0 + 4.39 * (C_1[n] - exp(nlp_threshold[n])) / (exp(nlp_width[n])));
  }
  // priors including all constants
  target += beta_lpdf(b_base | 250, 250)
    - 1 * log_diff_exp(beta_lcdf(0.75 | 250, 250), beta_lcdf(0.25 | 250, 250));
  target += normal_lpdf(b_lapse | -3, 10)
    - 1 * normal_lcdf(-1 | -3, 10);
  target += normal_lpdf(b_threshold[1] | 0, 3);
  target += normal_lpdf(b_threshold[2] | 0, 3);
  target += normal_lpdf(b_threshold[3] | 0, 3);
  target += normal_lpdf(b_threshold[4] | 0, 3);
  target += normal_lpdf(b_threshold[5] | 0, 3);
  target += normal_lpdf(b_threshold[6] | 0, 3);
  target += normal_lpdf(b_threshold[7] | 0, 3);
  target += normal_lpdf(b_threshold[8] | 0, 3);
  target += normal_lpdf(b_width[1] | 0, 3);
  target += normal_lpdf(b_width[2] | 0, 3);
  target += normal_lpdf(b_width[3] | 0, 3);
  target += normal_lpdf(b_width[4] | 0, 3);
  target += normal_lpdf(b_width[5] | 0, 3);
  target += normal_lpdf(b_width[6] | 0, 3);
  target += normal_lpdf(b_width[7] | 0, 3);
  target += normal_lpdf(b_width[8] | 0, 3);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  target += student_t_lpdf(sd_2 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_2[1] | 0, 1);
  target += student_t_lpdf(sd_3 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_3[1] | 0, 1);
  target += student_t_lpdf(sd_4 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_4[1] | 0, 1);
  target += student_t_lpdf(sd_5 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_5[1] | 0, 1);
  target += student_t_lpdf(sd_6 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_6[1] | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += bernoulli_lpmf(Y | mu);
  }
}
generated quantities {
}

