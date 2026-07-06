// Route C — Stage-2 censored factor model (dev/plan-route-c.md).
//
// Residuals r_ij (detected) / censoring bounds b_ij (BDL) from the Stage-1
// per-analyte GAM mean, modelled as a rank-K factor structure:
//   r_ij = Lambda[j] . f[i] + eps_ij,  f_i ~ N(0, I_K),  eps_ij ~ N(0, psi_j)
// BDL residual cells are left-censored data augmentation: the latent
// residual is a parameter constrained <= its censoring bound.
//
// Identifiability: the top K x K block of Lambda is fixed lower-triangular
// with a positive diagonal (standard factor-analysis constraint); rotation
// is otherwise irrelevant to prediction (see .factor_condition()).

data {
  int<lower=1> N;                       // samples
  int<lower=1> J;                       // analytes
  int<lower=1> K;                       // factors, K < J
  int<lower=0> N_obs;                   // detected residual cells
  int<lower=0> N_cens;                  // BDL residual cells
  array[N_obs] int<lower=1, upper=N> obs_row;
  array[N_obs] int<lower=1, upper=J> obs_col;
  vector[N_obs] r_obs;                  // detected residuals
  array[N_cens] int<lower=1, upper=N> cen_row;
  array[N_cens] int<lower=1, upper=J> cen_col;
  vector[N_cens] b_cens;                // residual upper bounds = log(DL) - mu
}
transformed data {
  int<lower=0> K_low = (K * (K - 1)) %/% 2; // free strictly-lower-triangular entries
  int<lower=0> J_rest = J - K;          // rows below the identifying top block
}
parameters {
  matrix[N, K] f;                       // latent per-sample factor scores
  vector<lower=0>[K] Lambda_diag;       // positive diagonal of the top block
  vector[K_low] Lambda_lower;           // free strictly-lower entries of the top block
  matrix[J_rest, K] Lambda_rest;        // free loadings for the remaining rows
  vector<lower=0>[J] psi_sd;            // idiosyncratic sd (Psi = diag(psi_sd^2))
  vector<upper=0>[N_cens] r_cens_raw;   // r_cens - b_cens <= 0 (left-censored)
}
transformed parameters {
  matrix[J, K] Lambda;
  vector[J] psi = square(psi_sd);
  {
    int idx = 1;
    for (i in 1:K) {
      for (j in 1:K) {
        if (j < i) {
          Lambda[i, j] = Lambda_lower[idx];
          idx += 1;
        } else if (j == i) {
          Lambda[i, j] = Lambda_diag[i];
        } else {
          Lambda[i, j] = 0;
        }
      }
    }
  }
  if (J_rest > 0) {
    Lambda[(K + 1):J, ] = Lambda_rest;
  }
}
model {
  to_vector(f) ~ std_normal();
  Lambda_diag  ~ std_normal();          // truncated at 0 by the <lower=0> declaration
  Lambda_lower ~ std_normal();
  to_vector(Lambda_rest) ~ std_normal();
  psi_sd ~ student_t(3, 0, 1);

  for (n in 1:N_obs) {
    real mean_n = dot_product(Lambda[obs_col[n]], f[obs_row[n]]);
    r_obs[n] ~ normal(mean_n, psi_sd[obs_col[n]]);
  }
  for (n in 1:N_cens) {
    real rc = b_cens[n] + r_cens_raw[n];
    real mean_n = dot_product(Lambda[cen_col[n]], f[cen_row[n]]);
    rc ~ normal(mean_n, psi_sd[cen_col[n]]);
  }
}
