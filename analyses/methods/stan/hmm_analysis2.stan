data {
  int<lower=0> I;  // Number of latent states
  int<lower=0> N;  // Number of observations

  // All observations
  array[N] int<lower=0> y;

  // Organization of observations across latent states
  array[I] int<lower=0, upper=N> state_N;
  array[I] int<lower=1, upper=N> state_start_idx;
  array[I] int<lower=1, upper=N> state_end_idx;
}

transformed data {
  int<lower=0> K = 4;  // Number of component observational models

  // z1: Detector status
  // z2: Active component observational model
  // z = 1: (z1 = 1, z2 = 1) Working detector, low environmental activity
  // z = 2: (z1 = 1, z2 = 2) Working detector, high environmental activity
  // z = 3: (z1 = 2, z2 = 1) Broken detector,  low environmental activity
  // z = 4: (z1 = 2, z2 = 2) Broken detector,  high environmental activity
}

parameters {
  // Component observational model configurations
  real<lower=0>       lambda1; // Low environmental activity
  real<lower=lambda1> lambda2; // High environmental activity

  // Hidden state dynamics
  simplex[K] rho;

  // Low activity to low activity transition probability
  real<lower=0, upper=1> tau11;

  // High activity to high activity transition probability
  real<lower=0, upper=1> tau22;

  // Working detector to working detector transition probability
  real<lower=0, upper=1> taud;
}

model {
  // Prior model
  target += normal_lpdf(lambda1 | 0,  5 / 2.57); // 0 <~ lambda1 <~ 5
  target += normal_lpdf(lambda2 | 0, 15 / 2.57); // 0 <~ lambda2 <~ 15

  target += dirichlet_lpdf(rho | rep_vector(1, K));

  target += beta_lpdf(tau11 | 1, 1);
  target += beta_lpdf(tau22 | 1, 1);
  target += beta_lpdf(taud  | 1, 1);

  // Observational model
  {
    // Construct transition matrix
    matrix[2, 2] Gamma2
      = [ [ tau11,     1 - tau11 ],
          [ 1 - tau22, tau22     ] ];

    matrix[K, K] Gamma;
    Gamma[1:2, 1:2] = taud       * Gamma2;
    Gamma[1:2, 3:4] = (1 - taud) * Gamma2;
    Gamma[3:4, 1:2] = taud       * Gamma2;
    Gamma[3:4, 3:4] = (1 - taud) * Gamma2;

    // Construct component observational model evaluations
    matrix[K, I] log_omega;
    for (i in 1:I) {
      log_omega[, i] = zeros_vector(K);

      if (state_N[i] > 0) {
        for (n in state_start_idx[i]:state_end_idx[i]) {
          log_omega[1, i] += poisson_lpmf(y[n] | lambda1);
          log_omega[2, i] += poisson_lpmf(y[n] | lambda2);
          if (y[n] != 0) {
            log_omega[3, i] = negative_infinity();
            log_omega[4, i] = negative_infinity();
          }
        }
      }
    }

    target += hmm_marginal(log_omega, Gamma, rho);
  }
}

generated quantities {
  // Posterior joint latent states
  array[I] int z;

  // Posterior environmental latent states
  array[I] int z_env;

  // Posterior predictive observations
  array[N] int<lower=0> y_pred;

  {
    // Construct transition matrix
    matrix[2, 2] Gamma2
      = [ [ tau11,     1 - tau11 ],
          [ 1 - tau22, tau22     ] ];

    matrix[K, K] Gamma;
    Gamma[1:2, 1:2] = taud       * Gamma2;
    Gamma[1:2, 3:4] = (1 - taud) * Gamma2;
    Gamma[3:4, 1:2] = taud       * Gamma2;
    Gamma[3:4, 3:4] = (1 - taud) * Gamma2;

    // Construct component observational model evaluations
    matrix[K, I] log_omega;
    for (i in 1:I) {
      log_omega[, i] = zeros_vector(K);

      if (state_N[i] > 0) {
        for (n in state_start_idx[i]:state_end_idx[i]) {
          log_omega[1, i] += poisson_lpmf(y[n] | lambda1);
          log_omega[2, i] += poisson_lpmf(y[n] | lambda2);
          if (y[n] != 0) {
            log_omega[3, i] = negative_infinity();
            log_omega[4, i] = negative_infinity();
          }
        }
      }
    }

    z = hmm_latent_rng(log_omega, Gamma, rho);

    for (i in 1:I) {
      z_env[i] = (z[i] - 1) / 2 + 1;

      if (state_N[i] > 0) {
        if (z[i] == 1) {
          y_pred[state_start_idx[i]:state_end_idx[i]]
            = poisson_rng(rep_vector(lambda1, state_N[i]));
        } else if (z[i] == 2) {
          y_pred[state_start_idx[i]:state_end_idx[i]]
            = poisson_rng(rep_vector(lambda2, state_N[i]));
        } else if (z[i] > 2) {
          y_pred[state_start_idx[i]:state_end_idx[i]]
            = rep_array(0, state_N[i]);
        }
      }
    }
  }
}
