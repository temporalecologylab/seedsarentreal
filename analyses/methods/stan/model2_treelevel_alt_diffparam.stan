functions {
  // Mean-dispersion parameterization of negative-binomial family
  // 0 < mu < +inf
  // 0 < psi    < +inf
  real neg_binomial_alt_lpmf(int n, real mu, real psi) {
    return neg_binomial_2_lpmf(n | mu, 1 / psi);
  }

  int neg_binomial_alt_rng(real mu, real psi) {
    real c = exp(-psi);
    return neg_binomial_2_rng(mu,  1 / psi);
  }
}

data {
  int<lower=1> N;           // Total number of observations
  int<lower=1> N_trees;     // Number of trees

  // Number of observations per tree
  array[N_trees] int<lower=1, upper=N> N_years;

  // Tree observation years
  array[N] real years;

  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;

  // Number of seeds in each tree each year
  array[N] int<lower=0> seed_counts;
}

parameters {
  array[N_trees] real<lower=0> lambda1;       // Non-masting intensity for tree 1
  array[N_trees] real<lower=0> psi1;          // Non-masting intensity for tree 1
  array[N_trees] real<lower=lambda1> lambda2; // Masting intensity for tree 1
  array[N_trees] real<lower=0> psi2;          // Masting dispersion for tree 1

  array[N_trees] real<lower=0, upper=1> rho0;  // Initial masting probability
  real<lower=0, upper=1> tau_nm_m; // No-masting to masting probability
  real<lower=0, upper=1> tau_m_nm; // Masting to no-masting probability
}

model {
  matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

  lambda1 ~ normal(0, 500 / 2.57); // 0 <~ lambda1 <~ 100
  psi1 ~ normal(0, 5 / 2.57); // 0 <~  psi1  <~ 2
  lambda2 ~ normal(0, 500 / 2.57); 
  psi2 ~ normal(0, 5 / 2.57); // 0 <~  psi2   <~ 5
  // Implicit uniform prior model over rho, tau_nm_m, tau_m_nm

  for (t in 1:N_trees) {

    matrix[2, N_years[t]] log_omega;
    for (n in 1:N_years[t]) {
      int y = seed_counts[tree_start_idxs[t] + n - 1];
      
      log_omega[1, n] = neg_binomial_alt_lpmf(y | lambda1[t], psi1[t]); // No-masting
      log_omega[2, n] = neg_binomial_alt_lpmf(y | lambda2[t], psi2[t]); // Masting
      
    }

    target += hmm_marginal(log_omega, Gamma, [1 - rho0[t], rho0[t]]');
  }
}

generated quantities {
  array[N] int<lower=1, upper=2> state_pred;
  array[N] int<lower=0> seed_count_pred;
  //array[N_years] real<lower=0, upper=1> p_masting;
  array[N] real p_masting;
  {
    matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                           [tau_m_nm, 1 - tau_m_nm] ];

    for (t in 1:N_trees) {
      array[N_years[t]] int tree_idxs
        = linspaced_int_array(N_years[t],
                              tree_start_idxs[t],
                              tree_end_idxs[t]);

      matrix[2, N_years[t]] log_omega;
      for (n in 1:N_years[t]) {
        int y = seed_counts[tree_start_idxs[t] + n - 1];
        log_omega[1, n] = neg_binomial_alt_lpmf(y | lambda1[t], psi1[t]); 
        log_omega[2, n] = neg_binomial_alt_lpmf(y | lambda2[t], psi2[t]); 
        
      }

      state_pred[tree_idxs]
        = hmm_latent_rng(log_omega, Gamma, [1 - rho0[t], rho0[t]]');

      for (n in 1:N_years[t]) {
        int idx = tree_start_idxs[t] + n - 1;
        if(state_pred[idx] == 1) {
          seed_count_pred[idx] = neg_binomial_alt_rng(lambda1[t], psi1[t]);
        } else if(state_pred[idx] == 2) {
          seed_count_pred[idx] = neg_binomial_alt_rng(lambda2[t], psi2[t]);
        }
      }

      p_masting[tree_idxs]
        = to_array_1d(hmm_hidden_state_prob(log_omega,
                                            Gamma,
                                            [1 - rho0[t], rho0[t]]')[2,]);
    }
  }
}
