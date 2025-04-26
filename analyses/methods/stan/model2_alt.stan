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
  int<lower=1> N_traps;     // Number of traps

  // Number of observations per trap
  array[N_traps] int<lower=1, upper=N> N_years;

  // Size of search area below trap
  array[N_traps] real<lower=0> sizes;

  // trap observation years
  array[N] real years;

  // Ragged array indexing
  array[N_traps] int<lower=1, upper=N> trap_start_idxs;
  array[N_traps] int<lower=1, upper=N> trap_end_idxs;

  // Number of seeds in each trap each year
  array[N] int<lower=0> seed_counts;
}

transformed data {
  // Defines interpretation of lambda1, lambda2, and psi2.
  real<lower=0> base_size = sizes[1];
}

parameters {
  real<lower=0> lambda1;       // Non-masting intensity for trap 1
  real<lower=0> psi1;          // Non-masting intensity for trap 1
  real<lower=lambda1> lambda2; // Masting intensity for trap 1
  real<lower=0> psi2;          // Masting dispersion for trap 1

  real<lower=0, upper=1> rho0;  // Initial masting probability
  real<lower=0, upper=1> tau_nm_m; // No-masting to masting probability
  real<lower=0, upper=1> tau_m_nm; // Masting to no-masting probability
}

model {
  matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

  lambda1 ~ normal(0, 10 / 2.57); // 0 <~ lambda1 <~ 10
  psi1 ~ normal(0, 5 / 2.57); // 0 <~  psi1  <~ 2
  lambda2 ~ normal(0, 500 / 2.57); 
  psi2 ~ normal(0, 5 / 2.57); // 0 <~  psi2   <~ 5
  // Implicit uniform prior model over rho, tau_nm_m, tau_m_nm

  for (t in 1:N_traps) {
    real trap_size = sizes[t];
    real l1 = lambda1 * trap_size / base_size;
    real l2 = lambda2 * trap_size / base_size;

    matrix[2, N_years[t]] log_omega;
    for (n in 1:N_years[t]) {
      int y = seed_counts[trap_start_idxs[t] + n - 1];
      
      log_omega[1, n] = neg_binomial_alt_lpmf(y | l1, psi1); // No-masting
      log_omega[2, n] = neg_binomial_alt_lpmf(y | l2, psi2); // Masting
      
    }

    target += hmm_marginal(log_omega, Gamma, [1 - rho0, rho0]');
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

    for (t in 1:N_traps) {
      array[N_years[t]] int trap_idxs
        = linspaced_int_array(N_years[t],
                              trap_start_idxs[t],
                              trap_end_idxs[t]);

      real trap_size = sizes[t];
      real l1 = lambda1 * trap_size / base_size;
      real l2 = lambda2 * trap_size / base_size;

      matrix[2, N_years[t]] log_omega;
      for (n in 1:N_years[t]) {
        int y = seed_counts[trap_start_idxs[t] + n - 1];
        log_omega[1, n] = neg_binomial_alt_lpmf(y | l1, psi1); 
        log_omega[2, n] = neg_binomial_alt_lpmf(y | l2, psi2); 
        
      }

      state_pred[trap_idxs]
        = hmm_latent_rng(log_omega, Gamma, [1 - rho0, rho0]');

      for (n in 1:N_years[t]) {
        int idx = trap_start_idxs[t] + n - 1;
        if(state_pred[idx] == 1) {
          seed_count_pred[idx] = neg_binomial_alt_rng(l1, psi1);
        } else if(state_pred[idx] == 2) {
          seed_count_pred[idx] = neg_binomial_alt_rng(l2, psi2);
        }
      }

      p_masting[trap_idxs]
        = to_array_1d(hmm_hidden_state_prob(log_omega,
                                            Gamma,
                                            [1 - rho0, rho0]')[2,]);
    }
  }
}
