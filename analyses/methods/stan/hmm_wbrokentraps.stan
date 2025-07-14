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

  // Size of trap
  // array[N_traps] real<lower=0> sizes;

  // Trap observation years
  array[N] real years;

  // Ragged array indexing
  array[N_traps] int<lower=1, upper=N> trap_start_idxs;
  array[N_traps] int<lower=1, upper=N> trap_end_idxs;

  // Number of seeds in each trap each year
  array[N] int<lower=0> seed_counts;
}

transformed data {
  int<lower=0> K = 4;  // Number of component observational models

  // z1: Detector status
  // z2: Active component observational model
  // z = 1: (z1 = 1, z2 = 1) Working detector, low environmental activity
  // z = 2: (z1 = 1, z2 = 2) Working detector, high environmental activity
  // z = 3: (z1 = 2, z2 = 1) Broken detector,  low environmental activity
  // z = 4: (z1 = 2, z2 = 2) Broken detector,  high environmental activity
  
  // Defines interpretation of lambda1, lambda2, and psi2.
  // real<lower=0> base_size = sizes[1];
}

parameters {
  // Component observational model configurations
  real<lower=0> lambda1;       // Non-masting intensity
  // real<lower=0> psi1;          // Non-masting dispersion
  real<lower=lambda1> lambda2; // Masting intensity 
  // real<lower=0> psi2;          // Masting dispersion

  // Hidden state dynamics
  simplex[K] rho;

  real<lower=0, upper=1> tau_nm_m; // No-masting to masting probability
  real<lower=0, upper=1> tau_m_nm; // Masting to no-masting probability

  // Working detector to working detector transition probability
  real<lower=0, upper=1> taud;
}

model {
  // Prior model
  target += normal_lpdf(lambda1 | 0,  500 / 2.57); // 0 <~ lambda1 <~ 5
  target += normal_lpdf(lambda2 | 0, 500 / 2.57); // 0 <~ lambda2 <~ 500

  target += dirichlet_lpdf(rho | rep_vector(1, K));

  target += beta_lpdf(tau_nm_m | 1, 1);
  target += beta_lpdf(tau_m_nm | 1, 1);
  target += beta_lpdf(taud  | 1, 1);
  
  // psi1 ~ normal(0, 5 / 2.57); // 0 <~  psi2   <~ 5
  // psi2 ~ normal(0, 5 / 2.57); // 0 <~  psi2   <~ 5

  // Observational model
  {
    // Construct transition matrix
    matrix[2, 2] Gamma2 = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

    matrix[K, K] Gamma;
    Gamma[1:2, 1:2] = taud       * Gamma2;
    Gamma[1:2, 3:4] = (1 - taud) * Gamma2;
    Gamma[3:4, 1:2] = taud       * Gamma2;
    Gamma[3:4, 3:4] = (1 - taud) * Gamma2;
    
    for (t in 1:N_traps) {
      //real trap_size = sizes[t];
      //real l1 = lambda1 * trap_size / base_size;
      //real l2 = lambda2 * trap_size / base_size;

      matrix[K, N_years[t]] log_omega;
      
      for (n in 1:N_years[t]) {
        
        log_omega[, n] = zeros_vector(K);
        
        int y = seed_counts[trap_start_idxs[t] + n - 1];
        log_omega[1, n] = poisson_lpmf(y | lambda1);
        log_omega[2, n] = poisson_lpmf(y | lambda2);
        
        if (y != 0) {
            log_omega[3, n] = negative_infinity();
            log_omega[4, n] = negative_infinity();
        }
        
      }

      target += hmm_marginal(log_omega, Gamma, rho);
    }
  }
}

generated quantities {
  array[N] int<lower=1, upper=K> state_pred;
  array[N] int<lower=0> seed_count_pred;
  //array[N_years] real<lower=0, upper=1> p_masting;
  array[N] real p_masting;
  {
    matrix[2, 2] Gamma2 = [ [1 - tau_nm_m, tau_nm_m],
                           [tau_m_nm, 1 - tau_m_nm] ];
                           
    matrix[K, K] Gamma;
    Gamma[1:2, 1:2] = taud       * Gamma2;
    Gamma[1:2, 3:4] = (1 - taud) * Gamma2;
    Gamma[3:4, 1:2] = taud       * Gamma2;
    Gamma[3:4, 3:4] = (1 - taud) * Gamma2;                       
                           

    for (t in 1:N_traps) {
      array[N_years[t]] int trap_idxs
        = linspaced_int_array(N_years[t],
                              trap_start_idxs[t],
                              trap_end_idxs[t]);
                              

      matrix[K, N_years[t]] log_omega;
      
       for (n in 1:N_years[t]) {
        
        log_omega[, n] = zeros_vector(K);
        
        int y = seed_counts[trap_start_idxs[t] + n - 1];
        log_omega[1, n] = poisson_lpmf(y | lambda1);
        log_omega[2, n] = poisson_lpmf(y | lambda2);
        
        if (y != 0) {
            log_omega[3, n] = negative_infinity();
            log_omega[4, n] = negative_infinity();
        }
        
      }

      state_pred[trap_idxs]
        = hmm_latent_rng(log_omega, Gamma, rho);

      for (n in 1:N_years[t]) {
        int idx = trap_start_idxs[t] + n - 1;
        if(state_pred[idx] == 1) {
          seed_count_pred[idx] = poisson_rng(lambda1);
        } else if(state_pred[idx] == 2) {
          seed_count_pred[idx] = poisson_rng(lambda2);
        } else if(state_pred[idx] > 2) {
          seed_count_pred[idx] = 0;
        }
        
      }

      p_masting[trap_idxs]
        = to_array_1d(hmm_hidden_state_prob(log_omega,
                                            Gamma,
                                            rho)[2,]);
    }
  }
}
