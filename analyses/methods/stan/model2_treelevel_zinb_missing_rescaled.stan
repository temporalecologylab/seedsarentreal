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
  int<lower=1> N_max_years; // total number of years 

  // Number of observations per tree
  array[N_trees] int<lower=1, upper=N> N_years;

  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;

  // Number of seeds in each tree each year
  array[N] int<lower=0> seed_counts; // -99 is a missing observation
  
  // Tree observation years
  array[N] int<lower=1, upper=N_max_years> years;
  
}

parameters {
  real<lower=0> lambda1; // Non-masting intensity for tree 1
  real<lower=0> psi1;          // Masting dispersion for tree 1
  real<lower=0, upper=1> theta1; // probability of drawing a zero (zero-inflation)
  
  real<lower=lambda1> lambda2; // Masting intensity for tree 1
  real<lower=0> psi2;          // Masting dispersion for tree 1

  real<lower=0, upper=1> rho0;  // Initial masting probability
  real<lower=0, upper=1> tau_nm_m; // No-masting to masting probability
  real<lower=0, upper=1> tau_m_nm; // Masting to no-masting probability
}

transformed parameters {
  
  simplex[2] rho = [1 - rho0, rho0]';
  
}

model {

  lambda1 ~ normal(0, 20 / 2.57); 
  psi1 ~ normal(0, 5 / 2.57); 
  lambda2 ~ normal(0, 500 / 2.57); 
  psi2 ~ normal(0, 5 / 2.57); 
  // Implicit uniform prior model over rho0 and theta1

  for (t in 1:N_trees) {
    
    int N_years_tree = N_years[t];
    array[N_years_tree] int seed_counts_tree = seed_counts[tree_start_idxs[t]:tree_end_idxs[t]];
    array[N_years_tree] int years_tree = years[tree_start_idxs[t]:tree_end_idxs[t]];
    
    
    // Observational model
    {
      // Construct transition matrix
      matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

      // Forward algorithm
      real norm;
      real log_norm = 0;

      vector[2] alpha = rho;

      norm = max(alpha);
      log_norm += log(norm);
      alpha /= norm;

      for (n in 1:N_years_tree) {
        
        int y = seed_counts_tree[n];
        
        vector[2] log_omega;
        // Non-masting
        if (y == 0){
          log_omega[1] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1) 
          + neg_binomial_alt_lpmf(y | lambda1, psi1));
        }else{
          log_omega[1] = bernoulli_lpmf(0 | theta1) + neg_binomial_alt_lpmf(y | lambda1, psi1);
        }
        // Masting
        log_omega[2] = neg_binomial_alt_lpmf(y | lambda2, psi2); 
        
        vector[2] omega;
        omega[1] = exp(log_omega[1]);
        omega[2] = exp(log_omega[2]);
        
        if (n > 1) {
          int delta = years_tree[n] - years_tree[n - 1];
          for (d in 1:delta)
            alpha = Gamma' * alpha;
          alpha = omega .* alpha;
        }
        else {
          for(i in 1:years_tree[1])
            alpha = Gamma' * alpha;
          alpha = omega .* alpha;
        }

        norm = max(alpha);
        log_norm += log(norm);
        alpha /= norm;
      }

      // Marginal observational model
      target += sum(alpha) + log_norm;
    }
    
  }
  
}

generated quantities {
  array[N_max_years*N_trees] int states_pred;       // Posterior latent states
  array[N_max_years*N_trees] real seed_counts_pred; // Posterior predictive observations
  
  for (t in 1:N_trees) {
    
    int N_years_tree = N_years[t];
    array[N_years_tree] int seed_counts_tree = seed_counts[tree_start_idxs[t]:tree_end_idxs[t]];
    array[N_years_tree] int years_tree = years[tree_start_idxs[t]:tree_end_idxs[t]];
    int global_start_idx = (1+N_max_years*(t-1));
    int global_end_idx = (N_max_years*t);
    
    matrix[2, N_max_years] omega = rep_matrix(1, 2, N_max_years);
    
    {
    // Construct transition matrix
    matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

    // Forward algorithm, jumping across any gaps of unobserved
    
    // Build omega with the observed latent states
    for(n in 1:N_years_tree){
      vector[2] log_omega;
      int idx = years_tree[n];
      int y = seed_counts_tree[n];
      
      // Non-masting
      if (y == 0){
        log_omega[1] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1) 
        + neg_binomial_alt_lpmf(y | lambda1, psi1));
      }else{
        log_omega[1] = bernoulli_lpmf(0 | theta1) + neg_binomial_alt_lpmf(y | lambda1, psi1);
      }
      // Masting
      log_omega[2] = neg_binomial_alt_lpmf(y | lambda2, psi2); 
      omega[1,idx] = exp(log_omega[1]);
      omega[2,idx] = exp(log_omega[2]);
    }
    
    // Loop over all years
    array[N_max_years] vector[2] alpha;
    for (n in 1:N_max_years) {
      if(n == 1){
        alpha[n] = omega[,n] .* (Gamma' * rho);
      }
      else{
        alpha[n] = omega[,n] .* (Gamma' * alpha[n - 1]);
      }
      alpha[n] /= max(alpha[n]);
    }

    // Sample final latent state
    vector[2] r = alpha[N_max_years];
    vector[2] lambda = r / sum(r);

    states_pred[global_end_idx] = categorical_rng(lambda);
    if(states_pred[global_end_idx] == 1) {
      if(bernoulli_rng(theta1)){
        seed_counts_pred[global_end_idx] = 0;
      }else{
        seed_counts_pred[global_end_idx] = neg_binomial_alt_rng(lambda1, psi1);
      }
    } 
    else if(states_pred[global_end_idx] == 2) {
        seed_counts_pred[global_end_idx] = neg_binomial_alt_rng(lambda2, psi2);
    }

    // Sample latent states while running the backward algorithm
    vector[2] beta = ones_vector(2);

    for (rn in 2:N_max_years) {
      int n = N_max_years + 1 - rn;
      int global_n = global_end_idx + 1 - rn;
      int states_pred_prev = states_pred[n + 1];

      vector[2] omega_prev = omega[, n + 1];

      // matrix[K, K] Gamma_local;
      // if (n > 1)
      //   Gamma_local = matrix_power(Gamma, state[n] - state[n - 1]);
      // else
      //   Gamma_local = Gamma;

      r =   beta[states_pred_prev] * omega_prev[states_pred_prev]
          * ( alpha[n] .* Gamma[,states_pred_prev] );
      lambda = r / sum(r);

      states_pred[global_n] = categorical_rng(lambda);
      
      if(states_pred[global_n] == 1) {
        if(bernoulli_rng(theta1)){
          seed_counts_pred[global_n] = 0;
        }else{
          seed_counts_pred[global_n] = neg_binomial_alt_rng(lambda1, psi1);
        }
      } 
      else if(states_pred[global_n] == 2) {
        seed_counts_pred[global_n] = neg_binomial_alt_rng(lambda2, psi2);
      }

      // Gamma_local = matrix_power(Gamma, state[n + 1] - state[n]);
      beta = Gamma * (omega_prev .* beta);
      beta /= max(beta);
      
    }
    }
  }
 
}

