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

  // Tree observation years
  array[N] real years;

  // Ragged array indexing
  array[N_trees] int<lower=1, upper=N> tree_start_idxs;
  array[N_trees] int<lower=1, upper=N> tree_end_idxs;

  // Number of seeds in each tree each year
  array[N] int<lower=0> seed_counts; // -99 is a missing observation
  
  // Index of observed years for each tree
  array[N] int<lower=1, upper=N_max_years> observed_years; 
}

transformed data {
  
  int<lower= N> N_pred = N_trees * N_max_years;
  
}

parameters {
  real<lower=0> lambda1; // Non-masting intensity for tree 1
  real<lower=0, upper=1> theta1; // probability of drawing a zero (zero-inflation)
  
  real<lower=lambda1> lambda2; // Masting intensity for tree 1
  real<lower=0> psi2;          // Masting dispersion for tree 1

  real<lower=0, upper=1> rho0;  // Initial masting probability
  real<lower=0, upper=1> tau_nm_m; // No-masting to masting probability
  real<lower=0, upper=1> tau_m_nm; // Masting to no-masting probability
}

model {
  matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
                         [tau_m_nm, 1 - tau_m_nm] ];

  lambda1 ~ normal(0, 500 / 2.57); 
  lambda2 ~ normal(0, 500 / 2.57); 
  psi2 ~ normal(0, 5 / 2.57); 
  // Implicit uniform prior model over rho, tau_nm_m, tau_m_nm, and theta1

  for (t in 1:N_trees) {

    matrix[2, N_max_years] log_omega = rep_matrix(0, 2, N_max_years);
    
    for (n in 1:N_years[t]) {
      
      int y = seed_counts[tree_start_idxs[t] + n - 1];
      int i = observed_years[tree_start_idxs[t] + n - 1];
      
      // Non-masting
      if (y == 0){
        log_omega[1, i] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1) + poisson_lpmf(y | lambda1));
      }else{
        log_omega[1, i] = bernoulli_lpmf(0 | theta1) + poisson_lpmf(y | lambda1);
      }
      
      log_omega[2, i] = neg_binomial_alt_lpmf(y | lambda2, psi2); // Masting
      
    }

    target += hmm_marginal(log_omega, Gamma, [1 - rho0, rho0]');
  }
}

// generated quantities {
//   
//   array[N_pred] int<lower=1, upper=2> state_pred;
//   array[N_pred] int<lower=0> seed_count_pred;
//   array[N_pred] real p_masting;
//   
//   {
//     matrix[2, 2] Gamma = [ [1 - tau_nm_m, tau_nm_m],
//                            [tau_m_nm, 1 - tau_m_nm] ];
// 
//     for (t in 1:N_trees) {
//       array[N_years[t]] int tree_idxs
//         = linspaced_int_array(N_years[t],
//                               tree_start_idxs[t],
//                               tree_end_idxs[t]);
// 
//       matrix[2, N_max_years] log_omega = rep_matrix(0, 2, N_max_years);
//       
//       for (n in 1:N_years[t]) {
//         int y = seed_counts[tree_start_idxs[t] + n - 1];
//         int i = observed_years[tree_start_idxs[t] + n - 1];
//         
//         if (y == 0){
//         log_omega[1, i] = log_sum_exp(bernoulli_lpmf(1 | theta1), bernoulli_lpmf(0 | theta1) + poisson_lpmf(y | lambda1));
//         }else{
//           log_omega[1, i] = bernoulli_lpmf(0 | theta1) + poisson_lpmf(y | lambda1);
//         }
//       
//         log_omega[2, i] = neg_binomial_alt_lpmf(y | lambda2, psi2); 
//         
//       }
//       
//       array[N_max_years] int tree_idxs_pred
//         = linspaced_int_array(N_max_years,
//                               1+(t-1)*N_max_years,
//                               t*N_max_years);
// 
//       state_pred[tree_idxs_pred]
//         = hmm_latent_rng(log_omega, Gamma, [1 - rho0, rho0]');
// 
//       for (n in 1:N_max_years) {
//         
//         int idx = 1+(t-1)*N_max_years + n - 1;
//         
//         if(state_pred[idx] == 1) {
//           
//           if(bernoulli_rng(theta1)){
//             seed_count_pred[idx] = 0;
//           }else{
//             seed_count_pred[idx] = poisson_rng(lambda1);
//           }
//           
//         } else if(state_pred[idx] == 2) {
//           seed_count_pred[idx] = neg_binomial_alt_rng(lambda2, psi2);
//         }
//       }
//     
//       p_masting[tree_idxs_pred]
//         = to_array_1d(hmm_hidden_state_prob(log_omega,
//                                             Gamma,
//                                             [1 - rho0, rho0]')[2,]);
//     }
//   }
// }
