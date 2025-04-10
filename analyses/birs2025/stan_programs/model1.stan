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
  array[N_traps] real<lower=0> sizes;

  // Trap observation years
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
  real<lower=0> psi1;          // Non-masting dispersion for trap 1
  real<lower=lambda1> lambda2; // Masting intensity for trap 1
  real<lower=0, upper=1> rho;  // Masting probability
}

model {
  lambda1 ~ normal(0, 500 / 2.57); // 0 <~ lambda1 <~ 500
  psi1 ~ normal(0, 5 / 2.57);      // 0 <~  psi2   <~ 5
  lambda2 ~ normal(0, 500 / 2.57); // 0 <~ lambda2 <~ 500
  // Implicit uniform prior model over rho

  for (t in 1:N_traps) {
    real trap_size = sizes[t];
    real l1 = lambda1 * trap_size / base_size;
    real l2 = lambda2 * trap_size / base_size;

    for (n in 1:N_years[t]) {
      int y = seed_counts[trap_start_idxs[t] + n - 1];
      target += log_mix(1 - rho,
                        neg_binomial_alt_lpmf(y | l1, psi1),
                        poisson_lpmf(y | l2));
    }
  }
}

generated quantities {
  array[N] int<lower=0> seed_count_pred;
  array[N] real<lower=0, upper=1> p_masting;

  for (t in 1:N_traps) {
    real trap_size = sizes[t];
    real l1 = lambda1 * trap_size / base_size;
    real l2 = lambda2 * trap_size / base_size;

    for (n in 1:N_years[t]) {
      int idx = trap_start_idxs[t] + n - 1;
      int y = seed_counts[idx];

      vector[2] xs
        = [ log(1 - rho) + neg_binomial_alt_lpmf(y | l1, psi1),
            log(rho)     + poisson_lpmf(y | l2)                 ]';
      p_masting[idx] = softmax(xs)[2];

      if(bernoulli_rng(1 - rho)) {
        seed_count_pred[idx] = neg_binomial_alt_rng(l1, psi1);
      } else {
        seed_count_pred[idx] = poisson_rng(l2);
      }
    }
  }
}
