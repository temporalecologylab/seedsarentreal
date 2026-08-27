library(cmdstanr)
library(future.apply)


x <- delta  
x_grid <- seq(min(x), max(x), length = 100)   
Ngrid <- length(x_grid)


stan_code <- "
data {
  int<lower=1> N;
  vector[N] x;
  vector[N] y;
  int<lower=1> Ngrid;
  vector[Ngrid] x_grid;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(alpha + beta * x, sigma);
  alpha ~ normal(0, 5);
  beta  ~ normal(0, 5);
  sigma ~ normal(0, 2);
}
generated quantities {
  vector[Ngrid] y_pred;
  y_pred = alpha + beta * x_grid;
}
"

mod <- cmdstan_model(write_stan_file(stan_code))  
sub_samples <- util$filter_expectands(samples, ynames)
idxs <- 1:length(sub_samples[[ynames[1]]])
plan(multisession, workers = 20)
y_pred <- future_lapply(idxs, function(s){
  y <- sapply(ynames, function(n) sub_samples[[n]][s])
  data_lm <- list(N = length(x), x = delta, y = y, Ngrid = Ngrid, x_grid = x_grid)
  fit_lm <- mod$sample(data = data_lm, chains = 4, iter_warmup = 1000, iter_sampling = 500, refresh = 0)
  samples_lm <- fit_lm$draws('y_pred', format = "matrix")
  return(samples_lm)
}, future.seed=TRUE)
plan(sequential);gc()
y_pred <- do.call(rbind, y_pred)   
dim(y_pred)
y_pred_mean <- colMeans(y_pred)
rm(y_pred);gc()
print(delta)
