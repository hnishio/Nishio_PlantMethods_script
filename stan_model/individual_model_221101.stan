
data {
  int N;
  int N_ex;
  vector[N] ex;
  vector[N] Y;
}


parameters {
  vector[N] w;
  vector[N_ex] b_ex;
  real<lower=0> s_w;
  real<lower=0> s_b_ex;
  real<lower=0> s_Y;
}


transformed parameters {
  vector[N] alpha;
  for (t in 1:29) {
    alpha[t] = w[t];
  }
  for (t in 30:119) {
    alpha[t] = w[t] + b_ex[t-29] * ex[t];
  }
}


model {
  // State equation of w
  for (t in 1:N) {
    w[t] ~ normal(0, s_w);
  }
  
  // Coefficient of ex
  b_ex[1] ~ normal(0, s_b_ex);
  for (t in 2:N_ex) {
    b_ex[t] ~ normal(b_ex[t-1], s_b_ex);
  }
  
  // Observation equation of Y
  for (t in 1:N) {
    Y[t] ~ normal(alpha[t], s_Y);
  }
}

