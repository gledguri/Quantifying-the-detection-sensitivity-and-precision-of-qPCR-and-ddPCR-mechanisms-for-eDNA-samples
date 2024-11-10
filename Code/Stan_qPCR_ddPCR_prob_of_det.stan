data { 
  //Numbers of dimentions
  int Nq; // Total number of observation in qPCR standard samples
  int Ndd; // Total number of observation in ddPCR standard samples
  int N_i; // Number of assays in both data (qPCR and ddPCR)
  //Indexes   
  array[Nq] int i_q_idx; // Species/assay index for qPCR 
  array[Ndd] int i_d_idx; // Species/assay index for ddPCR
  // Data
  array[Nq] int Z_qPCR; // Presence/Absence of targets in qPCR runs
  array[Ndd] int Z_ddPCR; // Presence/Absence of targets in ddPCR runs
  array[Nq] real C_qPCR; // Known concentration (log10) in qPCR data
  array[Ndd] real C_ddPCR; // Known concentration (log10) in qPCR data
}
parameters {
  vector[N_i] phi_qPCR_0; // Intercept of logistic regression for qPCR detection probability model
  vector[N_i] phi_qPCR_1; // Slope of logistic regression for qPCR detection probability model
  vector[N_i] phi_ddPCR_0; // Intercept of logistic regression for ddPCR detection probability model
  vector[N_i] phi_ddPCR_1; // Slope of logistic regression for ddPCR detection probability model
}
transformed parameters{
  vector[Nq] theta_qPCR; // probability of detection (logit-space) qPCR model
  vector[Ndd] theta_ddPCR; // probability of detection (logit-space) ddPCR model
  // Logistic regression for qPCR detection probability model
  for (i in 1:Nq){
    theta_qPCR[i] = phi_qPCR_0[i_q_idx[i]] + (phi_qPCR_1[i_q_idx[i]] * C_qPCR[i]);
  }
  // Logistic regression for ddPCR detection probability model
  for (i in 1:Ndd){
    theta_ddPCR[i] = phi_ddPCR_0[i_d_idx[i]] + (phi_ddPCR_1[i_d_idx[i]] * (C_ddPCR[i]+2));
    //' Since the standards of ddPCR range from 10^-3 - 10^4 while qPCR 
    //' standards range from 10^-1 - 10^6 we added 2 to C_ddPCR to make the 
    //' parameters phi_qPCR_0 and phi_qPCR_1 comparable with phi_ddPCR_0 and phi_ddPCR_1
  }
}
model {
  Z_qPCR ~ bernoulli(inv_logit(theta_qPCR)); 
  Z_ddPCR ~ bernoulli(inv_logit(theta_ddPCR));
  // Priors
  phi_qPCR_0 ~ normal(0, 1);
  phi_qPCR_1 ~ normal(0, 1);
  phi_ddPCR_0 ~ normal(0, 1);
  phi_ddPCR_1 ~ normal(0, 1);
}
