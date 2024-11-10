data {
  //Numbers of dimensions
  // // // qPCR
  int N_st_q; // Total number of observation in qPCR standard samples
  int N_en_q; // Total number of observation in qPCR environmental samples
  int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
  int N_en_qp; // Total number of observation in qPCR environmental samples for only detected samples
  // // // ddPCR
  int N_st_d; // Total number of observation in ddPCR standard samples
  int N_en_d; // Total number of observation in qPCR environmental samples
  // // Joined
  int N_i; // Number of species in both data
  int N_ij; // Number of species and stations in both data
  //
  //Indexes
  // // qPCR
  // // // Binomial model
  array[N_st_q] int i_qst_idx; // Species index for qPCR standard samples
  array[N_en_q] int i_qen_idx; // Species index for qPCR environmental samples
  array[N_en_q] int ij_qen_idx; // Species and standard index for qPCR environmental samples
  // // // Continious model
  array[N_st_qp] int i_qst_p_idx; // Species index for qPCR standard samples
  array[N_en_qp] int i_qen_p_idx; // Species index for qPCR environmental samples
  array[N_en_qp] int ij_qen_p_idx; // Species and standard index for qPCR environmental samples
  // // ddPCR
  array[N_st_d] int i_dst_idx; // Species index for ddPCR environmental samples
  array[N_en_d] int i_den_idx; // Species index for ddPCR standard samples
  array[N_en_d] int ij_den_idx; // Species and standard index for ddPCR environmental samples
  //
  // Data
  // // qPCR
  // // // Binomial model
  array[N_st_q] int Z_qPCR_st; // Presence/Absence response of qPCR standard data
  array[N_en_q] int Z_qPCR_en; // Presence/Absence response of qPCR environmental data
  array[N_st_q] real C_qPCR_st; // Known concentration (log10) in qPCR data
  // // // Continuous model
  array[N_st_qp] real Y_qPCR_st; // Ct values of qPCR standard data for only detected samples
  array[N_en_qp] real Y_qPCR_en; // Ct values of qPCR environmental data for only detected samples
  array[N_st_qp] real C_qPCR_st_continuous; // Known concentration (log10) in qPCR data for only detected samples
  // // ddPCR
  array[N_st_d] int W_st; // Observed positive droplets in ddPCR standard samples
  array[N_en_d] int W_en; // Observed positive droplets in ddPCR environmental samples
  array[N_st_d] int U_st; // Total droplets in ddPCR standard samples
  array[N_en_d] int U_en; // Total droplets in ddPCR environmental samples
  array[N_st_d] real C_ddPCR_st; // Known concentration (log10) in ddPCR data
}
parameters {
  // Parameters
  // // qPCR
  // // // Bernoulli model
  vector[N_i] phi_0;
  vector[N_i] phi_1;
  // // // Continous model
  vector[N_i] beta_0;
  vector[N_i] beta_1;
  vector[N_i] gamma_0;
  vector<upper=0>[N_i] gamma_1;
  vector[N_ij] C_qPCR;
  // // ddPCR
  vector[N_i] kappa_0;
  vector[N_i] kappa_1;
  vector<lower=-7>[N_ij] C_ddPCR;
}
transformed parameters{
  // Parameters
  // // qPCR
  // // // Bernoulli model
  vector[N_st_q] theta_st;
  vector[N_en_q] theta_un;
  // // // Continuous model
  vector[N_st_qp] mu_st;
  vector[N_en_qp] mu_en;
  vector[N_st_qp] sigma_st;
  vector[N_en_qp] sigma_en;
  // // ddPCR
  vector[N_st_d] omega_st;
  vector[N_en_d] omega_en;
  vector[N_ij] delta; // Difference between ddPCR and qPCR concentration estimates
  //
  // Model TP
  // // qPCR model
  // // // Bernuli module model compartment
  // // // // // Standard
  for (i in 1:N_st_q){
    theta_st[i] = phi_0[i_qst_idx[i]] + (phi_1[i_qst_idx[i]] * C_qPCR_st[i]);
  }
  // // // // // Unknown (Env samples)
  for (i in 1:N_en_q){
    theta_un[i] = phi_0[i_qen_idx[i]] + (phi_1[i_qen_idx[i]] * C_qPCR[ij_qen_idx[i]]);
  }
  // // // Continuous model compartment
  // // // // Standard
  for (i in 1:N_st_qp){
    mu_st[i] = beta_0[i_qst_p_idx[i]] + (beta_1[i_qst_p_idx[i]] * C_qPCR_st_continuous[i]);
    sigma_st[i] = exp(gamma_0[i_qst_p_idx[i]]+(gamma_1[i_qst_p_idx[i]] * C_qPCR_st_continuous[i]));
  }
  // // // // Unknown (Env samples)
  for (i in 1:N_en_qp){
    mu_en[i] = beta_0[i_qen_p_idx[i]] + (beta_1[i_qen_p_idx[i]] * C_qPCR[ij_qen_p_idx[i]]);
    sigma_en[i] = exp(gamma_0[i_qen_p_idx[i]]+(gamma_1[i_qen_p_idx[i]] * C_qPCR[ij_qen_p_idx[i]]));
  }
  // ddPCR Model
  // // Standard
  for (i in 1:N_st_d){
    omega_st[i] = kappa_0[i_dst_idx[i]]+(kappa_1[i_dst_idx[i]]*C_ddPCR_st[i]);
  }
  // // // Unknown (Env samples)
  for (i in 1:N_en_d){
    omega_en[i] = kappa_0[i_den_idx[i]]+(kappa_1[i_den_idx[i]]*C_ddPCR[ij_den_idx[i]]);
  }
  for (i in 1:N_ij){
    delta[i] = C_ddPCR[i]-C_qPCR[i];
  }
}
model {
  // Model
  // // qPCR
  // // // Bernoulli model
  Z_qPCR_st ~ bernoulli(inv_logit(theta_st)); //Standards
  Z_qPCR_en ~ bernoulli(inv_logit (theta_un)); //Unknown (Env samples)
  // // // Continuous (Ct) model compartment
  Y_qPCR_st ~ normal(mu_st,sigma_st);//Standards
  Y_qPCR_en ~ normal(mu_en,sigma_en);//Unknown (Env samples)
  // // // ddPCR
  W_st ~ binomial(U_st, inv_cloglog(omega_st)); //Standards
  W_en ~ binomial(U_en, inv_cloglog(omega_en)); //Unknown (Env samples)
  //
  // Priors
  // // qPCR
  // // // Bernoulli model
  phi_0 ~ normal(0, 2);
  phi_1 ~ normal(0, 2);
  // // // Continuous model
  beta_0 ~ normal(0, 3);
  beta_1 ~ normal(-3, 0.1);
  gamma_0 ~ normal(0, 0.1);
  gamma_1 ~ normal(-1, 0.1);
  C_qPCR ~ normal(0,3);
  // // ddPCR
  kappa_0 ~ normal(0, 1);
  kappa_1 ~ normal(0, 1);
  C_ddPCR ~ normal(0,3);
}
