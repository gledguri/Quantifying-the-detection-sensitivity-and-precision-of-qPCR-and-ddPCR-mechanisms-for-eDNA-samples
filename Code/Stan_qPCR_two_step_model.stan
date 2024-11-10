data { 
  //Numbers of dimentions
  // // // qPCR
  int N_st_q; // Total number of observation in qPCR standard samples
  int N_en_q; // Total number of observation in qPCR environmental samples
  int N_st_qp; // Total number of observation in qPCR standard samples for only detected samples
  int N_en_qp; // Total number of observation in qPCR environmental samples for only detected samples
  // // // ddPCR
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
  // 
  // Data
  // // qPCR
  // // // Binomial model
  array[N_st_q] int Z_qPCR_st; // Presence/Absence response of qPCR standard data
  array[N_en_q] int Z_qPCR_en; // Presence/Absence response of qPCR environmental data
  array[N_st_q] real C_qPCR_st; // Known concentration (log10) in qPCR data
  // // // Continious model
  array[N_st_qp] real Y_qPCR_st; // Ct values of qPCR standard data for only detected samples
  array[N_en_qp] real Y_qPCR_en; // Ct values of qPCR environmental data for only detected samples
  array[N_st_qp] real C_qPCR_st_cont; // Known concentration (log10) in qPCR data for only detected samples
}
parameters {
  // Parameters
  // // qPCR
  // // // Bernoulli model
  vector[N_i] phi_0;
  vector[N_i] phi_1;
  // // // Continous model (joint)
  vector[N_i] beta_joint_0;
  vector[N_i] beta_joint_1;
  vector[N_i] gamma_joint_0;
  vector<upper=0>[N_i] gamma_joint_1;
  vector[N_ij] C_qPCR_joint; 
  // // // Continous model (single)
  vector[N_i] beta_cont_0;
  vector[N_i] beta_cont_1;
  vector[N_i] gamma_cont_0;
  vector<upper=0>[N_i] gamma_cont_1;
  vector[N_ij] C_qPCR_continuous;
}
transformed parameters{
  // Parameters
  // // qPCR
  // // // Bernoulli model
  vector[N_st_q] theta_st;
  vector[N_en_q] theta_un;
  // // // Continious model (joint)
  vector[N_st_qp] mu_st;
  vector[N_en_qp] mu_en;
  vector[N_st_qp] sigma_st;
  vector[N_en_qp] sigma_en;
  // // // Continious model (single)
  vector[N_st_qp] mu_st_nor;
  vector[N_en_qp] mu_en_nor;
  vector[N_st_qp] sigma_st_nor;
  vector[N_en_qp] sigma_en_nor;
  // 
  // Model TP
  // // qPCR model
  // // // Bernuli module model compartment
  // // // // // Standard
  for (i in 1:N_st_q){
    theta_st[i] = phi_0[i_qst_idx[i]] + (phi_1[i_qst_idx[i]] * C_qPCR_st[i]);
  }
  // // // // // Unknown
  for (i in 1:N_en_q){
    theta_un[i] = phi_0[i_qen_idx[i]] + (phi_1[i_qen_idx[i]] * C_qPCR_joint[ij_qen_idx[i]]);
  }
  // // // Continious model compartment (joint)
  // // // // Standard
  for (i in 1:N_st_qp){
    mu_st[i] = beta_joint_0[i_qst_p_idx[i]] + (beta_joint_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]);
    sigma_st[i] = exp(gamma_joint_0[i_qst_p_idx[i]]+(gamma_joint_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]));
  }
  // // // // Unknown
  for (i in 1:N_en_qp){
    mu_en[i] = beta_joint_0[i_qen_p_idx[i]] + (beta_joint_1[i_qen_p_idx[i]] * C_qPCR_joint[ij_qen_p_idx[i]]);
    sigma_en[i] = exp(gamma_joint_0[i_qen_p_idx[i]]+(gamma_joint_1[i_qen_p_idx[i]] * C_qPCR_joint[ij_qen_p_idx[i]]));
  }
  // // // Continious model compartment (single)
  // // // // Standard
  for (i in 1:N_st_qp){
    mu_st_nor[i] = beta_cont_0[i_qst_p_idx[i]] + (beta_cont_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]);
    sigma_st_nor[i] = exp(gamma_cont_0[i_qst_p_idx[i]]+(gamma_cont_1[i_qst_p_idx[i]] * C_qPCR_st_cont[i]));
  }
  // // // // Unknown
  for (i in 1:N_en_qp){
    mu_en_nor[i] = beta_cont_0[i_qen_p_idx[i]] + (beta_cont_1[i_qen_p_idx[i]] * C_qPCR_continuous[ij_qen_p_idx[i]]);
    sigma_en_nor[i] = exp(gamma_cont_0[i_qen_p_idx[i]]+(gamma_cont_1[i_qen_p_idx[i]] * C_qPCR_continuous[ij_qen_p_idx[i]]));
  }
}
model {
  // Model 
  // // qPCR
  // // // Bernoulli model
  Z_qPCR_st ~ bernoulli(inv_logit(theta_st)); //Standards
  Z_qPCR_en ~ bernoulli(inv_logit (theta_un)); //Unknown (Env samples)
  // // // Continuous (Ct) model compartment (joint)
  Y_qPCR_st ~ normal(mu_st,sigma_st);//Standards
  Y_qPCR_en ~ normal(mu_en,sigma_en);//Unknown (Env samples)
  // // // Continuous (Ct) model compartment (single)
  Y_qPCR_st ~ normal(mu_st_nor,sigma_st_nor);//Standards
  Y_qPCR_en ~ normal(mu_en_nor,sigma_en_nor);//Unknown (Env samples)
  // 
  // Priors
  // // qPCR
  // // // Bernoulli model
  phi_0 ~ normal(0, 2);
  phi_1 ~ normal(0, 2);
  // // // Continious model (joint)
  beta_joint_0 ~ normal(0, 3);
  beta_joint_1 ~ normal(-3, 0.1);
  gamma_joint_0 ~ normal(0, 0.1);
  gamma_joint_1 ~ normal(-1, 0.1);
  C_qPCR_joint ~ normal(0,1);
  // // // Continious model (single)
  beta_cont_0 ~ normal(0, 3);
  beta_cont_1 ~ normal(-3, 0.1);
  gamma_cont_0 ~ normal(0, 0.1);
  gamma_cont_1 ~ normal(-1, 0.1);
  C_qPCR_continuous ~ normal(0,1);
}
