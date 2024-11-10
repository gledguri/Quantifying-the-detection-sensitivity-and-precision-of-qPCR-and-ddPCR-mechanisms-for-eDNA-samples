---
title: "Quantifying the detection sensitivity and precision of qPCR and ddPCR mechanisms for eDNA samples"
format: 
  html:
    toc: true
    toc-location: left
    toc-title: "Table of Contents"
    toc-depth: 2
    number-sections: true
    page-layout: full
    self-contained: true 
    collapsible: true
    collapse-level: 1
    grid:
      sidebar-width: 200px
      body-width: 880px
      margin-width: 200px
      gutter-width: 1.5rem
author: Gled Guri
date: "`r format(Sys.time(), '%d %B, %Y')`"
editor: source
editor_options: 
  chunk_output_type: console
---

<style>
body {
  font-family: Arial, sans-serif;
  font-size: 18px;
}

h1 {
  font-size: 22px;  /* Adjust the size of H1 headers */
  font-weight: bold;  /* Make H1 headers not bold */
}

h2 {
  font-size: 18px;  /* Adjust the size of H2 headers */
  font-weight: normal;  /* Make H2 headers not bold */
}

h3 {
  font-size: 18px;  /* Adjust the size of H3 headers */
  font-weight: normal;  /* Make H3 headers not bold */
}

h4 {
  font-size: 18px;  /* Adjust the size of H4 headers */
  font-weight: normal;  /* Make H4 headers not bold */
}
</style>

Load functions

```{r}
extract_list_param <- function(stanMod){
	l <- stanMod@model_pars
	x <- list(extract_param(stanMod,"lp__")) %>% setNames("lp__")
	for (i in l) {
		x <- c(x,list(extract_param(stanMod,i)) %>% setNames(i))
	}
	x <- x[-which(names(x) == c("lp__"))]
	return(x)
}

cloglog <- function(theta) log(-log(1 - theta))

extract_param <- function(model=stanMod,parmeter="alpha"){
	return(summary(model, par = parmeter)$summary %>% unlist()%>%as.data.frame%>%round(.,2))
}

inv.cloglog <- function(theta) 1-exp(-exp(theta))

logreg <- function(x,b0,b1,i) y=(1/(1+exp(-(b0+(b1*(x+i))))))

scientific_10 <- function(x) {
	c <- scales::scientific_format()(x)
	t <- gsub("1e", "10^", c)
	t2 <- gsub("10\\^\\+", "10\\^", t)
	str2expression(t2)}
```

Load libraries

```{r, message=FALSE}
library(ggplot2)
library(dplyr)
library(rstan)
library(here);options(mc.cores = parallel::detectCores())

```

Load data

```{r}
e_ddpcr <- read.csv(here('Data','ddPCR_environmental_samples.csv'))
e_qpcr <- read.csv(here('Data','qPCR_environmental_samples.csv'))
st_ddpcr <- read.csv(here('Data','ddPCR_standard_samples.csv'))
st_qpcr <- read.csv(here('Data','qPCR_standard_samples.csv'))

st_qpcr_cm <- st_qpcr %>% filter(!is.na(Ct))
e_qpcr_cm <- e_qpcr %>% filter(!is.na(Ct))

sample_metadata <- read.csv(here('Data','Sample_metadata.csv'))
stidx <- read.csv(here('Data','sample_idx.csv'))
```

# Probability of detection

## Create stan data

```{r}
stan_data_M1 <- list(
	N_i = length(unique(st_qpcr$Species)),

	Nq = nrow(st_qpcr),
	Ndd = nrow(st_ddpcr),

	i_q_idx = st_qpcr$species_idx,
	i_d_idx = st_ddpcr$species_idx,

	Z_qPCR = st_qpcr$pres,
	Z_ddPCR = st_ddpcr$pres,

	C_qPCR = log10(st_qpcr$int_concentation),
	C_ddPCR = log10(st_ddpcr$int_concentation)
);str(stan_data_M1)

initial_values <- replicate(4,list(
	alpha_0 = rep(1.0, 3),
	alpha_1 = rep(0.4, 3),
	beta_0 = rep(2.50, 3),
	beta_1 = rep(1.23, 3)
),simplify=F)
```

## The math

$$
\begin{equation}
\begin{aligned}
Z_{ijk} &\sim \text{Bernoulli}(\theta_{ij}) && \text{(1)} \\
\text{logit}(\theta_{ij}) &= \phi0_{i} + \phi1_{i} \times \log_{10}(C_{ij}) && \text{(2)}
\end{aligned}
\end{equation}
$$ \## The stan code

This stan code runs 2 independent logistic regression analyses (equation 1 and 2) of qPCR and ddPCR data, with separate sets of parameters (intercepts and slopes) for each technique. To ensure comparability between the parameters $φ0_{qPCR}$, $φ1_{qPCR}$ and $φ0_{ddPCR}$, $φ1_{ddPCR}$ , we adjusted $C_{ddPCR}$ by adding 2 in stan ($10^2$ copies/μL) since the concentration of ddPCR standards ranged from $10^{-3}$ to $10^4$ , while those for qPCR range from $10^{-1}$ to $10^6$ . The shift aligns the slope and intercept enabling to compare the two models through their parameters ($φ0$ and $φ1$).

```{cpp}
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
```

## Run stan Model

```{r, message=FALSE, results='hide'}
stanMod_M1 = stan(file = here('Code','Stan_qPCR_ddPCR_prob_of_det.stan'),
									chains = 4,
									model_name = 'qPCR_vs_ddPCR_probability_of_detection',
									iter = 5000,
									warmup = 2000,
									init = initial_values,
									data = stan_data_M1)
```

## Data for Figure 1

```{r, message=FALSE}
# Extract parameters from stan model to build the graph
mod_out <- extract_list_param(stanMod_M1)
mod_out <- mod_out[which(names(mod_out) %in% c("phi_qPCR_0","phi_qPCR_1","phi_ddPCR_0","phi_ddPCR_1"))]

# Add species names
for (i in names(mod_out)) {
	mod_out[[i]] <- mod_out[[i]] %>%
		cbind(.,tibble(sp_idx = 1:3, Species = c("Cod", "Herring", "Saithe")))
}

# Construct the logistic regression for qPCR probability of detection from phi_qPCR_0 & phi_qPCR_1 parameters
probability_of_detection_qpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$phi_qPCR_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$phi_qPCR_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,0))

# Construct the logistic regression for ddPCR probability of detection from alpha_0 & alpha_1 parameters
probability_of_detection_ddpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$phi_ddPCR_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$phi_ddPCR_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,2)) #added 2 to logistic regression (see why in stan model file)

# Combine the logistic regresion lines (prob of det) of both data (qPCR and ddPCR)
probability_of_detection_data <- cbind(probability_of_detection_qpcr %>% select(C,Species,prob) %>% rename(prob_q=prob),
										probability_of_detection_ddpcr %>% select(prob) %>% rename(prob_dd=prob)) %>%
	mutate(delta=prob_dd-prob_q)

probability_of_detection_data %>% as_tibble() %>% print(n=20)
```

## Plot Figure 1 - Detection probability (qPCR vs ddPCR)

```{r, message=FALSE, warning=FALSE}
#| fig-width: 16
#| fig-height: 8

source(here('Code','Surpressed','Figure_1.R'))
Figure_1
```

# Quantification precision

## Create stan data

```{r}
stan_data_M2 <- list(
	N_st_q = nrow(st_qpcr),
	N_en_q = nrow(e_qpcr),
	N_st_qp = nrow(st_qpcr_cm),
	N_en_qp = nrow(e_qpcr_cm),
	#
	N_st_d = nrow(st_ddpcr),
	N_en_d = nrow(e_ddpcr),
	#
	N_i = length(unique(st_qpcr$Species)),
	N_ij = length(unique(e_ddpcr$sp_st_idx)),
	#
	i_qst_idx = st_qpcr$species_idx,
	i_qen_idx = e_qpcr$species_idx,
	ij_qen_idx = e_qpcr$sp_st_idx,
	
	i_qst_p_idx = st_qpcr_cm$species_idx,
	i_qen_p_idx = e_qpcr_cm$species_idx,
	ij_qen_p_idx = e_qpcr_cm$sp_st_idx,
	
	Z_qPCR_st = st_qpcr$pres,
	Z_qPCR_en = e_qpcr$pres,
	C_qPCR_st = log10(st_qpcr$int_concentation),
	
	Y_qPCR_st = st_qpcr_cm$Ct,
	Y_qPCR_en = e_qpcr_cm$Ct,
	C_qPCR_st_continuous = log10(st_qpcr_cm$int_concentation),
	#
	W_st = st_ddpcr$Positives,
	W_en = e_ddpcr$Positives,
	#
	U_st = st_ddpcr$Tot_drop,
	U_en = e_ddpcr$Tot_drop,
	#
	C_ddPCR_st = log10(st_ddpcr$int_concentation),
	#
	i_dst_idx = st_ddpcr$species_idx,
	#
	i_den_idx = e_ddpcr$species_idx,
	ij_den_idx = e_ddpcr$sp_st_idx
);str(stan_data_M2)

initial_values <- replicate(4,list(
	kappa_0 = rep(-11, 3),
	kappa_1 = rep(0.48, 3)
),simplify=F)
```

## The math

$$
\begin{aligned}
&\text{qPCR model}\\
Z_{ijk} &\sim \text{Bernoulli}(\theta_{ij}) && \text{(1)} \\
\text{logit}(\theta_{ij}) &= \phi0_{i} + \phi1_{i} \times \log_{10}(C_{ij}) && \text{(2)} \\
Y_{ijk} &\sim \text{Normal}(\mu_{ij}, \sigma_{ij}) && \text{(3)} \\
\mu_{ij} &= \beta0_{i} + \beta1_{i} \times \log_{10}(C_{ij}) && \text{(4)} \\
\sigma_{ij} &= e^{(\gamma0_{i} + \gamma1_{i} \times \log_{10}(C_{ij}))} && \text{(5)} \\
\\
&\text{ddPCR model}\\
W_{ijk} &\sim \text{Binomial}(\omega_{ij}, U_{ijk}) && \text{(6)} \\
\text{cloglog}(\omega_{ij}) &= \kappa0_{i} + \kappa1_{i} \times \log_{10}(C_{ij}) && \text{(7)}
\end{aligned}
$$ \## The stan code This stan code runs 2 independent models (equation 1 - 5 and 6 - 7) for qPCR and ddPCR data. Both models are run simultaneously for standard and environmental samples where the known concentration from standards informs the intercept and the slope parameters which thereafter informs the initial concentration of unknown samples (C).

```{cpp}
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
```

## Run stan Model

```{r, message=FALSE, results='hide'}
stanMod_M2 = stan(file = here('Code','Stan_qPCR_ddPCR_quantification_precision.stan'),chains = 4,
									model_name='qPCR_vs_ddPCR_quantification_precision.stan',
									iter = 5000,
									warmup = 2000,
									init = initial_values,
									data = stan_data_M2)
```

## Data for Figure 2

```{r, message=FALSE}
# Extract DNA concentration estimation from qPCR model
q <- extract_param(stanMod_M2,"C_qPCR") 

# Extract DNA concentration estimation from ddPCR model
d <- extract_param(stanMod_M2,"C_ddPCR") 

# Combine the data together
diff <- 
	d %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_d","d2","d98")) %>%
	cbind(.,q %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_q","q2","q98"))) %>%
	cbind(.,e_ddpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx))
```

## Plot Figure 2 - Quantification precision (qPCR vs ddPCR)

```{r, message=FALSE,warning=FALSE}
#| fig-width: 18
#| fig-height: 11

source(here('Code','Surpressed','Figure_2.R'))
Figure_2
```

# Two-step model

## Create stan data

```{r}
# Remove an outlier that has an influence on posteriors
e_qpcr_cm <- e_qpcr_cm %>% 
	filter(!(Well=="A4"&
					 	Sample_name=="2021624_10"&
					 	Species =="Herring"&
					 	Ct==42.56510544&
					 	Conc==0.139179796&
					 	plate=="Plate_G"))

stan_data_M3 <- list(
	N_st_q = nrow(st_qpcr),
	N_en_q = nrow(e_qpcr),
	N_st_qp = nrow(st_qpcr_cm),
	N_en_qp = nrow(e_qpcr_cm),
	#
	N_st_d = nrow(st_ddpcr),
	N_en_d = nrow(e_ddpcr),
	#
	N_i = length(unique(st_qpcr$Species)),
	N_ij = length(unique(e_ddpcr$sp_st_idx)),
	#
	i_qst_idx = st_qpcr$species_idx,
	i_qen_idx = e_qpcr$species_idx,
	ij_qen_idx = e_qpcr$sp_st_idx,
	#
	i_qst_p_idx = st_qpcr_cm$species_idx,
	i_qen_p_idx = e_qpcr_cm$species_idx,
	ij_qen_p_idx = e_qpcr_cm$sp_st_idx,
	#
	Z_qPCR_st = st_qpcr$pres,
	Z_qPCR_en = e_qpcr$pres,
	C_qPCR_st = log10(st_qpcr$int_concentation),
	#
	Y_qPCR_st = st_qpcr_cm$Ct,
	Y_qPCR_en = e_qpcr_cm$Ct,
	C_qPCR_st_cont = log10(st_qpcr_cm$int_concentation)
	);str(stan_data_M3)
```

## The math

$$
\begin{aligned}
&\text{qPCR model joint}\\
Z_{ijk} &\sim \text{Bernoulli}(\theta_{ij}) && \text{(1)} \\
\text{logit}(\theta_{ij}) &= \phi0_{i} + \phi1_{i} \times \log_{10}(C_{ij}) && \text{(2)} \\
Y_{ijk} &\sim \text{Normal}(\mu_{ij}, \sigma_{ij}) && \text{(3)} \\
\mu_{ij} &= \beta0_{i} + \beta1_{i} \times \log_{10}(C_{ij}) && \text{(4)} \\
\sigma_{ij} &= e^{(\gamma0_{i} + \gamma1_{i} \times \log_{10}(C_{ij}))} && \text{(5)} \\
\\
&\text{qPCR model continuous (single)}\\
Y_{ijk} &\sim \text{Normal}(\mu_{ij}, \sigma_{ij}) && \text{(3)} \\
\mu_{ij} &= \beta0_{i} + \beta1_{i} \times \log_{10}(C_{ij}) && \text{(4)} \\
\sigma_{ij} &= e^{(\gamma0_{i} + \gamma1_{i} \times \log_{10}(C_{ij}))} && \text{(5)} \\
\end{aligned}
$$

## The stan code

```{cpp}
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
```

## Run stan Model

```{r, message=FALSE, results='hide'}
stanMod_M3 = stan(file = here('Code','Stan_qPCR_two_step_model.stan'),
									chains = 4,
									model_name='qPCR_two_step_model',
									iter = 5000,
									warmup=2000,
									data = stan_data_M3)
```

## Data for Figure 3

```{r, message=FALSE}
q_bin <-
	extract_param(stanMod_M3,"C_qPCR_joint") %>% cbind(.,e_qpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_qPCR_joint","q98","q2","Species")) %>% 
	slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort())
q_nor <-
	extract_param(stanMod_M3,"C_qPCR_continuous") %>% slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort()) %>% 
	cbind(.,e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_qPCR_continuous","q98","q2","Species"))



diff <- cbind(q_bin %>% setNames(c("C_qPCR_joint","q98_bin","q2_bin","Species")) %>% 
								select(-Species),
							q_nor %>% setNames(c("C_qPCR_continuous","q98_nor","q2_nor","Species"))) %>% 
	mutate(diff=(q98_nor-q2_nor)-(q98_bin-q2_bin))
```

## Plot Figure 3 - Two-step qPCR model quantification precision

```{r, message=FALSE,warning=FALSE}
#| fig-width: 20
#| fig-height: 10

source(here('Code','Surpressed','Figure_3.R'))
Figure_3
```

# Lower threshold

## The math

$$
\begin{aligned}
\log_{10}(C_{lt}) = \frac{\text{cloglog}\left(\frac{1}{U}\right) - \kappa_0}{\kappa_1}
\end{aligned}
$$ \## Extract the data from stan Model 2

```{r}
kappa_0 <- extract_param(stanMod_M2,"kappa_0") %>%
	cbind(.,e_ddpcr %>%
					filter(!duplicated(species_idx)) %>%
					arrange(species_idx) %>%
					select(Species,species_idx))
kappa_1 <- extract_param(stanMod_M2,"kappa_1") %>%
	cbind(.,e_ddpcr %>%
					filter(!duplicated(species_idx)) %>%
					arrange(species_idx) %>%
					select(Species,species_idx))
```

```{r}
nrep_st <- 1
nrep_end <- 10
repl <- rep(c(nrep_st:nrep_end),3)
minimum_threshold <- repl %>% as.data.frame() %>%
	cbind(.,kappa_0 %>% pull(mean) %>% rep(.,each=nrep_end)) %>%
	cbind(.,kappa_1 %>% pull(mean) %>% rep(.,each=nrep_end)) %>%
	cbind(.,kappa_0 %>% pull(Species) %>% rep(.,each=nrep_end)) %>%
	setNames(c("replicates","kappa_0","kappa_1","Species")) %>%
	mutate(C_lt=(cloglog(1/(20000*replicates))-kappa_0)/kappa_1) %>%
	mutate(C_ut=(cloglog(((20000*replicates)-1)/(20000*replicates))-kappa_0)/kappa_1)
```

## Plot Figure 4 - ddPCR minimum threshold

```{r, message=FALSE,warning=FALSE}
#| fig-width: 12
#| fig-height: 6

source(here('Code','Surpressed','Figure_4.R'))
Figure_4
```

# Supplementary plots

Load additional data

```{r}
qpcr_optimisation_cod_herring_plate <- read.csv(here('Data','qPCR_opt_Cod_Herring.csv'))
qpcr_optimisation_cod_saithe_plate <- read.csv(here('Data','qPCR_opt_Cod_Saithe.csv'))

ddpcr_amplitude_cod_herring_standard_samples <- read.csv(here('Data','ddpcr_amplitude_cod_herring_standard_samples.csv'))
ddpcr_amplitude_cod_saithe_standard_samples <- read.csv(here('Data','ddpcr_amplitude_cod_saithe_standard_samples.csv'))
```

## Figure S1 - Optimization of qPCR plates

```{r, warning=FALSE, message=FALSE}
opt_dat <- rbind(qpcr_optimisation_cod_saithe_plate,
						 qpcr_optimisation_cod_herring_plate) %>% 
	filter(SuperMix!="OUT") %>% 
	mutate(Conc=if_else(Sample=="ST1",6,0)) %>% 
	mutate(Conc=if_else(Sample=="ST3",4,Conc)) %>% 
	mutate(Conc=if_else(Sample=="ST5",2,Conc)) %>% 
	mutate(Conc=if_else(Sample=="Blank",-Inf,Conc)) %>% 
	mutate(SM=substr(SuperMix,start = 3,stop = 3))

pp1 <-
	opt_dat %>% filter(Conc!=-Inf) %>% 
	ggplot(aes(x=SM,y=Ct,color=Species))+
	geom_boxplot()+
	theme_bw()+
	facet_wrap(~Conc ~ Species,scales="free_y",ncol=3,
						 labeller=label_bquote(cols = 10 ^ .(Conc) ~ "(copies/µL)"))+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	ylab("Cycle threshold (Ct)")+
	xlab("SuperMix (SM)")+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_text(size=17),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=15),
				axis.text.y=element_text(size=15))

p1_leg <-
	opt_dat %>% filter(SuperMix!="OUT"&Conc!=-Inf) %>% 
	ggplot(aes(x=SuperMix,y=Ct,color=Species))+
	geom_boxplot()+
	theme_bw()+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	theme(legend.title = element_text(size = 20),
				legend.text = element_text(size = 16))  # Adjust the size as needed

mat <- matrix(paste0("SM", c(1:3,5:7,9,8,4)), nrow = 3, byrow = TRUE)

df <- expand.grid(row = 1:3, col = 1:3)
df$label <- as.vector(mat)
df$F_p <- rep(c(2.5,10,17),3)
df$R_p <- rep(c(2.5,10,17),each=3)

pp2 <-
	ggplot(df, aes(x = col, y = row, label = label, fill = label)) +
	geom_tile(color = "white", size = 1) +
	geom_text(size = 5, color = "black") +
	geom_text(aes(x = 0.2, y = row, label = F_p), hjust = 0.1,size=5) +
	geom_text(aes(x = col, y = 0.15, label = R_p), vjust = -0.5,size=5) +
	scale_fill_manual(values = rep("grey", 9)) +
	ggtitle("SuperMix (SM) \nconcentration") +  # Add your desired title here
	theme_minimal() +
	ylab("Forward primer \nconcentration (nM)")+
	xlab("Reverse primer \nconcentration (nM)")+
	theme(
		legend.position = "none",
		axis.title.y = element_text(size = 14),
		axis.title.x = element_text(size = 14),
		axis.title = element_blank(),
		plot.margin = margin(1, 0, 0.1, 0.5, "cm"),
		axis.text = element_blank(),
		axis.ticks = element_blank(),
		panel.grid = element_blank(),
		plot.title = element_text(hjust = 0.5,size = 20)
	)

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0.4,0.0)))

p1 <-
	cowplot::plot_grid(p1_leg_leg,pp2, nrow = 3)

Figure_S1 <- cowplot::plot_grid(pp1,p1,rel_widths = c(5,1.5))
```

```{r}
#| fig-width: 16
#| fig-height: 11

Figure_S1
```

## Figure S2 - ddPCR droplet proportion

```{r, warning=FALSE, message=FALSE}
pp1 <-
	ggplot(data=st_ddpcr %>%
				 	group_by(Sample_name, Species) %>%
				 	summarise(nominal_c = mean(int_concentation),
				 						pos_agg = sum(Positives),
				 						tot_agg = sum(Tot_drop),
				 						sp_idx = mean(species_idx),
				 						Conc = mean(Conc))) +
	geom_point(aes(x=Conc, y=cloglog(pos_agg/tot_agg),
								 color=Species),size=2)+
	geom_point(aes(x=nominal_c, y=cloglog(pos_agg/tot_agg),
								 color=Species),shape=23,size=2)+
	geom_smooth(aes(x=Conc, y=cloglog(pos_agg/tot_agg),
									color=Species),
							method=lm, se=F, fullrange=T,lty=3,size=0.5)+
	stat_smooth(aes(x=nominal_c, y=cloglog(pos_agg/tot_agg),
									color=Species),
							method = lm, formula = y ~ poly(x, 3),lty=1,size=0.5, se=F)+
	theme_bw()+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	facet_grid(~Species)+
	scale_x_log10(labels=scientific_10,lim=c(1e-3,1e5),breaks=c(10^seq(-3,5,by=2)))+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0, 0.3, 0, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size = 15),
				axis.text.y=element_text(size=15))+
	ylab(expression(atop("cloglog("*omega*")",paste(omega*" = positive droplets / total droplets"))))+
	xlab("Nominal concentration (copies/µL)")

p1_leg <- ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Gadus morhua"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Clupea harengus"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Pollachius virens"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Gadus morhua"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Clupea harengus"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Pollachius virens"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Assay',
		breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
		values = c("tomato2", "deepskyblue2","orange"),
		guide = guide_legend(
			override.aes = list(lty = c(1,1,1),
													size = c(3,3,3)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,0.0)))

p2_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "ddPCR (in-built estimation)"),size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Nominal concentration"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "ddPCR (in-built estimation)"), size = 1,shape=23) +
	geom_point(aes(x = NA, y = NA, color = "Nominal concentration"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Method of measurement',
		breaks=c("ddPCR (in-built estimation)","Nominal concentration"),
		values = c("black","black"),
		guide = guide_legend(
			override.aes = list(lty = c(3,1),
													size = c(3,3),
													shape= c(19,23)))) +
	theme(plot.margin = margin(1, 0.1, 0, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))
p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,1.0)))

legend <- cowplot::plot_grid(p1_leg_leg, p2_leg_leg, nrow = 2)

Figure_S2 <- cowplot::plot_grid(pp1,legend,nrow=1,rel_widths = c(4,1.1))
```

```{r}
#| fig-width: 16
#| fig-height: 8

Figure_S2
```

## Figure S3 - Possion statistics between assays

```{r, warning=FALSE, message=FALSE}
kappa_0 <- extract_param(stanMod_M2,"kappa_0") %>%
	cbind(.,e_ddpcr %>%
					filter(!duplicated(species_idx)) %>%
					arrange(species_idx) %>%
					select(Species,species_idx))
kappa_1 <- extract_param(stanMod_M2,"kappa_1") %>%
	cbind(.,e_ddpcr %>%
					filter(!duplicated(species_idx)) %>%
					arrange(species_idx) %>%
					select(Species,species_idx))

conc_st <- -2
conc_end <- 5
repl <- rep(seq(conc_st,conc_end,by=0.1),3)
logist_diff <-
	repl %>% as.data.frame() %>%
	cbind(.,kappa_0 %>% pull(mean) %>% rep(.,each=length(repl)/3)) %>%
	cbind(.,kappa_1 %>% pull(mean) %>% rep(.,each=length(repl)/3)) %>%
	cbind(.,kappa_0 %>% pull(Species) %>% rep(.,each=length(repl)/3)) %>%
	setNames(c("C","intercept","slope","Species")) %>%
	mutate(omega=intercept+(slope*C))

sim <- seq(-2,4.9,by=0.1) %>% as.data.frame() %>% setNames("C") %>%
	mutate(omega=inv.cloglog(-7.07+(2.3*C)))

pp1 <-
	logist_diff %>%
	ggplot()+
	geom_line(aes(y=inv.cloglog(omega),x=10^C,color=Species),lwd=1)+
	geom_line(data=sim,aes(y=(omega),x=10^C),color="black",lty=2,lwd=0.3)+
	scale_x_log10(labels=scientific_10,breaks=c(10^seq(conc_st,conc_end)))+ #
	scale_y_sqrt()+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	theme_bw()+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=14),
				axis.text.y=element_text(size=14))+
	ylab(expression("Proportion of positive droplets"))+
	xlab("Nominal concentration (copies/µL)")


p1_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "Gadus morhua"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Clupea harengus"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "Pollachius virens"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Bayesian estimated\nassay',
		breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
		values = c("tomato2", "deepskyblue2","orange"),
		guide = guide_legend(
			override.aes = list(lty = c(1,1,1),
													size = c(2,2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,0.1)))
p2_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "ddPCR (in-built estimation)"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Poisson statistics',
		breaks=c("ddPCR (in-built estimation)"),
		values = c("black"),
		guide = guide_legend(
			override.aes = list(lty = c(2),
													size = c(1)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,1.0)))

legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)
Figure_S3 <- cowplot::plot_grid(pp1,legend,nrow=1,ncol=2,rel_widths = c(4,1.3))
```

```{r}
#| fig-width: 16
#| fig-height: 11

Figure_S3
```

## Figure S4 - Assay performance amplitude

```{r, warning=FALSE, message=FALSE}
amp_dat <- rbind(ddpcr_amplitude_cod_herring_standard_samples,
								 ddpcr_amplitude_cod_saithe_standard_samples)

pp1 <-
	amp_dat %>% 
	ggplot() +
	geom_jitter(aes(y=MeanAmplitudeOfPositives,x=Species,shape = assay),size=2,color="forestgreen")+
	geom_boxplot(aes(x = Species, y = MeanAmplitudeOfPositives),fill="forestgreen",
							 alpha=0.5)+
	geom_jitter(aes(y=MeanAmplitudeOfNegatives,x=Species,shape = assay),size=2,color="tomato2")+
	geom_boxplot(aes(x = Species, y = MeanAmplitudeOfNegatives),fill="tomato2",
							 alpha=0.5)+
	theme_bw()+
	scale_shape_manual(values = c(17, 19)) +
	labs(x = "Species",
			 y = "Mean amplitude of droplets")+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 22),
				axis.title.x = element_text(size = 22),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.text.x = element_text(size=18),
				axis.text.y=element_text(size=16))

p1_leg <-
	ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Positive droplet"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Negative droplet"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Detection',
		breaks=c("Positive droplet", "Negative droplet"),
		values = c("forestgreen", "tomato2"),
		guide = guide_legend(
			override.aes = list(lty = c(1,1),
													size = c(2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))

p2_leg <-
	ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Cod + herring"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Cod + saithe"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="Assay",
		breaks=c("Cod + herring", "Cod + saithe"),
		values = c("black", "black"),
		guide = guide_legend(
			override.aes = list(shape = c(17,19),
													size = c(2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,0.0)))
p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,1.0)))

legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)

Figure_S4 <- cowplot::plot_grid(pp1,legend,ncol = 2, rel_widths = c(4,1))
```

```{r, warning=FALSE}
#| fig-width: 16
#| fig-height: 11

Figure_S4
```

## Figure S5 - LoD and LoQ

```{r, warning=FALSE, error=FALSE, message=FALSE,results='hide'}
# Filter out 2 outliers
DAT <- st_qpcr %>% 
	filter(!(Well=='B2'&Ct==16.80747223)) %>% 
	filter(!(Well=='G1'&Ct==44.88344574))
source(here('Code','Surpressed','Figure_S5_LoD_LoQ.R'))
```

```{r}
#| fig-width: 10
#| fig-height: 20

Figure_S5
```
