# Functions -----------------------------------------------------------------------------------
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


# Load data -----------------------------------------------------------------------------------

e_ddpcr <- read.csv(here('Data','ddPCR_environmental_samples.csv'))
e_qpcr <- read.csv(here('Data','qPCR_environmental_samples.csv'))
st_ddpcr <- read.csv(here('Data','ddPCR_standard_samples.csv'))
st_qpcr <- read.csv(here('Data','qPCR_standard_samples.csv'))

st_qpcr_cm <- st_qpcr %>% filter(!is.na(Ct))
e_qpcr_cm <- e_qpcr %>% filter(!is.na(Ct))

sample_metadata <- read.csv(here('Data','Sample_metadata.csv'))
stidx <- read.csv(here('Data','sample_idx.csv'))

# qpcr_optimisation_cod_herring_plate <- read.csv(here('Data','qPCR_opt_Cod_Herring.csv'))
# qpcr_optimisation_cod_saithe_plate <- read.csv(here('Data','qPCR_opt_Cod_Saithe.csv'))

# Stan Model - probability of detection -------------------------------------------------------
stan_data_M1 <- list(
	N_i = length(unique(st_qpcr$Species)),

	Nq = nrow(st_qpcr),
	Ndd = nrow(st_ddpcr),

	i_q_idx = st_qpcr$species_idx,
	i_d_idx = st_ddpcr$species_idx,

	Z_q = st_qpcr$pres,
	Z_dd = st_ddpcr$pres,

	S_q = log10(st_qpcr$int_concentation),
	S_dd = log10(st_ddpcr$int_concentation)
);str(stan_data_M1)

initial_values <- replicate(4,list(
	alpha_0 = rep(1.0, 3),
	alpha_1 = rep(0.4, 3),
	beta_0 = rep(2.50, 3),
	beta_1 = rep(1.23, 3)
),simplify=F)

stanMod_v1 = stan(file = here('Code','Stan_qPCR_ddPCR_prob_of_det.stan'), 
									chains = 4, 
									model_name = 'qPCR_vs_ddPCR_probability_of_detection', 
									iter = 10000, 
									warmup = 6000, 
									init = initial_values, 
									data = stan_data_M1)

# Figure 1 - Detection probability (qPCR vs ddPCR) ----------------------------------
# Extract parameters to build the graph
mod_out <- extract_list_param(stanMod_v1)
mod_out <- mod_out[which(names(mod_out) %in% c("alpha_0","alpha_1","beta_0","beta_1"))]

# Add species names
for (i in names(mod_out)) {
	mod_out[[i]] <- mod_out[[i]] %>%
		cbind(.,tibble(sp_idx = 1:3, Species = c("Cod", "Herring", "Saithe")))
}

# Construct the qPCR probability of detection line from alpha_0 & alpha_1 parameters
theta_line_qpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$alpha_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$alpha_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,0))

# Construct the ddPCR probability of detection line from beta_0 & beta_1 parameters
theta_line_ddpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$beta_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$beta_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,2)) #added 2 to logistic regression (see why in stan model file)

# Combine the logistic regresion lines (prob of det) of both data (qPCR and ddPCR)
theta_line <- cbind(theta_line_qpcr %>% rename(prob_q=prob),
										theta_line_ddpcr %>% select(prob) %>% rename(prob_dd=prob)) %>%
	mutate(delta=prob_dd-prob_q)

# Partial plot 1
pp1 <-
	theta_line %>%
	ggplot()+
	geom_line(aes(x=10^C,y=prob_dd,color=Species),lty=2,size=0.9)+
	geom_line(aes(x=10^C,y=prob_q,color=Species),lty=1,size=0.9)+
	geom_point(data=st_qpcr %>% group_by(Species,Sample_name) %>%
						 	summarise(int_concentation=mean(int_concentation),
						 						pres=mean(pres)),
						 aes(x=int_concentation,y=pres,color=Species),shape=15,size=3)+
	geom_point(data=st_ddpcr %>% group_by(Species,Sample_name) %>%
						 	summarise(int_concentation=mean(int_concentation),
						 						pres=mean(pres)),
						 aes(x=int_concentation,y=pres,color=Species),shape=19,size=3)+
	scale_x_log10(labels=NULL,lim=c(1e-2,1e6),breaks=10^c(-2,-1,0,1,2,4,6))+
	ylim(0,1)+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	facet_grid(~Species)+
	theme_bw()+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, -0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=15),
				axis.text.y=element_text(size=15))+
	ylab("Detection probability")+
	xlab("")

# Partial plot 2
pp2 <-
	theta_line %>%
	ggplot()+
	geom_line(aes(x=10^C,y=delta,color=Species),lty=6,size=0.9)+
	scale_x_log10(labels=scientific_10,lim=c(1e-2,1e6),breaks=10^c(-2,-1,0,1,2,4,6))+
	ylim(-0.1,0.8)+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	facet_grid(~Species)+
	theme_bw()+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=15),
				axis.text.y=element_text(size=15))+
	ylab("ddPCR - qPCR \ndetection probability")+
	xlab("Nominal DNA concentration (copies/µL)")

# Legend 1
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
																		theme(legend.justification = c(0,-1.3)))

# Legend 2
p2_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "Sensitivity difference\n(ddPCR - qPCR)"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="",
		breaks=c("Sensitivity difference\n(ddPCR - qPCR)"),
		values = c("black"),
		guide = guide_legend(
			override.aes = list(lty = c(6),
													size = c(1)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=1),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))
p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,0.95)))

# Legend 3
p3_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "ddPCR sensitivity"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "qPCR sensitivity"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="Modelled sensitivity",
		breaks=c("ddPCR sensitivity","qPCR sensitivity"),
		values = c("black","black"),
		guide = guide_legend(
			override.aes = list(lty = c(2,1),
													size = c(1,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))
p3_leg_leg <- cowplot::get_legend(p3_leg+
																		theme(legend.justification = c(0,-0.1)))

# Legend 4
p4_leg <-
	ggplot() +
	geom_point(aes(x = NA, y = NA, color = "ddPCR standard samples"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "qPCR standard samples"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="Measured sensitivity",
		breaks=c("ddPCR standard samples","qPCR standard samples"),
		values = c("black","black"),
		guide = guide_legend(
			override.aes = list(shape = c(19,15),
													size = c(3)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))
p4_leg_leg <- cowplot::get_legend(p4_leg+
																		theme(legend.justification = c(0,2.2)))

# Combine legends
legend <- cowplot::plot_grid(p1_leg_leg,p3_leg_leg,p2_leg_leg,p4_leg_leg,nrow = 4)

# Combine partial plots
p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,rel_heights = c(4.7,5.3), labels = c("a","b"),label_size = 20, align = "v")

# Combine partial plots and legends
Figure_1 <- cowplot::plot_grid(p1,legend,nrow = 1,ncol = 2,rel_widths = c(7,1.8))
Figure_1


# Stan Model - quantification precision -------------------------------------------------------
# Create stan data
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
	
	Z_qst = st_qpcr$pres,
	Z_qen = e_qpcr$pres,
	S_q = log10(st_qpcr$int_concentation),
	
	R_qst = st_qpcr_cm$Ct,
	R_qen = e_qpcr_cm$Ct,
	S_q_p = log10(st_qpcr_cm$int_concentation),
	#
	W_st = st_ddpcr$Positives,
	W_en = e_ddpcr$Positives,
	#
	U_st = st_ddpcr$Tot_drop,
	U_en = e_ddpcr$Tot_drop,
	#
	S_d = log10(st_ddpcr$int_concentation),
	#
	i_dst_idx = st_ddpcr$species_idx,
	#
	i_den_idx = e_ddpcr$species_idx,
	ij_den_idx = e_ddpcr$sp_st_idx
);str(stan_data_M2)

initial_values <- replicate(4,list(
	kappa_0 = rep(-11, 3),
	kappa_1 = rep(0.48, 3),
	kappa_2 = rep(7.93, 3)
),simplify=F)

# Run stan model
stanMod_M2 = stan(file = here('Code','Stan_qPCR_ddPCR_quantification_precision.stan'),chains = 4,
									model_name='qPCR_vs_ddPCR_quantification_precision.stan',
									iter = 5000,
									warmup=2000,
									init = initial_values,
									data = stan_data_M2)

# Figure 2 - Quantification precision (qPCR vs ddPCR) -------------------------------

# Extract DNA concentration estimation from qPCR model
q <- extract_param(stanMod_M2,"C_q") 

# Extract DNA concentration estimation from ddPCR model
d <- extract_param(stanMod_M2,"C_d") 

# Combine the data together
diff <- 
	d %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_d","d2","d98")) %>%
	cbind(.,q %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_q","q2","q98"))) %>%
	cbind(.,e_ddpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx))

# Partial plot 1
pp1 <-
	diff %>%
	ggplot()+
	geom_errorbar(aes(y=10^C_d,x=10^C_q,ymin=10^d2,ymax=10^d98),color="grey")+
	geom_errorbar(aes(y=10^C_d,x=10^C_q,xmin=10^q2,xmax=10^q98),color="grey")+
	geom_point(aes(y=10^C_d,x=10^C_q,color=Species))+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	geom_abline(intercept = 0,slope=1,lty=2)+
	theme_bw()+
	scale_y_log10(labels=scientific_10,lim=c(1e-4,1e4),breaks=10^c(-4,-2,0,2,4))+
	scale_x_log10(labels=scientific_10,lim=c(1e-5,2e2),breaks=10^c(-4,-2,0,2))+
	facet_grid(~Species)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=18),
				axis.text.y=element_text(size=14))+
	ylab("ddPCR modelled concentration \n(copies/µL)")+
	xlab("qPCR modelled concentration (copies/µL)")

# Partial plot 2
pp2 <-
	diff %>% filter(C_d>-2) %>%
	filter(C_q>-2) %>%
	ggplot()+
	geom_smooth(aes(y=10^(d98-d2),x=10^C_d,color=Species),lty=2,se=T)+
	geom_smooth(aes(y=10^(q98-q2),x=10^C_q,color=Species),lty=1,se=T)+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	scale_y_log10(labels=scientific_10,lim=c(1e0,1e3))+
	scale_x_log10(labels=scientific_10,lim=c(1e-5,1e2),breaks=10^c(-4,-2,0,2))+
	theme_bw()+
	facet_grid(~Species)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=18),
				axis.text.y=element_text(size=14))+
	ylab(expression("95% Credible Interval range"))+
	xlab("Modelled concentration (copies/µL)")

# Legend 1
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
			override.aes = list(#lty = c(1,1,1),
				size = c(3,3,3)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,0.1)))

# Legend 2
p2_leg <- ggplot() +
	geom_line(aes(x = NA, y = NA, color = "ddPCR"), size = 1) +
	geom_line(aes(x = NA, y = NA, color = "qPCR"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="Quantification precision \n(expressed as variance)",
		breaks=c("ddPCR","qPCR"),
		values = c("black","black"),
		guide = guide_legend(
			override.aes = list(lty = c(2,1),
													size = c(1,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))
p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,1.0)))

# Combine legend
legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)

# Combine partial plots
p1 <- cowplot::plot_grid(pp1,pp2,nrow=2,align = "v",rel_heights = c(4,3),label_size = 20,labels = c("a","b"),label_x = -0.003)

# Combine partial plots and legend
Figure_2 <- cowplot::plot_grid(p1,legend,ncol=2,rel_widths = c(4,1.0))
Figure_2

# Stan Model - two step model -----------------------------------------------------------------

# Remove an outlier that has an influence on posteriors
e_qpcr_cm <- e_qpcr_cm %>% 
	filter(!(Well=="A4"&
					 	Sample_name=="2021624_10"&
					 	Species =="Herring"&
					 	Ct==42.56510544&
					 	Conc==0.139179796&
					 	plate=="Plate_G"))

# Create stan data
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
	Z_qst = st_qpcr$pres,
	Z_qen = e_qpcr$pres,
	S_q = log10(st_qpcr$int_concentation),
	#
	R_qst = st_qpcr_cm$Ct,
	R_qen = e_qpcr_cm$Ct,
	S_q_p = log10(st_qpcr_cm$int_concentation)
	);str(stan_data_M3)

# Run stan model
stanMod_M3 = stan(file = here('Code','Stan_qPCR_two_step_model.stan'),
									chains = 4,
									model_name='qPCR_two_step_model',
									iter = 5000,
									warmup=2000,
									data = stan_data_M3)

# Figure 3 - Two-step qPCR model quantification precision ---------------------------
q_bin <-
	extract_param(stanMod_M3,"C_q") %>% cbind(.,e_qpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_bin","q98","q2","Species")) %>% 
	slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort())
q_nor <-
	extract_param(stanMod_M3,"C_nor") %>% slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort()) %>% 
	cbind(.,e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_nor","q98","q2","Species"))



diff <- cbind(q_bin %>% setNames(c("C_bin","q98_bin","q2_bin","Species")) %>% 
								select(-Species),
							q_nor %>% setNames(c("C_nor","q98_nor","q2_nor","Species"))) %>% 
	mutate(diff=(q98_nor-q2_nor)-(q98_bin-q2_bin))

# Partial plot 1
pp1 <-
	diff %>%
	ggplot()+
	geom_point(aes(y=(q98_bin-q2_bin),x=10^C_bin,color=Species),size=2,shape=16)+
	geom_point(aes(y=(q98_nor-q2_nor),x=10^C_nor,color=Species),size=2,shape=3)+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	scale_x_log10()+
	theme_bw()+
	facet_grid(~Species)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 1, "cm"),
				axis.title.x = element_blank(),
				axis.text.x = element_blank(),
				axis.text.y=element_text(size=14))+
	ylab(expression(atop("Quantification variance \n(measured as 95% CI)",paste("log"[10]*"(copies/µL)"))))

# Partial plot 2
pp2 <-
	diff %>%
	ggplot()+
	geom_smooth(aes(y=diff,x=10^C_bin,color=Species),lty=2,se=F)+
	geom_point(aes(y=diff,x=10^C_bin,color=Species),size=2,shape=4)+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	theme_bw()+
	scale_x_log10(labels=scientific_10)+
	facet_grid(~Species)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 1, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=18),
				axis.text.y=element_text(size=14))+
	ylab(expression(atop("Quantification variance difference \n     (measured as CI difference)",paste("(CI"[continuous]*" - CI"[logistic*" + "* continuous]*")"))))+
	ylim(-0.3,0.8)+
	xlab("Modelled qPCR concentration (copies/µL)")

# Legend 1
p1_leg <- ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Gadus morhua"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Clupea harengus"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Pollachius virens"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Species',
		breaks=c("Gadus morhua", "Clupea harengus","Pollachius virens"),
		values = c("tomato2", "deepskyblue2","orange"),
		guide = guide_legend(
			override.aes = list(#lty = c(1,1,1),
				size = c(2,2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,-1.6)))

# Legend 2
p2_leg <- ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Joint model"), size = 1) +
	geom_point(aes(x = NA, y = NA, color = "Continuous model"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name='Model used',
		breaks=c("Joint model", "Continuous model"),
		values = c("black", "black"),
		guide = guide_legend(
			override.aes = list(shape = c(16,3),
													size = c(2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))

p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,-0.4)))

# Legend 3
p3_leg <- ggplot() +
	geom_point(aes(x = NA, y = NA, color = "Precision difference"), size = 1) +
	theme_bw()+
	scale_color_manual(
		name="",
		breaks=c("Precision difference"),
		values = c("black"),
		guide = guide_legend(
			override.aes = list(shape = c(4),
													size = c(2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=1),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16))

p3_leg_leg <- cowplot::get_legend(p3_leg+
																		theme(legend.justification = c(0,0.8)))

# Combine legends
legend <-
	cowplot::plot_grid(p1_leg_leg,
										 p2_leg_leg,
										 p3_leg_leg,
										 nrow = 4)

# Combine partial plots
p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,align = "v",rel_heights = c(3.5,4),labels = c("a","b"),label_size = 20, label_x = 0.08)

# Combine partial plots and legend
Figure_3 <- cowplot::plot_grid(p1,legend,ncol=2,rel_widths = c(4,0.8))
Figure_3
