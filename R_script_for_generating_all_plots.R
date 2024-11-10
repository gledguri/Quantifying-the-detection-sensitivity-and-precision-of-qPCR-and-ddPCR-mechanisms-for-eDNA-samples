# Load libraries and data -----------------------------------------------------------
rm(list=ls())
library(ggplot2)
library(dplyr)
library(rstan)
library(here);options(mc.cores = parallel::detectCores())

load("All_data_for_reviewers.RData")

# Sensitivity model (Figure 1) is stanMod_M1
# Quantification precision model (Figure 2) is stanMod_M2 
# qPCR two-step model (Figure 3) is stanMod_M3
# For retrieving the stan code of each model type "stanMod_M1@stanmodel"
# For retrieving the arguments of each model type "stanMod_M1@stan_args" or just pressing stanMod_M1 would give a summary of arguments used

# Figure 1 - Detection probability (qPCR vs ddPCR) ----------------------------------
mod_out <- extract_list_param(stanMod_M1)
mod_out <- mod_out[which(names(mod_out) %in% c("alpha_0","alpha_1","beta_0","beta_1"))]


for (i in names(mod_out)) {
	mod_out[[i]] <- mod_out[[i]] %>%
		cbind(.,tibble(sp_idx = 1:3, Species = c("Cod", "Herring", "Saithe")))
}

theta_line_qpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$alpha_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$alpha_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,0))

theta_line_ddpcr <-
	as.data.frame(matrix(rep(seq(-3,6,by=0.1),3),ncol=1)) %>% setNames("C") %>%
	mutate(Species=rep(c("Cod","Herring","Saithe"),each=nrow(.)/3)) %>%
	left_join(.,mod_out$beta_0 %>% select(Species,mean),by="Species") %>%
	rename(intercept=mean) %>%
	left_join(.,mod_out$beta_1 %>% select(Species,mean),by="Species") %>%
	rename(slope=mean) %>%
	mutate(prob=logreg(C,intercept,slope,2))

theta_line <- cbind(theta_line_qpcr %>% rename(prob_q=prob),
										theta_line_ddpcr %>% select(prob) %>% rename(prob_dd=prob)) %>%
	mutate(delta=prob_dd-prob_q)

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
	# scale_x_log10(labels=scientific_10,lim=c(1e-4,1e6))+
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

legend <- cowplot::plot_grid(p1_leg_leg,p3_leg_leg,p2_leg_leg,p4_leg_leg,nrow = 4)

p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,rel_heights = c(4.7,5.3), labels = c("a","b"),label_size = 20, align = "v")
Figure_1 <- cowplot::plot_grid(p1,legend,nrow = 1,ncol = 2,rel_widths = c(7,1.8))
Figure_1


# Figure 2 - Quantification precision (qPCR vs ddPCR) -------------------------------
q <- extract_param(stanMod_M2,"C_q")
d <- extract_param(stanMod_M2,"C_d")
de <- extract_param(stanMod_M2,"delta")

diff <- de %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("delta","de2","de98")) %>%
	cbind(.,d %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_d","d2","d98"))) %>%
	cbind(.,q %>% select(mean,`2.5%`,`97.5%`) %>% setNames(c("C_q","q2","q98"))) %>%
	cbind(.,e_ddpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx))

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

legend <- cowplot::plot_grid(p1_leg_leg,p2_leg_leg,nrow = 2)
p1 <- cowplot::plot_grid(pp1,pp2,nrow=2,align = "v",rel_heights = c(4,3),label_size = 20,labels = c("a","b"),label_x = -0.003)
Figure_2 <- cowplot::plot_grid(p1,legend,ncol=2,rel_widths = c(4,1.0))
Figure_2


# Figure 3 - Two-step qPCR model quantification precision ---------------------------
q_bin <-
	extract_param(stanMod_v4,"C_q") %>% cbind(.,e_qpcr %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_bin","q98","q2","Species")) %>% 
	slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort())
q_nor <-
	extract_param(stanMod_v4,"C_nor") %>% slice(e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% pull(sp_st_idx) %>% sort()) %>% 
	cbind(.,e_qpcr_cm %>% filter(!duplicated(sp_st_idx)) %>% arrange(sp_st_idx) %>% select(-Ct,-Conc,-pres)) %>% 
	select(mean,`97.5%`,`2.5%`,Species) %>% setNames(c("C_nor","q98","q2","Species"))



diff <- cbind(q_bin %>% setNames(c("C_bin","q98_bin","q2_bin","Species")) %>% 
								select(-Species),
							q_nor %>% setNames(c("C_nor","q98_nor","q2_nor","Species"))) %>% 
	mutate(diff=(q98_nor-q2_nor)-(q98_bin-q2_bin))

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

p1 <- cowplot::plot_grid(pp1,pp2,nrow = 2,align = "v",rel_heights = c(3.5,4),labels = c("a","b"),label_size = 20, label_x = 0.08)

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


p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,-1.6)))
p2_leg_leg <- cowplot::get_legend(p2_leg+
																		theme(legend.justification = c(0,-0.4)))
p3_leg_leg <- cowplot::get_legend(p3_leg+
																		theme(legend.justification = c(0,0.8)))

legend <-
	cowplot::plot_grid(p1_leg_leg,
										 p2_leg_leg,
										 p3_leg_leg,
										 nrow = 4)
Figure_3 <- cowplot::plot_grid(p1,legend,ncol=2,rel_widths = c(4,0.8))
Figure_3


# Figure 4 - ddPCR minimum threshold ------------------------------------------------
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

nrep_st <- 1
nrep_end <- 10
repl <- rep(c(nrep_st:nrep_end),3)
minimum_threshold <- repl %>% as.data.frame() %>%
	cbind(.,kappa_0 %>% pull(mean) %>% rep(.,each=nrep_end)) %>%
	cbind(.,kappa_1 %>% pull(mean) %>% rep(.,each=nrep_end)) %>%
	cbind(.,kappa_0 %>% pull(Species) %>% rep(.,each=nrep_end)) %>%
	setNames(c("replicates","intercept","slope","Species")) %>%
	mutate(C_th=(cloglog(1/(20000*replicates))-intercept)/slope) %>%
	mutate(C_th_up=(cloglog(((20000*replicates)-1)/(20000*replicates))-intercept)/slope)

pp1 <-
	minimum_threshold %>%
	ggplot()+
	geom_line(aes(y=C_th,x=replicates,color=Species))+
	geom_point(aes(y=C_th,x=replicates,color=Species))+
	scale_x_continuous(breaks=c(1,3,5,7,10))+
	scale_y_continuous(breaks=seq(-2.5,-1,by=0.2))+
	scale_color_manual(values=c("tomato2", "deepskyblue2","orange2"))+
	theme_bw()+
	facet_grid(~Species)+
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.title.y = element_text(size = 19),
				plot.margin = margin(0.1, 0.1, 0.5, 0.5, "cm"),
				axis.title.x = element_text(size = 19),
				axis.text.x = element_text(size=14),
				axis.text.y=element_text(size=14))+
	ylab(expression(atop("Lower threshold",paste("Log"[10]*" concentration (copies/µL)"))))+
	xlab("Number of replicates")

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
													size = c(2,2,2)))) +
	theme(plot.margin = margin(1, 4, 1, 0, "cm"),
				legend.key.width = unit(1.2,"cm"),
				legend.title = element_text(size=18),
				legend.key.height = unit(0.8, 'cm'),
				legend.text = element_text(size=16,face = "italic"))

p1_leg_leg <- cowplot::get_legend(p1_leg+
																		theme(legend.justification = c(0,0.5)))


Figure_4 <- cowplot::plot_grid(pp1,p1_leg_leg,nrow=1,ncol=2,rel_widths = c(4,1))
Figure_4

# Figure S1 - Optimization of qPCR plates -------------------------------------------
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
Figure_S1

# Figure S2 - ddPCR droplet proportion ----------------------------------------------
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
Figure_S2

# Figure S3 - Possion statistics between assays ------------------------------------
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
Figure_S3


# Figure S4 - Assay performance amplitude -------------------------------------------
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
Figure_S4

