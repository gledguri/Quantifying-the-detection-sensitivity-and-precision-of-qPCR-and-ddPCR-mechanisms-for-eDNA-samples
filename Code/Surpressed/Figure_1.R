# Partial plot 1
pp1 <-
	probability_of_detection_data %>%
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
	probability_of_detection_data %>%
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
