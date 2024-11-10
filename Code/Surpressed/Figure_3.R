# Partial plot 1
pp1 <-
	diff %>%
	ggplot()+
	geom_point(aes(y=(q98_bin-q2_bin),x=10^C_qPCR_joint,color=Species),size=2,shape=16)+
	geom_point(aes(y=(q98_nor-q2_nor),x=10^C_qPCR_continuous,color=Species),size=2,shape=3)+
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
	geom_smooth(aes(y=diff,x=10^C_qPCR_joint,color=Species),lty=2,se=F)+
	geom_point(aes(y=diff,x=10^C_qPCR_joint,color=Species),size=2,shape=4)+
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