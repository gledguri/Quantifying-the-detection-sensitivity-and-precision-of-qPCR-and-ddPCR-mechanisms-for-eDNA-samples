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