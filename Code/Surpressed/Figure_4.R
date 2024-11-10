pp1 <-
	minimum_threshold %>%
	ggplot()+
	geom_line(aes(y=C_lt,x=replicates,color=Species))+
	geom_point(aes(y=C_lt,x=replicates,color=Species))+
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
