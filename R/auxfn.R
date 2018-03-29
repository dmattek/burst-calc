# some defs
myGGtheme = theme_bw(base_size = 18, base_family = "Helvetica") +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major=element_blank(),
        panel.border=element_blank(),
        axis.line = element_line(colour="black"),
        strip.text.x = element_text(size=16, face="bold"),
        strip.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size=16, face='bold'),
        legend.key.height = unit(1, "lines"),
        legend.key.width = unit(1, "lines"),
        legend.position = "right")

#define root finding f-n
f.inv = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}
