# fy1--surface potential of particle 
# fy2--surface potential of collector surface
# r--radius of particle
# A-- Hamaker constant
library(readxl)
library(gstat)
library(ggplot2)
library(stringr)
library(ggpmisc)
library(scales) 
library(ggpubr)
library(multcompView)
library(car)

setwd('/Users/teagank/Desktop/WSU/NSF/Column\ Experiments')

l<-seq(0.001,100,0.1)
dlvo_NaCl <-function(fy1,fy2,r,A){
  
  vdw <- (- A * r / ( 6 * l * 10 ^ -9 ) ) / (( 1 + 14 * l * 10 ^ -9 / (100 * 10 ^ -9) )) 
  # vdw <- (- A * r / ( 6 * l * 10 ^ -9 ) ) * (1 - ((5.32 * l * 10 ^ -9)/ (100 * 10 ^ -9) * log(1 + (100 * 10 ^ -9) / (l * 10 ^ -9) )))
  ee = 80.1 * 8.854 * 10 ^ -12
  temp =  293
  I = 10              
  e = 1.6 * 10 ^ -19
  k = 1.3806503 * 10 ^ -23
  kt = k * temp
  na = 6.02214179 * 10 ^ 23
  kapa = sqrt((2 * na * I * e ^ 2) / (ee * k * temp))
  ele <- pi * r * ee * (2 * fy1 * fy2 * log((1 + exp(-kapa * l * 10 ^ -9)) / (1 - exp(-kapa * l * 10 ^ -9))) + (fy1 ^ 2 + fy2 ^ 2) * log((1 - exp(-2 * kapa * l * 10 ^ -9))))
  
  Classical_DLVO <- vdw / kt + ele / kt
  Classical_DLVO
}

PS_COOH_sand = dlvo_NaCl(fy1 = -57.86875 * (10 ^ -3), fy2 = -45.2 * 10 ^ -3,r = 221.4166667 / 2 * 10 ^ -9, 4.04 * 10 ^ -21)
PS_COOH_sand_BSA = dlvo_NaCl(fy1 = -40 * (10 ^ -3), fy2 = -45.2 * 10 ^ -3,r = 221.4166667 / 2 * 10 ^ -9, 4.04 * 10 ^ -21)
PS_COOH_sand_LSZ = dlvo_NaCl(fy1 = -30 * (10 ^ -3), fy2 = -45.2 * 10 ^ -3,r = 221.4166667 / 2 * 10 ^ -9, 4.04 * 10 ^ -21)
PS_COOH_sand_DOM = dlvo_NaCl(fy1 = -20 * (10 ^ -3), fy2 = -45.2 * 10 ^ -3,r = 221.4166667 / 2 * 10 ^ -9, 4.04 * 10 ^ -21)




data_PS_COOH_sand = cbind.data.frame(l, PS_COOH_sand, "PS-COOH", "1")
data_PS_COOH_LSZ_sand = cbind.data.frame(l, PS_COOH_sand_LSZ, "PS-COOH + LSZ", "2")
data_PS_COOH_BSA_sand = cbind.data.frame(l, PS_COOH_sand_BSA, "PS-COOH + BSA", "3")
data_PS_COOH_DOM_sand = cbind.data.frame(l, PS_COOH_sand_DOM, "PS-COOH + Compost Extract", "4")

colnames(data_PS_COOH_sand) = c("l","Classical_DLVO", "type", "order")
colnames(data_PS_COOH_LSZ_sand) = c("l","Classical_DLVO", "type", "order")
colnames(data_PS_COOH_BSA_sand) = c("l","Classical_DLVO", "type", "order")
colnames(data_PS_COOH_DOM_sand) = c("l","Classical_DLVO", "type", "order")

data_PS_COOH = rbind.data.frame(data_PS_COOH_sand,data_PS_COOH_LSZ_sand,data_PS_COOH_BSA_sand,data_PS_COOH_DOM_sand)


pdf("PS_COOH_DLVO_classical.pdf",width = 4,height = 3)

ggplot(data_PS_COOH, aes(x = l, y = Classical_DLVO, colour = reorder(type, order), linetype = reorder(type, order))) +
  geom_line(size = 0.5) + 
  scale_color_manual(values= rep(c("red","black","#FF00FF","blue"), 2)) +
  scale_x_continuous(expand = c(0,0), limits = c(-0.0010,50), name = "Separation Distance (nm)",  breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(expand = c(0,0),limits = c(-20,400),name =  "Total Interaction Energy (kT)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") + 
  #annotate(geom="text", x=2, y=280, label="A") +
  ggtitle("Classical DLVO") + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0),
    plot.margin=unit(c(0.2,0.2,0.2,0.2),"cm"),
    panel.border = element_rect(fill=NA, colour = "black", size=1),
    legend.title = element_blank(),
    legend.position = c(0.6,0.75),
    legend.background = element_rect(fill=alpha('black', 0)),
    #legend.key.size = unit(1.5, 'cm'),
    legend.spacing.y = unit(0.1, 'cm'),
    legend.spacing.x = unit(0.2, 'cm'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size= 12,color="black",face="plain"),
    axis.title.y = element_text(size= 12,color="black",face="plain"),
    axis.text.x = element_text(size= 10,color="black",face="plain"),
    axis.text.y = element_text(size= 10,color="black",face="plain")  
  )
dev.off()
