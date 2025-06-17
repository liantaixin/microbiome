

KEGG <- read.table("L2.txt",header = TRUE,sep = "\t")

library(ggplot2)

library(reshape2)


swr = function(string, nwrap = 12){paste(strwrap(string,width = nwrap),collapse = "\n")}


swr = Vectorize(swr)

KEGG$L1 <- swr(KEGG$L1)

p <- ggplot(KEGG,aes(Abundance,L2)) +
  geom_bar(aes(fill = L1),stat = "identity",width = 0.6) +
  xlab("Relative abundance (%)") +
  ylab("KEGG Pathway") +
  theme(panel.background = element_rect(fill = "white",colour='black'),
        panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
        panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=20,face = "bold"),
        axis.title.y=element_text(colour='black', size=20),
        axis.text.x=element_text(colour='black',size=20),
        axis.text.y = element_text(color = "black",size = 12),
        legend.position = "none",
        strip.text.y = element_text(angle = 0,size = 12,face = "bold")) +
  facet_grid(L1~.,space = "free_y",scales = "free_y")


png(filename="L2.png",res=600,height=7200,width=6000)
p
dev.off()