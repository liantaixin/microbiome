
rm(list=ls())



library(ggplot2)
library(ggExtra)
library(vegan)
library(ggthemes)
library(reshape2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(gghalves)
library(reshape2)


data<-read.csv("D:/phD/16s/16s_files/β-diversity/β2.csv", header=TRUE,row.names=1)
head(data)

data<-t(data)


group<-read.csv('D:/phD/16s/16s_files/β-多样性/group2.csv',header=TRUE,row.names=1,stringsAsFactors=FALSE)

bray<-vegdist(data,              
              method='bray',#"manhattan","euclidean","bray","jaccard",
                            #"clark", "kulczynski","gower", "altGower", "morisita","canberra",
                            #"horn","mountford","raup","binomial","chao","cao","mahalanobis",
                            #"chisq","chord","hellinger","aitchison","robust.aitchison".
              binary=T,## Standardization. Since I was using absolute abundance, I adopted T. Relative abundance was replaced by F              
              diag=FALSE, #Calculate the diagonal              
              upper=FALSE,#Only return to the upper diagonal              
              na.rm = FALSE#If there are missing values, fill in T              
              )

#jaccard distance
jaccard<-vegdist(data,method="jaccard",binary=T,diag=FALSE,upper=FALSE,na.rm = FALSE)
#euclidean distance
euclidean<-vegdist(data,method="euclidean",binary=T,diag=FALSE,upper=FALSE,na.rm = FALSE)
#Take the Brecurtis distance matrix as an example to sort out the data
bray<-as.matrix(bray)
bray

#PCoA
pcoa<-cmdscale(bray,k=2,eig=T)

pcoa_data<-data.frame({pcoa$point})
pcoa_data$Sample_ID<-rownames(pcoa_data)
names(pcoa_data)[1:3]<-paste0("PCoA",1:3)
eig=pcoa$eig
eig_percent<-round(pcoa$eig/sum(pcoa$eig)*100,3)
eig_percent
poi=pcoa$points
poi=as.data.frame(poi)

#Add grouping information for the sample point coordinates
pcoa_data

pcoa_result<-cbind(pcoa_data,group)#PCoA values and tax table/ASV table are merged


head(pcoa_result)

#Conduct permutation multivariate (factor) analysis of variance (PERMANOVA)
#Based on the Bray-Curtis distance

dune.div<-adonis2(data~group,data=group,permutations=999,method="bray")
dune.div

dune_adonis<-paste0("adonis R2:",round(dune.div$R2,2),";P-value:",dune.div$`Pr(>F)`)

#plot

p=ggplot(pcoa_result,aes(x=PCoA1,y=PCoA2,color=group))+
  geom_point(aes(color=group),size=4)+
  labs(x=paste("PCoA1(",eig_percent[1],"%)",sep=""),       
       y=paste("PCoA2(",eig_percent[2],"%)",sep=""),
       caption=dune_adonis)+
  scale_colour_manual(values=c("#FF0000","#800080","#FFA500","#0000FF"))+  
  stat_ellipse(data=pcoa_result,geom="polygon",level=0.9,                
               linetype=2,linewidth=0.5,aes(fill=group),                
               alpha=0.3,show.legend=T)+
  scale_fill_manual(values=c("#FF0000","#800080","#FFA500","#0000FF"))+
  theme(legend.position=c(0.12,0.85),
        legend.title=element_blank(),
        panel.grid=element_blank(),plot.title=element_text(hjust=0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_rect(color='black',fill='transparent'),
        axis.text=element_text(color="black",size=10))+
  geom_hline(aes(yintercept=0),colour="#BEBEBE",linetype="dashed")+
  geom_vline(aes(xintercept=0),colour="#BEBEBE",linetype="dashed")
p


##NMDS
nmds<-metaMDS(bray, k = 2)

stress.score <- nmds$stress

stress.score
nmds <- data.frame(nmds$point)
nmds
nmdsdata <- cbind(nmds, group)
nmdsdata

p1<-ggplot(data=nmdsdata,aes(x=MDS1,y=MDS2))+#Specify data, X-axis, Y-axis, and color
  geom_point(aes(color = group), shape = 19, size=4)+#Draw a dot plot and set its size
  geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
  geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+#The dotted line in the picture
  theme_bw()+#theme settings
  stat_ellipse(data=nmdsdata,#Add the confidence ellipse
               geom = "polygon",level=0.95,#What confidence interval is displayed
               linetype = 2,size=0.5,#Coil type size
               aes(fill=group),
               alpha=0.2)+            
  scale_color_manual(values = c("#FF0000","#800080","#FFA500","#0000FF")) +#dot
  scale_fill_manual(values = c("#FF0000","#800080","#FFA500","#0000FF"))+#ellipse
  theme(axis.title.x=element_text(size=12),#Modify the title text of the X-axis
        axis.title.y=element_text(size=12,angle=90),#Modify the Y-axis title text
        axis.text.y=element_text(size=10),#Modify the text of the X-axis scale label
        axis.text.x=element_text(size=10),#Modify the text of the Y-axis scale label
        panel.grid=element_blank())+#Hidden grid lines
  ggtitle(paste('Stress=',round(stress.score, 3)))#Add the value of the stress function
p1


#box plots
be <- as.data.frame(bray)
head(be)
data <- cbind(be, group)
head(data)

beta <- melt(data, 
             id.vars = "group",
             variable.name = "Comparison",
             value.name = "Beta"
              )
head(beta)


p2<-ggplot(beta,aes(group,Beta,fill=group))+
  geom_half_violin(position = position_nudge(x=0.25),
                   side = "r",width=0.6,color=NA)+
  geom_boxplot(width=0.4,size=1,outlier.color =NA)+
  geom_jitter(aes(fill=group),shape=21,size=2.5,width=0.2)+
  geom_signif(comparisons = list(c("Content","Hemolymph"),
                                 c("Content","Sediment"),
                                 c("Content","Water"),
                                 c("Hemolymph","Sediment"),
                                 c("Hemolymph","Water"),
                                 c("Sediment","Water")),
              map_signif_level = TRUE,
              test = t.test,
              y_position = c(0.9,1,1.1),#这里写显著性标记的位置和高度
              size=1,color="black",textsize = 4)+#这里写显著性的字体大小
  scale_y_continuous(limits = c(0.4,1.2),#这里设置一下y轴最低和最高的刻度
                     breaks = c(0.6,0.8,1,1.2))+#这里设置一下y轴刻度增量
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size = 1),
        axis.text.x = element_text(color = "black", size = 13),
        axis.text.y = element_text(color = "black",size = 13),
        legend.position = "none",
        axis.ticks = element_line(color="black",linewidth = 1))+
  labs(x=NULL,y=NULL)+
  scale_fill_manual(values = c("#FF0000","#800080","#FFA500","#0000FF"))
p2
