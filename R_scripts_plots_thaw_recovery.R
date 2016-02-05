#version 2: changed to pdfs and using source(themePM1)

library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/in_vitro_individual_prefreeze')



   
####### Recovery plots: prefreeze vs post-thaw  #######

dat<-read.csv("plerixafor_pre_post_correl.csv",header=T,check.names=T)
#check that the columns line up
#dat$Sample == dat$Sample.post
#just take some columns
#names(dat)
dat<-dat[,c(3,5,6,9:12,28:32)]
#add recovery %
dat$recovCFC <- 100*dat$PB.post.Total.CFC.x10e3.mL/dat$pre.Total.CFC.x10e3.mL
dat$recovCD34 <- 100*dat$PB.post.34.per.uL.of.original/dat$pre.CD34.x10e3.mL
dat$BMrecovCFC <- 100*dat$BM.post.Total.CFC.per.10e5/dat$pre.Total.CFC.per.10e5
#dat$BMrecovCD34 <- 100*dat$BM.post.34.percent.of.45/dat$pre.CD34.percent #this incorrectly looks at post-CD34%, not actual CD34 numbers which is more appropriate for recovery calculations
dat$BMrecovCD34 <- 100*dat$BM.post.34.per.original.TNC/dat$pre.CD34.percent
#log10 transform some columns
dat[4:11]<-log10(dat[4:11])
#rename group factors
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"
#rename timepoint factors
levels(dat$Timepoint)[levels(dat$Timepoint)=="BL"] <- "BL"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day -1"] <- "-24 h"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 0"] <- "+4 h"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 1"] <- "+24 h"


#correlation
#out<-cor.test(x=dat$post.Total.CFC.x10e3.mL,y=dat$pre.Total.CFC.x10e3.mL,use="complete.obs",method="spearman",exact=FALSE)
#dput(out,file="150220 correl prepost CFC.txt")



#barplot:
## PB ##
PB <- dat[dat$Sample=="Blood",]
PB <- PB[,c(1,3,13,14)]
PB <- melt(PB,id.vars=c("Group", "Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
PB <- summarySE(PB, measurevar="value", groupvars=c("Timepoint","assay"),na.rm=TRUE)
recovCD34<-PB[PB$assay=="recovCD34",]
recovCFC<-PB[PB$assay=="recovCFC",]


plot <- ggplot(recovCD34, aes(x=Timepoint, y=value))+
    geom_bar(stat="identity", fill="#000000") +
    geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.75) +
    ggtitle("PB CD34") +
    ylab("% recovery")+
    coord_cartesian(ylim=c(2.5,65)) + #, xlim=c(-0.6,1.2)) +
  #	annotation_logticks(sides = "lb") +
    scale_x_discrete("", labels=c("D0", "Post-G", "P+4 h", "P+24 h"))  +
    themePM1()
plot + 	theme(axis.text.x = element_text(angle = 45, hjust = 1, , vjust = 1))
ggsave(filename="PB CD34 prepost bar.pdf",width=4,height=5, units="cm")
#ggsave(filename="150720 PB CFC prepost bar.tiff",width=6,height=8.1,units="cm",dpi=600, compression="lzw")


plot <- ggplot(recovCFC, aes(x=Timepoint, y=value))+
	geom_bar(stat="identity", fill="#000000") +
	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.75) +
	ggtitle("PB CFC") +
	ylab("% recovery")+
	coord_cartesian(ylim=c(4,90)) + #, xlim=c(-0.6,1.2)) +
#	annotation_logticks(sides = "lb") +
    scale_x_discrete("", labels=c("D0", "Post-G", "P+4 h", "P+24 h"))  +
	themePM1()
plot + 	theme(axis.text.x = element_text(angle = 45, hjust = 1, , vjust = 1))
ggsave(filename="PB CFC prepost bar.pdf",width=4,height=5, units="cm")
#ggsave(filename="150720 PB CFC prepost bar.tiff",width=6,height=8.1,units="cm",dpi=600, compression="lzw")


## BM ##
BM <- dat[dat$Sample=="Marrow",]
BM <- BM[,c(1,3,15,16)]
BM <- melt(BM,id.vars=c("Group", "Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
BM <- summarySE(BM, measurevar="value", groupvars=c("Timepoint","assay"),na.rm=TRUE)
#omit Day-1:
BM <- BM[BM$Timepoint!="-24 h",]

recovCFC<-BM[BM$assay=="BMrecovCFC",]
recovCD34<-BM[BM$assay=="BMrecovCD34",]

plot <- ggplot(recovCD34, aes(x=Timepoint, y=value))+
  geom_bar(stat="identity", fill="#000000") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.75) +
  ggtitle("BM CD34") +
  ylab("% recovery")+
  coord_cartesian(ylim=c(2.5,65)) + #, xlim=c(-0.6,1.2)) +zz
  #	annotation_logticks(sides = "lb") +
  scale_x_discrete("", labels=c("D0", "P+4 h", "P+24 h"))  +
  themePM2()
plot + 	theme(axis.text.x = element_text(angle = 45, hjust = 1, , vjust = 1))
ggsave(filename="BM CD34 prepost bar.pdf",width=4,height=5, units="cm")
#ggsave(filename="150720 BM CFC prepost bar1.tiff",width=5,height=8.1,units="cm",dpi=600, compression="lzw")


plot <- ggplot(recovCFC, aes(x=Timepoint, y=value))+
	geom_bar(stat="identity", fill="#000000") +
	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.75) +
	ggtitle("BM CFC") +
	ylab("% recovery")+
	coord_cartesian(ylim=c(4,90)) + #, xlim=c(-0.6,1.2)) +
#	annotation_logticks(sides = "lb") +
    scale_x_discrete("", labels=c("D0", "P+4 h", "P+24 h"))  +
	themePM2()
plot + 	theme(axis.text.x = element_text(angle = 45, hjust = 1, , vjust = 1))
ggsave(filename="BM CFC prepost bar.pdf",width=4,height=5, units="cm")
#ggsave(filename="150720 BM CFC prepost bar1.tiff",width=5,height=8.1,units="cm",dpi=600, compression="lzw")














## PB ##
#scatterplot:
## CFC ##
dat<-melt(dat,id.vars=c("Group","Sample","Timepoint","post.Total.CFC.x10e3.mL"),variable.name="assay",value.name="value",na.rm=TRUE) #post.Total.CFC.x10e3.mL added to make plot between pre & post CFC
dat<-summarySE(dat, measurevar="value", groupvars=c("Group","Sample","Timepoint","post.Total.CFC.x10e3.mL","assay"),na.rm=TRUE)
PB<-dat[dat$Sample=="Blood",]
CFC<-PB[PB$assay=="pre.Total.CFC.x10e3.mL",]
tiff(filename="150220 PB CFC prepost xy.tiff",width=11.5,height=8.1,units="cm",res=600, compression="lzw")
ggplot(CFC, aes(x=post.Total.CFC.x10e3.mL, y=value, shape=Timepoint, color=Group))+
	scale_colour_manual(values=paletteB) +
	geom_point(size=3)+ 	# shape=1
	geom_errorbar(data=CFC,aes(ymin=value-se, ymax=value+se), lty=1, width=0.1, size=0.75) +
	scale_shape_manual(values=c(23,19,15,17))+ #specify shapes
	ggtitle("PB CFC") +
	ylab(expression(paste("Pre-freeze CFC x 10"^"3","/mL"))) +
	xlab(expression(paste("Post-thaw CFC x 10"^"3","/mL"))) +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
	scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
	coord_cartesian(ylim=c(-1,2), xlim=c(-0.6,1.2)) +
	annotation_logticks(sides = "lb") +
#    geom_smooth(method=lm,se=FALSE)+   # Add linear regression line, but not shaded confidence region
	theme_pres()
dev.off()

	
	
	
	
	
	
	
	
## BM ##
## CFC ##
tiff(filename="141123 BM CFC prepost xy.tiff",width=12.5,height=9,units="cm",res=600, compression="lzw")
ggplot(BM, aes(x=post.Total.CFC.per.10e5, y=pre.Total.CFC.per.10e5, shape=Timepoint, color=Group))+
	scale_colour_manual(values=paletteB) +
	geom_point(size=3)+ 	# shape=1
#	scale_shape_manual(values=c(23,19,15,17))+ #to specify which shapes to use
	ggtitle("BM CFC") +
	ylab(expression(paste("Pre-freeze CFC per 10"^"5"," cells"))) +
	xlab(expression(paste("Post-thaw CFC per 10"^"5"," cells"))) +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
	scale_x_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) +
#	coord_cartesian(ylim=c(-1,2), xlim=c(-0.6,1.2)) +
	annotation_logticks(sides = "lb") +
#    geom_smooth(method=lm,se=FALSE)+   # Add linear regression line, but not shaded confidence region
	theme_pres()
dev.off()

	
	
	
	



	
	
	
	
	
	
	
	
	


###dump
pre<-read.csv("plerixafor_prefreeze.csv",header=T,check.names=F)
post<-read.csv("plerixafor_pools_vitro.csv",header=T,check.names=F)

#take assay names
#assays<-names(dat[6:21])

#Convert to tall. http://www.cookbook-r.com/Manipulating_prea/Converting_prea_between_wide_and_long_format/
pre<-melt(pre,id.vars=c("SCA ID","Group","Patient","Sample","Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
post<-melt(post,id.vars=c("Pool","Sample","Timepoint","Group"),variable.name="assay",value.name="value",na.rm=TRUE) 
#log10 transform
pre$value<-log10(pre$value)
#rename group factors
levels(pre$Group)[levels(pre$Group)=="A"] <- "P"
levels(pre$Group)[levels(pre$Group)=="B"] <- "G+P"
#re-order factor levels:
#pre$Timepoint<-factor(pre$Timepoint,levels=c("BL","Day -1","Day 0","Day 1")) 

preBM<-pre[pre$Sample=="Marrow",]
prePB<-pre[pre$Sample=="Blood",]

#switch between BM & PB
preCFC

BLP<-dat[dat$Sample=="Marrow" & dat$Timepoint=="BL",] 


postBM<-post[post$Sample=="Marrow",]
postPB<-post[post$Sample=="Blood",]


##########Means with SEM#############

###BONE MARROW###
BMi<-BM[BM$Sample=="Marrow",]
BMs<-summarySE(BMi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval. #measurevar is the x-axis
#to only plot means, change plots to plot XXXs only.
#s means stats, i means individual values

#NCCi<-BMi[BMi$assay=="NCC x10e6/mL",]
NCCs<-BMs[BMs$assay=="NCC x10e6/mL",]
#CD34i<-BMi[BMi$assay=="CD34 percent",]
CD34s<-BMs[BMs$assay=="CD34 percent",]
#CFCti<-BMi[BMi$assay=="Total CFC per 10e5",]
CFCts<-BMs[BMs$assay=="Total CFC per 10e5",]









































################ Individual dots with group means - Timecourse plots ######################
may want to not pool baselines, or change shapes so that one group is obscured

###BONE MARROW###
BMi<-BM[BM$Sample=="Marrow",]
BMs<-summarySE(BMi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval. #measurevar is the x-axis
#to only plot means, change plots to plot XXXs only.
#s means stats, i means individual values

NCCi<-BMi[BMi$assay=="NCC x10e6/mL",]
NCCs<-BMs[BMs$assay=="NCC x10e6/mL",]
CD34i<-BMi[BMi$assay=="CD34 percent",]
CD34s<-BMs[BMs$assay=="CD34 percent",]
CFCti<-BMi[BMi$assay=="Total CFC per 10e5",]
CFCts<-BMs[BMs$assay=="Total CFC per 10e5",]


p1<-
#tiff(filename="BM_all_points_NCC.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(NCCi,aes(x=factor(Timepoint), y=as.numeric(value), group=Group))+
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(NCCs,aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=NCCs,aes(linetype=Group)) +
	ylab(expression(paste("TNC x 10"^"6","/mL"))) +
	ggtitle("BM") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
	coord_cartesian(ylim=c(0.8,2.1)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()
  
p2<-
#tiff(filename="BM_all_points_CD34.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(CD34i,aes(x=factor(Timepoint), y=as.numeric(value), group=Group)) + 
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=CD34s,aes(linetype=Group)) +
	ylab(expression(paste("CD34"^"+"," percent"))) +
	ggtitle("BM") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
	coord_cartesian(ylim=c(-1,0.4)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()

p3<-
#tiff(filename="BM_all points_tCFC.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(CFCti,aes(x=factor(Timepoint), y=as.numeric(value), group=Group)) + 
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=CFCts,aes(linetype=Group)) +
	ylab(expression(paste("Total CFC per 10"^"5"," cells"))) +
	ggtitle("BM") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
	coord_cartesian(ylim=c(1.8,2.7)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()



###PERIPHERAL BLOOD###
PBi<-dat[dat$Sample=="Blood",]
PBs<-summarySE(PBi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval. #measurevar is the x-axis
#to only plot means, change plots to plot XXXs only.
#s means stats, i means individual values

NCCi<-PBi[PBi$assay=="NCC x10e6/mL",]
NCCs<-PBs[PBs$assay=="NCC x10e6/mL",]
CD34i<-PBi[PBi$assay=="CD34 x10e3/mL",]
CD34s<-PBs[PBs$assay=="CD34 x10e3/mL",]
CFCti<-PBi[PBi$assay=="Total CFC x10e3/mL",]
CFCts<-PBs[PBs$assay=="Total CFC x10e3/mL",]


p4<-
#tiff(filename="PB_all_points_NCC.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(NCCi,aes(x=factor(Timepoint), y=as.numeric(value), group=Group))+
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(NCCs,aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=NCCs,aes(linetype=Group)) +
	ylab(expression(paste("TNC x 10"^"6","/mL"))) +
	ggtitle("PB") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
#	coord_cartesian(ylim=c(0.8,2.1)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()

p5<-
#tiff(filename="PB_all_points_CD34.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(CD34i,aes(x=factor(Timepoint), y=as.numeric(value), group=Group)) + 
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=CD34s,aes(linetype=Group)) +
	ylab(expression(paste("CD34 x 10"^"3","/mL"))) +
	ggtitle("PB") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
#	coord_cartesian(ylim=c(0.8,2.1)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()

p6<-
#tiff(filename="PB_all points_tCFC.tiff",width=12,height=12,units="cm",res=1000,compression="lzw")
ggplot(CFCti,aes(x=factor(Timepoint), y=as.numeric(value), group=Group)) + 
	geom_point(aes(shape=Group,colour=Patient),size=4) +
#	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	#this width is the width of the error bar handles
	geom_line(data=CFCts,aes(linetype=Group)) +
	ylab(expression(paste("Total CFC x 10"^"3","/mL"))) +
	ggtitle("PB") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + #sets y-range and override axis label (ENSURE THAT LIMITS AND BREAKS ARE EQUAL)
	coord_cartesian(ylim=c(-1,2)) +
	annotation_logticks(sides = "l") +
	scale_x_discrete("Timepoint")+ #, labels=c("??"))  +
	theme_pres()  
#graphics.off()


#letter size: 8.5 by 11 inches (216 mm x 279 mm)
pdf("140724b_all_points.pdf", width=35/2.54, height=21/2.54)
multiplot(p1, p4, p2, p5, p3, p6, cols=3)
dev.off()
   
  

  
  
  
  
  
  

#read in file
dat <- read.csv("plerixafor_prefreeze.csv",header=T,check.names=F) #,row.names=1)
#take assay names
assays<-names(dat[6:21])

#Convert to tall. http://www.cookbook-r.com/Manipulating_data/Converting_data_between_wide_and_long_format/
dat <- melt(dat,id.vars=c("ID","Group","Donor","Sample","Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
#log10 transform
dat$value <- log10(dat$value)
#rename group factors
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"
#re-order factor levels:
dat$Timepoint<-factor(dat$Timepoint,levels=c("BL","Day -1","Day 0","Day 1")) 
#re-name time points
levels(dat$Timepoint)[levels(dat$Timepoint)=="BL"] <- "BL"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day -1"] <- "-24 h"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 0"] <- "+4 h"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 1"] <- "+24 h"
#duplicate and pool baseline data
	#marrow
BLP<-dat[dat$Sample=="Marrow" & dat$Timepoint=="BL",]
BLP$Group<-"P"
BLG<-dat[dat$Sample=="Marrow" & dat$Timepoint=="BL",]
BLG$Group<-"G+P"
notBL<-dat[dat$Sample=="Marrow" & dat$Timepoint!="BL",]
BM<-rbind(BLP,BLG,notBL)
	#blood
BLP<-dat[dat$Sample=="Blood" & dat$Timepoint=="BL",]
BLP$Group<-"P"
BLG<-dat[dat$Sample=="Blood" & dat$Timepoint=="BL",]
BLG$Group<-"G+P"
notBL<-dat[dat$Sample=="Blood" & dat$Timepoint!="BL",]
PB<-rbind(BLP,BLG,notBL)



##########Kinetics of pre-freeze, Means with SEM#############

###BONE MARROW###
BMi<-BM[BM$Sample=="Marrow",]

#change scale of CFC from per 10e5
BMi[BMi$assay=="Total CFC per 10e5","value"]<-BMi[BMi$assay=="Total CFC per 10e5","value"]-2

BMs<-summarySE(BMi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval. #measurevar is the x-axis
#to only plot means, change plots to plot XXXs only.
#s means stats, i means individual values

#NCCi<-BMi[BMi$assay=="NCC x10e6/mL",]
NCCs<-BMs[BMs$assay=="NCC x10e6/mL",]
#CD34i<-BMi[BMi$assay=="CD34 percent",]
CD34s<-BMs[BMs$assay=="CD34 percent",]
#CFCti<-BMi[BMi$assay=="Total CFC per 10e5",]
CFCts<-BMs[BMs$assay=="Total CFC per 10e5",]


plot<-ggplot(data=NCCs,aes(x=factor(Timepoint), y=as.numeric(value), colour=Group, group=Group))+#, linetype=Group))+
	geom_point(aes(shape=Group),size=4) +
	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +	
	geom_line() +
	scale_colour_manual(values=paletteA) +
	ylab(expression(paste("TNC x 10"^"6","/mL"))) +
	ggtitle("") +
	scale_y_continuous(breaks=c(-2, -1, 0, 1, 2, 3),labels = c(0.01, 0.1, 1, 10, 100, 1000)) + 
	coord_cartesian(ylim=c(0,2)) +
	annotation_logticks(sides = "l") +
	guides(colour=F, shape=F)+
	scale_x_discrete(" ")+ #, labels=c("??"))  +
	themePM1()  
ggsave(filename="150720 BM individual NCC.tiff",width=8,height=8.1,units="cm",dpi=600, compression="lzw")
    
plot %+% CD34s +
	coord_cartesian(ylim=c(-1,1))+       
	ylab(expression(paste("CD34"^"+"," percent"))) 
ggsave(filename="150720 BM individual CD34.tiff",width=8,height=8.3,units="cm",dpi=600, compression="lzw")

plot %+% CFCts +
	coord_cartesian(ylim=c(-1,1))+ 
	ylab(expression(paste("Total CFC per 10"^"3"," cells"))) 
ggsave(filename="150720 BM individual tCFC.tiff",width=8,height=8.3,units="cm",dpi=600, compression="lzw")


	
	
	


###PERIPHERAL BLOOD###
PBi<-PB[PB$Sample=="Blood",]
PBs<-summarySE(PBi, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE) # summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval. #measurevar is the x-axis
#to only plot means, change plots to plot XXXs only.
#s means stats, i means individual values

#NCCi<-PBi[PBi$assay=="NCC x10e6/mL",]
NCCs<-PBs[PBs$assay=="NCC x10e6/mL",]
#CD34i<-PBi[PBi$assay=="CD34 x10e3/mL",]
CD34s<-PBs[PBs$assay=="CD34 x10e3/mL",]
#CFCti<-PBi[PBi$assay=="Total CFC x10e3/mL",]
CFCts<-PBs[PBs$assay=="Total CFC x10e3/mL",]


plot<-ggplot(data=NCCs,aes(x=factor(Timepoint), y=as.numeric(value), colour=Group, group=Group))+#, linetype=Group))+
	geom_point(aes(shape=Group),size=4) +
	scale_colour_manual(values=paletteA) +
	geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.2, size=0.75) +
	geom_line()+#data=NCCs,aes(linetype=Group)) +
	ylab(expression(paste("TNC x 10"^"6","/mL"))) +
	ggtitle(" ") +
	scale_y_continuous(breaks=c(-2, -1, -0.30103, 0, 1, 1.90309, 2, 3),labels = c(0.01, 0.1, 0.5, 1, 10, 80, 100, 1000)) + 
	coord_cartesian(ylim=c(-0.30103,1.90309))+
	annotation_logticks(sides = "l") +
	scale_x_discrete(" ")+ #, labels=c("??"))  +
	guides(colour=F, shape=F)+
	themePM1()  
graphics.off()
plot
ggsave(filename="150424 PB individual NCC.tiff",width=8,height=8.3,units="cm",dpi=600, compression="lzw")

plot %+% CD34s +
	ylab(expression(paste("CD34 x 10"^"3","/mL"))) +
	scale_y_continuous(breaks=c(-2, -1, -0.30103, 0, 1, 1.90309, 2, 3),labels = c(0.01, 0.1, 0.5, 1, 10, 80, 100, 1000)) + 
	coord_cartesian(ylim=c(-0.30103,1.90309))+         # (0.30103,2)) +
	ggtitle(" ")
ggsave(filename="150424 PB individual CD34.tiff",width=8,height=8.3,units="cm",dpi=600, compression="lzw")

plot %+% CFCts +
	ylab(expression(paste("CFC x 10"^"3","/mL"))) +
	scale_y_continuous(breaks=c(-2, -1, -0.30103, 0, 1, 1.90309, 2, 3),labels = c(0.01, 0.1, 0.5, 1, 10, 80, 100, 1000)) + 
	coord_cartesian(ylim=c(-0.30103, 1.90309)) +
	ggtitle(" ")	
ggsave(filename="150424 PB individual tCFC.tiff",width=8,height=8.3,units="cm",dpi=600, compression="lzw")



   





   
   