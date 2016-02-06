# R scripts for CD34 and CFC recovery, for plerixafor +/-G-CSF paper, 2016.
# Shows barplots representing % of post-thaw results compared to pre-freeze. 
# Paul Miller, paulhmiller@gmail.com

library(reshape2)
library(ggplot2)
library(grid)
library(plyr)

source("C:/Users/paulm/Documents/R/source/functions.R")
source("C:/Users/paulm/Documents/R/source/plotting_themes.R")
setwd('C:/Users/paulm/CRC Paul/PROJECTS/Plerixafor/figures/in_vitro_individual_prefreeze')

dat <- read.csv("plerixafor_pre_post_correl.csv", header=T, check.names=T)
# Check that the columns line up
dat$Sample == dat$Sample.post
# Just take some columns
dat <- dat[,c(3,5,6,9:12,28:32)]

# Add % recovery
dat$recovCFC <- 100*dat$PB.post.Total.CFC.x10e3.mL/dat$pre.Total.CFC.x10e3.mL
dat$recovCD34 <- 100*dat$PB.post.34.per.uL.of.original/dat$pre.CD34.x10e3.mL
dat$BMrecovCFC <- 100*dat$BM.post.Total.CFC.per.10e5/dat$pre.Total.CFC.per.10e5
#dat$BMrecovCD34 <- 100*dat$BM.post.34.percent.of.45/dat$pre.CD34.percent #this incorrectly looks at post-CD34%, not actual CD34 numbers which is more appropriate for recovery calculations
dat$BMrecovCD34 <- 100*dat$BM.post.34.per.original.TNC/dat$pre.CD34.percent

# Log10 transform some columns
dat[4:11] <- log10(dat[4:11])

# Rename group factors
levels(dat$Group)[levels(dat$Group)=="A"] <- "P"
levels(dat$Group)[levels(dat$Group)=="B"] <- "G+P"

# Pool BL values
levels(dat$Group) <- c(levels(dat$Group), "BL")
dat[dat$Timepoint=="BL", "Group"] <- "BL"

# Rename timepoint and Group factors
levels(dat$Timepoint)[levels(dat$Timepoint)=="BL"] <- "D0"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day -1"] <- "Post-G"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 0"] <- "P+4 h"
levels(dat$Timepoint)[levels(dat$Timepoint)=="Day 1"] <- "P+24 h"
levels(dat$Group)[levels(dat$Group)=="BL"] <- "D0"

# Extract PB data
PB <- dat[dat$Sample=="Blood",]
PB <- PB[,c(1,3,13,14)]
PB <- melt(PB,id.vars=c("Group", "Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
PB <- summarySE(PB, measurevar="value", groupvars=c("Group","Timepoint","assay"),na.rm=TRUE)
PBrecovCD34<-PB[PB$assay=="recovCD34",]
PBrecovCFC<-PB[PB$assay=="recovCFC",]

# Extract BM data
BM <- dat[dat$Sample=="Marrow",]
BM <- BM[,c(1,3,15,16)]
BM <- melt(BM,id.vars=c("Group", "Timepoint"),variable.name="assay",value.name="value",na.rm=TRUE) 
BM <- summarySE(BM, measurevar="value", groupvars=c("Group", "Timepoint","assay"),na.rm=TRUE)
BM <- BM[BM$Timepoint!="Post-G",]  # To omit Day-1, as we don't trust this data
BMrecovCFC <- BM[BM$assay=="BMrecovCFC",]
BMrecovCD34 <- BM[BM$assay=="BMrecovCD34",]

# Global plotting specs
PBcols <- c("#CD0000","#A020F0","white") 
PBerrcols <- c("black", "#4169E1", "#CD0000", "#A020F0", "#CD0000", "#A020F0")
BMcols <- c("#CD0000","#A020F0","white") 
BMerrcols <- c("black", "#CD0000", "#A020F0", "#CD0000")

# Plots
p1 <- ggplot(PBrecovCD34, aes(x=Group, fill=Group, y=value))+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.2, position=position_dodge(0.9), colour=PBerrcols) +
    geom_bar(stat="identity",position=position_dodge(), colour="black", size=0.2)+
	ggtitle("PB CD34") +
    ylab("% recovery")+
    coord_cartesian(ylim=c(3.2,70)) + 
    scale_x_discrete(breaks=NULL) +
	guides(fill=F, color=F) +
	scale_fill_manual(values=PBcols) + 
	facet_grid(.~Timepoint, scale="free_x", space="free_x", switch = "x") +
    themePM2()

p2 <- p1 %+% PBrecovCFC +
	ggtitle("PB CFC") +
	ylab("% recovery")+
	coord_cartesian(ylim=c(4,91))

p3 <- ggplot(BMrecovCD34, aes(x=Group, fill=Group, y=value))+
    geom_errorbar(aes(ymin=value-se, ymax=value+se), lty=1, width=0.5, size=0.2, position=position_dodge(0.9), colour=BMerrcols) +
    geom_bar(stat="identity", position=position_dodge(),colour="black", size=0.2)+
    ggtitle("BM CD34") +
    ylab("% recovery")+
    coord_cartesian(ylim=c(3.2,70)) +
    scale_x_discrete(breaks=NULL) +
	guides(fill=F, color=F) +
	scale_fill_manual(values=BMcols) + 
	facet_grid(.~Timepoint, scale="free_x", space="free_x", switch = "x") +
    themePM2()

p4 <- p3 %+% BMrecovCFC +
	ggtitle("BM CFC") +
	ylab("% recovery")+
	coord_cartesian(ylim=c(4,91)) )  

# Make single PDF
pdf("recovery_plots.pdf", width=16/2.54, height=11.5/2.54)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2,widths=c(0.58,0.42))))
print(p1, vp=viewport(layout.pos.row=1,layout.pos.col=1))
print(p2, vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(p3, vp=viewport(layout.pos.row=1,layout.pos.col=2))
print(p4, vp=viewport(layout.pos.row=2,layout.pos.col=2))
dev.off()
