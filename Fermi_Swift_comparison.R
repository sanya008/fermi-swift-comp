
rm(list=ls()) # Clear workspace of RStudio

# libraries needed
library(FNN)
library(ggplot2)
library(cowplot)
library(GGally)
library(ggExtra)
library(ggforce)
library(celestial)
library(MASS)
library(xtable)
library(CCA)
library(CCP)
library(corrplot)
library(mclust)
library(bootstrap)


# parameters for plotting
fontsize=18/0.8


# Loading data  
drive <- "E:" 
directory <- "/Documents/2023/mnras/a_ref2/"
#
datafile <- "BAT23.dat"
BAT <- read.table(paste(drive,directory,datafile,sep=""),header=T) # Input BAT data
BAT$Flu <- BAT$Flu/1e7
BAT$Flue <- BAT$Flue/1e7

datafile <- "GBM23.dat"
GBM <- read.table(paste(drive,directory,datafile,sep=""),header=T,sep="|") # Input GBM data



############################################
# Preparing Swift BAT data                 #


# Omit GRBs before Fermi launch
nsw0 <- 335
nsw <- dim(BAT)[1] - nsw0
BAT <- BAT[1:nsw,]

# Converting trigger dates & times to POSIX
syear  <- substring(BAT$SName,1,2)
syear  <- paste("20",syear,sep="")
smonth <- substring(BAT$SName,3,4)
sday   <- substring(BAT$SName,5,6)
sdate  <- paste(syear,"-",smonth,"-",sday,sep="")

# Full date of trigger
swTrdate <- paste(sdate,BAT$Trtime,sep=" ")

# Swift satellite launch date
start <- "2004-11-20 17:16:01"

# Elapsed time of the trigger from the start of Swift in days
swTrtime <- rep(0.0,nsw)
for(i in 1:nsw){
  swTrtime[i] <- as.numeric(difftime(swTrdate[i],start,units="days"))
  }

# Converting RA, Dec into Descartes coordinates
sx <- cos(BAT$RAsw*pi/180)*cos(BAT$Desw*pi/180)
sy <- sin(BAT$RAsw*pi/180)*cos(BAT$Desw*pi/180)
sz <- sin(BAT$Desw*pi/180)

# Creating a position--time dataframe 
sxyzt <- data.frame(sx,sy,sz,swTrtime)

# Main characteristics of the data.frame
summary(sxyzt)
head(sxyzt)


############################################
# Preparing Fermi GBM data                 #


# Elapsed time of the trigger from the start of Swift in days
feTrtime <- rep(0,dim(GBM)[1])
for(i in 1:dim(GBM)[1]){
  feTrtime[i] <-as.numeric(difftime(GBM$Ftri_tim[i],start,units="days"))
  }

# Converting RA, Dec into Descartes coordinates
FRa  <- hms2deg(as.character(GBM$Fra),sep=" ")  # hms to degree
FDec <- dms2deg(as.character(GBM$Fdec),sep=" ") # dms to degree
fx <- cos(FRa*pi/180)*cos(FDec*pi/180)
fy <- sin(FRa*pi/180)*cos(FDec*pi/180)
fz <- sin(FDec*pi/180)

# Creating a position-time data.frame 
fxyzt <- data.frame(fx,fy,fz,feTrtime)

# Main charateristics of the data.frame
summary(fxyzt)
head(fxyzt)


# Searching for the next neighbors in the space-trigger_time frame
fsknn1 <- get.knnx(fxyzt,sxyzt,k=1)
nndist <- fsknn1$nn.dist
nnind <-  fsknn1$nn.index
nnres <- data.frame(nndist,nnind)
 
# Plotting histogram of nndist
status <- rep("single",nsw)
status[log10(nndist) < -2.5] <- "coinc"

pFScoinc1 <- ggplot(nnres, aes(x=log10(nndist))) +
    geom_histogram(bins=40,fill="#00bfc4",alpha=0.5,color="black") +
    ggtitle("Coincidence between Swift and Fermi GRBs") +
    theme_bw() + theme(plot.title = element_text(size=15)) +
    geom_vline(aes(xintercept=-2.5),color="#F8766D",linetype="dashed") +
    geom_text(x=-4.25,y=200,label="Swift-Fermi coincidences",color="#F8766D",size=fontsize*0.8*0.8/.pt) +
    labs(subtitle="", y="count",x="lg(NN distances)",title="") +
    theme(legend.text = element_text(size = fontsize*0.8),
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66),
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black')     )

#pfile1 <- paste(drive,directory,"fig1a-knn.pdf",sep="")
#ggsave(paste(drive,directory,"fig1a-knn.pdf",sep=""))

# Compute next neighbor angular difference between BAT and GBM
dpos <- sqrt((sx - fx[nnind])^2 + (sy - fy[nnind])^2 + (sz - fz[nnind])^2)
dang <- asin(dpos/2)*180/pi

# Compute next neighbor trigger time difference
dtim <- abs(swTrtime - feTrtime[nnind])

# Display coincidences in position - time plane
dfanti <- data.frame(dtim,dang,status)
pFScoinc2 <- ggplot(dfanti,aes(x=log10(dtim),y=log10(dang))) + 
  geom_point(aes(color=status)) + 
  theme_bw() + theme(legend.position="bottom") + 
  labs(subtitle="", y="lg(angular distance)",x="lg(trigger time)",title="") +
  theme(legend.text = element_text(size = fontsize*0.8), 
          axis.title = element_text(size = fontsize*0.8),
          text = element_text(size = fontsize,colour='black'),
          legend.title=element_blank(),
          legend.position = 'none',# legend.position=c(0.97,0.66), 
          legend.justification=c(1,0),
          strip.text = element_text(size = fontsize*0.8, colour='black'),
          axis.ticks = element_line(colour = 'black'),
          axis.text = element_text(size = fontsize*0.8, colour='black'))

plot_grid(pFScoinc1,pFScoinc2,nrow=2,ncol=1,rel_heights=c(3,5))

# Save last plot in pdf file
pfile1 <- paste(drive,directory,"FS_new_fig1.pdf",sep="")
ggsave(pfile1,width = 6.22, height = 6.22*1.3)




############################################
# Creating dataframes for couples & widows #


# Create Swift couple and widow
nbat <- dim(BAT)[1]
batst <- rep("widows",nbat)
batst[log10(nndist) < -2.5] <- "couples"

# Create Fermi couple and widow
ngbm <- dim(GBM)[1]
gbmst <- rep("widows",ngbm)
gbmst[nnind[log10(nndist) < -2.5]] <- "couples"

# Show results
table(batst)
table(gbmst)

# Create data frame for BAT-GBM next neighbors
swfenn <- data.frame(BAT,GBM[nnind,1:12])
 
# Create data frames for BAT-GBM couples
swfecoup <- swfenn[batst=="couples",]
names(swfecoup)[c(6,7,9,15,18,20)] <- c("ST90","Sflu","SPeak","FT90","Fflu","Fpeak")

# Display Swift -- Fermi  T90, Fluence, Peak flux scatter plots
lmpT <- lm(log10(swfecoup$FT90)~log10(swfecoup$ST90))
pT <- ggplot(swfecoup,aes(log10(ST90),log10(FT90))) + 
  stat_smooth(method = "lm", col = "red",alpha=0.5,level=erf(5/sqrt(2)),fullrange=T,linewidth=1) +
  geom_point(shape=19, color="black",size=1.5) + 
  theme_bw() + xlim(-2, 3) + ylim(-2, 3) + 
  ggtitle(expression(paste("Swift – Fermi ",T[90]," durations",sep=""))) + 
  geom_abline(intercept = 0, slope = 1, color="blue",  linetype="dashed", linewidth=1) + 
  annotate(geom = 'text', label = '- - -: y = x', x = -2, y = Inf, hjust = 0.0, vjust = 2, 
            color="blue", size = fontsize/.pt*0.8) + 
  annotate(geom = 'text', 
           label=paste("slope = ",round(lmpT$coefficients[2],2),'\u00B1',round(coef(summary(lmpT))[2,"Std. Error"],2),sep=""),hjust=1.0,x=3,y=-Inf,vjust=-2.4,color="red",size=fontsize/.pt*0.8) +
  annotate(geom = 'text', 
           label=paste("intercept = ",round(lmpT$coefficients[1],2),'\u00B1',round(coef(summary(lmpT))[1,"Std. Error"],2),sep=""),hjust=1.0,x=3,y=-Inf,vjust=-0.9,color="red",size=fontsize/.pt*0.8) +
  theme(plot.title = element_text(size=fontsize*0.8,hjust = 0.5)) + 
  xlab(expression(paste("Swift lg(",T[90],")",sep=""))) + 
  ylab(expression(paste("Fermi lg(",T[90],")",sep=""))) +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize*0.8*0.8,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black'))
 
#pT
 

   
 
lmpF <- lm(log10(swfecoup$Fflu)~log10(swfecoup$Sflu))
pF <- ggplot(swfecoup,aes(log10(Sflu),log10(Fflu))) + 
  stat_smooth(method = "lm", col="red",alpha=0.5,level=erf(5/sqrt(2)),fullrange=T,linewidth=1) +
  geom_point(shape=19, color="black",size=1.5) + 
  theme_bw() + ggtitle("Swift – Fermi fluences") + 
  geom_abline(intercept = 0, slope = 1, color="blue",  linetype="dashed", linewidth=1) + 
  annotate(geom = 'text', label = '- - -: y = x', x = -8, y = Inf, hjust = 0, vjust = 1.75, 
            color="blue", size=fontsize/.pt*0.8) + xlim(-8, -2.5) + ylim(-8, -2.5) + 
  annotate(geom = 'text', 
           label=paste("slope = ",round(lmpF$coefficients[2],2),'\u00B1',round(coef(summary(lmpF))[2,"Std. Error"],2),sep=""),hjust=1.0,x=-2.5,y=-Inf,vjust=-2.4,color="red",size=fontsize/.pt*0.8) +
  annotate(geom = 'text', 
           label=paste("intercept = ",round(lmpF$coefficients[1],2),'\u00B1',round(coef(summary(lmpF))[1,"Std. Error"],2),sep=""),hjust=1.0,x=-2.5,y=-Inf,vjust=-0.9,color="red",size=fontsize/.pt*0.8) +
  theme(plot.title = element_text(size=fontsize*0.8,hjust = 0.5)) + 
  xlab(expression(paste("Swift lg(","fluence",")",sep=""))) + 
  ylab(expression(paste("Fermi lg(","fluence",")",sep=""))) +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize*0.8*0.8,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black'))

#pF

 
lmpP <- lm(log10(swfecoup$Fpeak)~log10(swfecoup$SPeak))
pP <- ggplot(swfecoup,aes(log10(SPeak),log10(Fpeak))) + 
  stat_smooth(method = "lm", col="red",alpha=0.5,level=erf(5/sqrt(2)),fullrange=T,linewidth=1) +
  geom_point(shape=19, color="black",size=1.5) + 
  theme_bw() + ggtitle("Swift – Fermi 1024 ms peak flux") + 
  geom_abline(intercept = 0, slope = 1, color="blue",  linetype="dashed", linewidth=1) + 
  annotate(geom = 'text', label = '- - -: y = x', x = -1, y = Inf, hjust = 0, vjust = 1.5, 
            color="blue", size=fontsize/.pt*0.8)  + xlim(-1, 2.5) + ylim(-1, 2.5) + 
  annotate(geom = 'text', 
           label=paste("slope = ",round(lmpP$coefficients[2],2),'\u00B1',round(coef(summary(lmpP))[2,"Std. Error"],2),sep=""),hjust=1.0,x=2.5,y=-Inf,vjust=-2.4,color="red",size=fontsize/.pt*0.8) +
  annotate(geom = 'text', 
           label=paste("intercept = ",round(lmpP$coefficients[1],2),'\u00B1',round(coef(summary(lmpP))[1,"Std. Error"],2),sep=""),hjust=1.0,x=2.5,y=-Inf,vjust=-0.9,color="red",size=fontsize/.pt*0.8) +
  theme(plot.title = element_text(size=fontsize*0.8,hjust = 0.5)) + 
  xlab(expression(paste("Swift lg(","peak flux",")",sep=""))) + 
  ylab(expression(paste("Fermi lg(","peak flux",")",sep=""))) +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize*0.8*0.8,colour='black'),
        legend.title=element_blank(),
        legend.position = 'none',# legend.position=c(0.97,0.66), 
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black'))

#pP

plot_grid(pT,pF,pP,ncol=1,nrow=3)

pfile1 <- paste(drive,directory,"FS_new_fig2.pdf",sep="")
ggsave(pfile1,width = 6.22, height = 6.22*1.3)



############################################
# Performing LDA on BAT                    #


# Prepare input for LDA
 lBAT <- data.frame(log10(BAT[,6]),log10(BAT[,7]),log10(BAT[,9]),batst)
 nlBAT <- na.omit(lBAT)
 names(nlBAT) <- c("nlST90","nlSFlu","nlSP1","batst")
#
  Swlda <- lda(nlBAT[,1:3],nlBAT[,4])
 Swlda
#
# Computing best discriminating variable
#
 LD1 <- Swlda$scaling
 SLD1 <- as.matrix(nlBAT[1:3]) %*% LD1

# Testing the discriminant power of LD1, using ANOVA

 dfbat <-data.frame(nlBAT,SLD1)
 dfcor <- dfbat[,c(4,5)]

 status <- nlBAT[,4]
 aov.res <- aov(SLD1 ~ status,data=dfcor)
 summary(aov.res)

 scale(Swlda$scaling, F, sqrt(colSums(Swlda$scaling^2)))
 
 
# Displaying discriminant power of LD1
pLD1 <- ggplot(dfbat, aes(0-LD1, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("Greatest BAT couples–widows difference") + theme_bw() + 
   ylab("Density") + xlab("LD1") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

# Display discriminant power on T90, Fluence, Peak flux
pT90 <- ggplot(dfbat, aes(nlST90, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle(expression(paste("BAT ",T[90]," couples–widows difference"))) + theme_bw() + 
   ylab("Density") + xlab(expression(paste("lg(",T[90],")",sep=""))) +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

pFlu <- ggplot(dfbat, aes(nlSFlu, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("BAT fluence couples–widows difference") +theme_bw() + 
   ylab("Density") + xlab("lg(fluence)") +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize*0.8*0.8,colour='black'),
        legend.title=element_blank(),
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black'),
        plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

pPeak <- ggplot(dfbat, aes(nlSP1, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("BAT peak flux couples–widows difference") +theme_bw() + 
   ylab("Density") + xlab("lg(peak flux)") +
  theme(legend.text = element_text(size = fontsize*0.8), 
        axis.title = element_text(size = fontsize*0.8),
        text = element_text(size = fontsize*0.8*0.8,colour='black'),
        legend.title=element_blank(),
        legend.justification=c(1,0),
        strip.text = element_text(size = fontsize*0.8, colour='black'),
        axis.ticks = element_line(colour = 'black'),
        axis.text = element_text(size = fontsize*0.8, colour='black'),
        plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

plot_grid(pLD1,pT90,pFlu,pPeak,nrow=4,ncol=1)

pfile3 <- paste(drive,directory,"FS_new_fig3.pdf",sep="")
ggsave(pfile3,width = 6.22, height = 6.22*1.6)



# Create tables for imputing LateX
aggbat <- aggregate(dfbat[,c(5,1:3)],list(dfbat[,4]),mean)
names(aggbat) <- c("Group","LD1","log10(T90)","log10(Flu)","log10(Peak)")

xtable(aggbat,caption="Differences  between coup and widow groups in BAT LDA.",label="batlda")


############################################
# Performing LDA on GBM                    #

# Prepare input for LDA
 lGBM <- data.frame(log10(GBM[,5]),log10(GBM[,8]),log10(GBM[,10]),gbmst)
 nlGBM <- na.omit(lGBM)
 names(nlGBM) <- c("nlFT90","nlFFlu","nlFP1","gbmst")

 Felda <- lda(nlGBM[,1:3],nlGBM[,4])
 Felda

# Computing best discriminating variable

 LD1 <- Felda$scaling
 FLD1 <- as.matrix(nlGBM[1:3]) %*% LD1

 dfgbm <-data.frame(nlGBM,FLD1)
 dfcor <- dfgbm[,c(4,5)]

 status <- nlGBM[,4]
 aov.res <- aov(FLD1 ~ status,data=dfcor)
 summary(aov.res)

 scale(Felda$scaling, F, sqrt(colSums(Felda$scaling^2)))
 
# Displaying discriminant power of LD1

 pFLD1 <- ggplot(dfgbm, aes(0-FLD1, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("Greatest GBM couples–widows difference") + theme_bw() + 
   ylab("Density") + xlab("LD1") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))


# Display discriminant power on T90, Fluence, Peakflux

 pFT90 <- ggplot(dfgbm, aes(nlFT90, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle(expression(paste("GBM ",T[90]," couples–widows difference"))) + theme_bw() + 
   ylab("Density") + xlab(expression(paste("lg(",T[90],")",sep=""))) +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

 pFFlu <- ggplot(dfgbm, aes(nlFFlu, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("GBM fluence couples–widows difference") + theme_bw() + 
   ylab("Density") + xlab("lg(fluence)") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

 pFPeak <- ggplot(dfgbm, aes(nlFP1, fill = status)) + geom_density(alpha = 0.5) + 
   ggtitle("GBM peak flux couples–widows difference") + theme_bw() + 
   theme(plot.title = element_text(size=12,hjust = 0.5)) + theme(legend.title=element_blank()) + 
   ylab("Density") + xlab("lg(peak flux)") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

 plot_grid(pFLD1,pFT90,pFFlu,pFPeak,nrow=4,ncol=1)

 pfile4 <- paste(drive,directory,"FS_new_fig4.pdf",sep="")
 ggsave(pfile4,width = 6.22, height = 6.22*1.6)

# Create table for LateX
 agggbm <- aggregate(dfgbm[,c(5,1:3)],list(dfgbm[,4]),mean)
 names(agggbm) <- c("Group","LD1","log10(T90)","log10(Flu)","log10(Peak)")

 xtable(agggbm,caption="Differences  between coup and widow groups in GBM LDA.",label="gbmlda")


#####################################################
# Canonical correlation between BAT and GBM couples #


# Prepare input for canonical correlation

 nSwFecoup <- na.omit(log10(swfecoup[,c(6,7,9,15,18,20)]))
 inSw <- nSwFecoup[,c(1:3)]
 inFe <- nSwFecoup[,c(4:6)]

# Perform canonical correlations

cc.res <- cc(inSw,inFe)

# Retrieve canonical scores for BAT (xscores) and GBM (yscores)

 U <- cc.res$scores$xscores
 V <- cc.res$scores$yscores

# Estimate significance with Wilks' Lambda

 rho <- cc.res$cor
 N <- dim(inSw)[1]
 p <- dim(inSw)[2]
 q <- dim(inFe)[2] 

 p.asym(rho, N, p, q, tstat = "Wilks")

 
 
# Display results

 
# Plot U-BAT, pairwise

 UBAT <- data.frame(U,inSw)
 names(UBAT) <- c("U1","U2","U3","T90","fluence","peak flux")

 pUBAT <-  ggpairs(UBAT, columns = c(1:6),
  lower = list(continuous=wrap(ggally_density,size =0.1,color="blue")),
  upper = list(continuous=wrap(ggally_cor,size=fontsize*0.8*0.8*0.8/.pt,color="black")),
  diag  = list(continuous=wrap(ggally_barDiag,bins=20,fill="lightblue",color="blue"))) + 
   theme_bw() + removeGrid() + ggtitle("BAT & U data") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8*0.5, colour='black'),
         plot.title = element_text(size=fontsize*0.8*0.8,hjust = 0.5))
 
 
 pUBAT
 pfile5 <- paste(drive,directory,"FS_new_fig5.pdf",sep="")
 ggsave(pfile5,width=6.22,height=6.22)

 
# Plot V-BAT, pairwise 

 VBAT <- data.frame(V,inSw)
 names(VBAT) <- c("V1","V2","V3","T90","fluence","peak flux")

 pVBAT <- ggpairs(VBAT, columns = c(1:6),
  lower = list(continuous=wrap(ggally_density,size=0.1,color="blue")),
  upper = list(continuous=wrap(ggally_cor,size=fontsize*0.8*0.8*0.8/.pt,color="black")),
  diag  = list(continuous=wrap(ggally_barDiag,bins=20,fill="lightblue",color="blue"))) + 
   theme_bw() + removeGrid() + ggtitle("BAT & V data") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8*0.5, colour='black'),
         plot.title = element_text(size=fontsize*0.8*0.8,hjust = 0.5))

 pVBAT
 pfile6 <- paste(drive,directory,"FS_new_fig6.pdf",sep="")
 ggsave(pfile6,width=6.22,height=6.22)

 
# Plot U-GBM, pairwise 
 
 UGBM <- data.frame(U,inFe)
 names(UGBM) <- c("U1","U2","U3","T90","fluence","peak flux")

 pUGBM <- ggpairs(UGBM, columns = c(1:6),
  lower = list(continuous = wrap(ggally_density,size=0.1,color="darkred")),
  upper = list(continuous=wrap(ggally_cor,size=fontsize*0.8*0.8*0.8/.pt,color="black")),
  diag  = list(continuous = wrap(ggally_barDiag,bins=20,fill="pink",color="darkred")))  + 
   theme_bw() + removeGrid() + ggtitle("GBM & U data") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8*0.5, colour='black'),
         plot.title = element_text(size=fontsize*0.8*0.8,hjust = 0.5))

 pUGBM
 pfile7 <- paste(drive,directory,"FS_new_fig7.pdf",sep="")
 ggsave(pfile7,width=6.22,height=6.22)

# Plot V-GBM, pairwise 
 
  VGBM <- data.frame(V,inFe)
  names(VGBM) <- c("V1","V2","V3","T90","fluence","peak flux")
#
 pVGBM <- ggpairs(VGBM, columns = c(1:6),
  lower = list(continuous=wrap(ggally_density, size =0.1, color = "darkred")),
  upper = list(continuous=wrap(ggally_cor,size=fontsize*0.8*0.8*0.8/.pt,color="black")),
  diag  = list(continuous=wrap(ggally_barDiag,bins=20,fill="pink",color="darkred")))  + 
   theme_bw() + removeGrid() + ggtitle("GBM & V data") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8*0.5, colour='black'),
         plot.title = element_text(size=fontsize*0.8*0.8,hjust = 0.5))

 pVGBM
 pfile8 <- paste(drive,directory,"FS_new_fig8.pdf",sep="")
 ggsave(pfile8,width=6.22,height=6.22)


######################################
# Color display of correlations      #


# Compute structure matrices for u and v values

  UV  <- data.frame(U,V)
 cxuv <- cor(inSw,UV,use="pairwise.complete")
 cyuv <- cor(inFe,UV,use="pairwise.complete")

# Cross correlations of U and V with BAT and GBM data

 colnames(cxuv) <- c("U1","U2","U3","V1","V2","V3")
 rownames(cxuv) <- c("$Swift~T[90]","Swift fluence","Swift peak flux")
 colnames(cyuv) <- c("U1","U2","U3","V1","V2","V3")
 rownames(cyuv) <- c("$Fermi~T[90]","Fermi fluence","Fermi peak flux")
 cxyuv <- round(rbind(cxuv,cyuv),3)
 cxyuv

# Estimate significance

 psig <- 0.003 # 3sigma

 colx <- dim(inSw)[2]
 coly <- dim(inFe)[2]
 coluv <- dim(UV)[2]
 pmatx <- matrix(rep(0,colx*coluv),colx,coluv)
 pmaty <- matrix(rep(0,coly*coluv),coly,coluv)

 for(i in 1:colx){
  for(j in 1:coluv){pmatx[i,j] <- cor.test(inSw[,i],UV[,j])$p.value
   cxuv[i,j] <- ifelse(pmatx[i,j]>psig,0,cxuv[i,j])}
  }
 for(i in 1:coly){
  for(j in 1:coluv){pmaty[i,j] <- cor.test(inFe[,i],UV[,j])$p.value
   cyuv[i,j] <- ifelse(pmaty[i,j]>psig,0,cyuv[i,j])}
  }
 
 # Structure matrix resulted
 
  cxyuv <- round(rbind(cxuv,cyuv),3)
  cxyuv

# Colored display of the structure matrix
pfile9 <- paste(drive,directory,"FS_new_fig9.pdf",sep="")
pdf(file=pfile9)
  corrplot(cxyuv[,c(1,4,2,5,3,6)],is.corr=T,tl.cex=1.2,tl.col='black',
                     cl.cex=1.2,cl.ratio=0.2,cl.align.text='r',cl.offset=0.1,cl.length=5)
dev.off()

#ggsave(filename=pfile9,plot=replayPlot(cplot))



#############################################
# Compare Swift and Fermi T90 distributions #

 T90 <- c(inSw[,1],inFe[,1])
 satellite <- c(rep("Swift",length(inSw[,1])),rep("Fermi",length(inFe[,1])))
 dfT90 <- data.frame(T90,satellite)
#
# Display jointly detected Swift and Fermi T90 distributions
#
pT90 <- ggplot(dfT90,aes(x=T90,color=satellite,fill=satellite)) +   
   geom_density(linewidth=1,alpha=0.4) + theme_bw() + 
   ggtitle(expression(paste("Swift–Fermi ",T[90]," distributions"))) + 
   ylab("Density") + xlab(expression(paste("lg(",T[90],")",sep=""))) +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.position=c(0.25,0.75), 
         legend.justification=c(0.5,0.5),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5)) + 
     scale_fill_manual(values=c("darkred","blue")) +
     scale_color_manual(values=c("darkred","blue"))



 pT90
 
 pfile10 <- paste(drive,directory,"FS_new_fig10.pdf",sep="")
 ggsave(pfile10,width=6.22)

 
 
 
############################################### 
# Perform Mclust on Fermi and Swift T90 data
 
 Stest <- as.data.frame(NA)
 for(i in 1:7){
   Sclust <- Mclust(inSw[,1],G=c(i),modelNames = "V")
   Stest[i,1] <- Sclust$loglik
   Stest[i,2] <- Sclust$BIC[1]
 }
 colnames(Stest) <- c("loglik","BIC")
 
 Ftest <- as.data.frame(NA)
 for(i in 1:7){
   Fclust <- Mclust(inFe[,1],G=c(i),modelNames = "V")
   Ftest[i,1] <- Fclust$loglik
   Ftest[i,2] <- Fclust$BIC[1]
 }
 colnames(Ftest) <- c("loglik","BIC")
 
 
# Data frame for Swift and Fermi clusters
 Sbic <- Stest$BIC#[,1]
 Fbic <- Ftest$BIC#[,1]
 
 Sdf <- data.frame(nclust=c(1:7),BIC=Sbic)
 Fdf <- data.frame(nclust=c(1:7),BIC=Fbic)
 
 
 # Error estimates with jackknife - Swift
 sds <- as.data.frame(NA)
 for (k in 1:7) {
   theta <- function(x){mclustBIC(x,modelNames = "V", G=k,verbose=F)}
   jack_data <- jackknife(inSw[,1],theta)
   sds[k,1] <- sd(jack_data$jack.values)
   }
 Sdf <- data.frame(Sdf,sds)
 colnames(Sdf)[3] <- "dev_jack"
 
  
 # Error estimates with jackknife - Fermi
 sds <- data.frame(NA)
 for (k in 1:7) {
   theta <- function(x){mclustBIC(x,modelNames = "V", G=k,verbose=F)}
   jack_data <- jackknife(inFe[,1],theta)
   sds[k,1] <- sd(jack_data$jack.values)
 }
 Fdf <- data.frame(Fdf,sds)
 colnames(Fdf)[3] <- "dev_jack"
 
 
# Display the frames with jackknife error estimates
 
 Fncp <- ggplot(Fdf,aes(x=nclust,y=BIC)) + 
   geom_line(color="darkred") + geom_point(size=3,color="darkred") + 
   geom_pointrange(aes(ymin=BIC-dev_jack, ymax=BIC+dev_jack),color="darkred") + 
   xlim(0,8) +  
   scale_x_continuous(breaks=seq(0,8,1)) + 
   scale_y_continuous(limits=c(-895,-755),breaks=c(-875,-825,-775))+
   theme_bw() + ggtitle(expression(paste("Fermi ",T[90]," clusters",sep=""))) + 
   xlab("# of clusters") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

 Sncp <- ggplot(Sdf,aes(x=nclust,y=BIC)) + 
   geom_line(color="blue") + geom_point(size=3,color="blue") +
   geom_pointrange(aes(ymin=BIC-dev_jack, ymax=BIC+dev_jack),color="blue") + 
     xlim(0,8) +  
   scale_x_continuous(breaks=seq(0,8,1)) + 
   scale_y_continuous(limits=c(-990,-840),breaks=c(-875,-925,-975)) +
   theme_bw() + ggtitle(expression(paste("Swift ",T[90]," clusters",sep=""))) + 
   xlab("# of clusters") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         #legend.position = 'none',# legend.position=c(0.97,0.66), 
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))

 plot_grid(Sncp,Fncp,nrow=2,ncol=1)

 pfile11 <- paste(drive,directory,"FS_new_fig11d.pdf",sep="")
 ggsave(pfile11,width=6.22,height=6.22)

 # Display the frames without errors 
 
 Fncp <- ggplot(Fdf,aes(x=nclust,y=BIC)) + 
   geom_line(color="darkred") + geom_point(size=3,color="darkred") + 
   xlim(0,8) +  
   scale_x_continuous(breaks=seq(0,8,1)) + scale_y_continuous(limits=c(-905,-755),breaks=c(-875,-825,-775))+
   theme_bw() + ggtitle(expression(paste("Fermi ",T[90]," clusters",sep=""))) + 
   xlab("# of clusters") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))
 
 Sncp <- ggplot(Sdf,aes(x=nclust,y=BIC)) + 
   geom_line(color="blue") + geom_point(size=3,color="blue") +
   xlim(0,8) +  
   scale_x_continuous(breaks=seq(0,8,1)) + 
   scale_y_continuous(limits=c(-1000,-820),breaks=c(-850,-900,-950,-1000)) +
   theme_bw() + ggtitle(expression(paste("Swift ",T[90]," clusters",sep=""))) + 
   xlab("# of clusters") +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         #legend.position = 'none',# legend.position=c(0.97,0.66), 
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))
 
 plot_grid(Sncp,Fncp,nrow=2,ncol=1)
 
 pfile11 <- paste(drive,directory,"FS_new_fig11.pdf",sep="")
 ggsave(pfile11,width=6.22,height=6.22)
 
 
 
 
 
 ######
 # Bootstrap Likelihood Ratio Test statistic
 
 
 # Swift
 
 X <- inSw[,1]
 BICx <- mclustBIC(X)
 plot(BICx)

 summary(BICx)
 
 mod1 <- Mclust(X, x = BICx, modelName = "V")
 summary(mod1, parameters = TRUE)
 
 plot(mod1, what = "classification")
 plot(mod1, what = "uncertainty")
 
 LRT_Swift <- mclustBootstrapLRT(X,modelName="V",nboot=100000,
                                 verbose=T,maxG=3,level=erf(5/sqrt(2)))
 LRT_Swift
 
 
 # Fermi
 
 Y <- inFe[,1]
 BICy <- mclustBIC(Y)
 plot(BICy)
 
 summary(BICy)
 
 mod2 <- Mclust(Y, x = BICy, modelName = "V")
 summary(mod2, parameters = TRUE)
 
 plot(mod2, what = "classification")
 plot(mod2, what = "uncertainty")
 

 LRT_Fermi <- mclustBootstrapLRT(Y,modelName="V",nboot=100000,
                                 verbose=T,maxG=3,level=erf(5/sqrt(2)))
 LRT_Fermi
 
 
 ######
 # Other plots
 ###
 
 ######
 # T90 distributions
 ###
 
ggplot(BAT, aes(log10(T90))) + geom_density(alpha = 0.5,linewidth=1,fill="blue") + 
   ggtitle(expression(paste("BAT ",T[90]))) + theme_bw() + #xlim(1.5,2) +
   ylab("Density") + xlab(expression(paste("lg(",T[90],")",sep=""))) +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         #legend.position = 'none',# legend.position=c(0.97,0.66), 
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))
 
 ggplot(GBM, aes(log10(FT90))) + geom_density(alpha = 0.5,linewidth=1,fill="darkred") + 
   ggtitle(expression(paste("GBM ",T[90]))) + theme_bw() + 
   ylab("Density") + xlab(expression(paste("lg(",T[90],")",sep=""))) + #xlim(1.5,1.7) +
   theme(legend.text = element_text(size = fontsize*0.8), 
         axis.title = element_text(size = fontsize*0.8),
         text = element_text(size = fontsize*0.8*0.8,colour='black'),
         legend.title=element_blank(),
         #legend.position = 'none',# legend.position=c(0.97,0.66), 
         legend.justification=c(1,0),
         strip.text = element_text(size = fontsize*0.8, colour='black'),
         axis.ticks = element_line(colour = 'black'),
         axis.text = element_text(size = fontsize*0.8, colour='black'),
         plot.title = element_text(size=fontsize*0.8,hjust = 0.5))
 
 

