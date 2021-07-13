library(ggplot2)
library(gplots)
library(scatterplot3d)
library(reshape)
library(gridExtra)
library(car)
library(vegan)
library(reshape)
library(reshape2)
library(broom)
library(PMCMR)
rm(list=ls())


#######   Read in data  - - - - - - - - - - - - - - - -
periods<-read.csv("Porites-GCB-21.csv")[,-1]
summary(periods)
dim(periods)

#######   stats on day 14 physiological metrics ----------------------------

periods<-subset(periods, X1=="C")
periods<-periods[,1:21]

har<-t(periods[,9:21])
colnames(har)<-periods$local
stat<-har
rs<-rownames(stat)
cls<-colnames(stat)
Boo<-dim(stat)[[1]]
datum <- list()
atum <- list()
Tukey<-list()

for (i in 1:Boo)
{		H<-rownames(stat)[i]
		if (as.matrix(shapiro.test(resid(aov(stat[1,]~periods$X4)))[[2]])>0.05)
		{	atum[[H]]<-tidy(aov(stat[i,]~periods$X4))[1,6]
			if(tidy(aov(stat[i,]~periods$X4))[1,6]<0.05)
			{
				Tukey[[H]]<-TukeyHSD(aov(stat[i,]~periods$X4))
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))	
			}
		}
		else if (as.matrix(shapiro.test(resid(aov(log(stat[i,])~periods$X4)))[[2]])>0.05)
		{	atum[[H]]<-tidy(aov(log(stat[i,])~periods$X4))[1,6]
			if(tidy(aov(log(stat[i,])~periods$X4))[1,6]<0.05)
			{
				Tukey[[H]]<-TukeyHSD(aov(log(stat[i,])~periods$X4))
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))	
			}
		}
		else if (as.matrix(shapiro.test(resid(aov(sqrt(stat[i,])~periods$X4)))[[2]])>0.05)
		{	atum[[H]]<-tidy(aov(sqrt(stat[i,])~periods$X4))[1,6]
			if(tidy(aov(sqrt(stat[i,])~periods$X4))[1,6]<0.05)
			{
				Tukey[[H]]<-TukeyHSD(aov(sqrt(stat[i,])~periods$X4))
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
		else
		{	atum[[H]]<-kruskal.test(stat[i,]~factor(periods$X4))[[3]]
			if(kruskal.test(stat[i,]~factor(periods$X4))[[3]]<0.05)
			{
				Tukey[[H]]<-posthoc.kruskal.nemenyi.test(x=stat[i,], g=as.factor(periods$X4), method="Tukey", p.adjust.method="bonferonni")
				datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
			}
		}
}
allerie2<-data.matrix(do.call(rbind, datum))
stater2<-data.matrix(do.call(rbind, atum))
#stater1<-data.matrix(do.call(rbind, Tukey))
Tukey
colnames(allerie2)<-cls
rownames(allerie2)<-rs
summary(allerie2)
dim(allerie2)


###### Heatmap and clutering---------------------

library(pvclust)
library(dendextend)
library(vegan)
cb<-pvclust(allerie2, method.dist="euclidean", method.hclust="ward.D", nboot=1000)
ca<-hclust(dist(allerie2, method="manhattan"), method="ward.D2")
summary(cb$hclust$order)


par(mfrow=c(2,1), oma = c(0.5, 4.5, 0.75, 0.1), mar = c(4, 0.5, 0.1, 0.1))

cb$hclust$order
plot(cb, hang=0.1)
ords<-c(9, 10, 12, 11, 14, 8, 13, 5, 1, 3, 4, 6, 2, 7, 19, 20, 16, 17, 18, 15, 29, 21, 22, 24, 27, 28, 23, 25, 26)
#dd <- reorder(as.dendrogram(cb$hclust), ords, agglo.fun=mean)
dd<-flip_leaves(as.dendrogram(cb$hclust), c(15), c(29))
df<-flip_leaves(dd, c(24,27,25,26,23,28), c(19,20,16,17,18,21,29,15,22))
da<-flip_leaves(df, c(1,6,3,4), c(5,8,13))
db<-flip_leaves(da, c(5), c(8,13))
dc<-flip_leaves(db, c(25,26,23,28), c(24,27))
dg<-flip_leaves(dc, c(25,26), c(23,28))
eg<-flip_leaves(dg, c(19,20,16,17,18), c(21,29,15,22))
ef<-flip_leaves(eg, c(21), c(29,15,22))
plot(ef)
order.dendrogram(ef)
pp<-flip_leaves(ef, c(9, 10, 11, 12, 14, 8, 13, 5, 1, 6, 3, 4, 2, 7), c(23, 28, 25, 26, 24, 27, 29, 15, 22, 21, 19, 20, 16, 17, 18))
plot(pp)


plot(ef)
pvrect(ef$hclust)

pdf(file = "heatmap-char.pdf", width = 7, height = 5, bg="transparent")
col_breaks = c(seq(-3,-0.2,length.out=115), seq(-0.1,0.1,length.out=275), seq(0.2,3,length.out=115))
my_palette <- colorRampPalette(c("#d73027","black","#2166ac"))(length(col_breaks)-1)
par(oma=c(0.1,0.1,0.1,0.1))
heatmap.2(allerie2, col=my_palette, breaks=col_breaks, margins = c(5,5), density.info="none", trace="none", dendrogram = "column", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F,na.color = "pink", Colv=as.dendrogram(pp), Rowv=as.dendrogram(ca), labRow=ca$labels, scale="none", key=F, key.xlab="Z-score", key.title=NA, cexRow=0.000075, cexCol=0.00000075, offsetRow=3, srtRow=180, lmat=rbind(c(4,3), c(2,1)), lhei=c(1.55, 3.5), lwid=c(0.5, 8.5))
dev.off()

######### data normalization and MDS plot ---------------

har<-t(periods[,9:21])
colnames(har)<-periods$local
stat<-har
rs<-rownames(stat)
cls<-colnames(stat)
Boo<-dim(stat)[[1]]
datum <- list()
for (i in 1:Boo)
{	
	H<-rownames(stat)[i]
    #datum[[H]]<-stat[i,]
    #datum[[H]]<-((as.numeric(stat[i,]) - mean(as.numeric(stat[i,]), na.rm=TRUE))/sd(as.numeric(stat[i,]), na.rm=TRUE))
    datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
}
allerie2<-data.matrix(do.call(rbind, datum))
colnames(allerie2)<-cls
rownames(allerie2)<-rs
summary(allerie2)
dim(allerie2)
Acrodata1<-data.frame(periods[,1:8], t(allerie2))
summary(Acrodata1)
dim(Acrodata1)
acfac<-Acrodata1[,1:8]
acdata.mds<-metaMDS(log(Acrodata1[,9:21]+1), distance='euclidean', k=2, trymax=40, autotransform=FALSE, na.rm=TRUE)
acdata.mds2<-metaMDS(log(Acrodata1[,9:21]+1), previous.best = acdata.mds, distance='euclidean', k=2, trymax=40, autotransform=FALSE, na.rm=TRUE)
acord <- acdata.mds2$points
Acroneword <- cbind(acfac, acord)
summary(Acroneword)
acell<-paste(Acroneword$X1, Acroneword$local, sep="")
acsym <- c(1,19)
acroda<-Acrodata1[,9:21]
summary(acroda)
acrovec<-envfit(acdata.mds2$points, acroda, perl = 999, na.rm=TRUE)
scores(acrovec, "vectors")
acrovecT<-acrovec
Acro<-Acroneword

pdf(file = "MDSplot.pdf", width = 6, height = 6, bg="transparent")
par(mar = c(2.2, 0.5, 0.4, 3.7))
plot(Acro[,9], Acro[,10], col = c("#1f78b4", "#e31a1c", "#ff7f00", "#b2df8a", "#33a02c")[as.factor(Acro$X4)], bg = c("#1f78b4", "#e31a1c", "#ff7f00", "#b2df8a", "#33a02c")[as.factor(Acro$X4)], pch = c(2,2,0,0,0)[as.factor(Acro$X4)], cex = 2, xlab='', ylab='', axes=FALSE, lwd=3, ylim=c(-0.4,0.4), xlim=c(-0.6,0.6))
axis(4, at=seq(-0.4,0.4,0.2), col.axis="black", las=2, cex.axis=1.7)
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
grid(col = "grey")
box(col = "black")
plot(acrovecT, p.max=0.01, col="black", cex=0.0000005, axes=F)
dev.off()

########ANOSIM##############

C15a<-subset(Acrodata1, X4=="C15a")
C15b<-subset(Acrodata1, X4=="C15b")[-4,]
C15c<-subset(Acrodata1, X4=="C15c")
C15d<-subset(Acrodata1, X4=="C15d")
C15e<-subset(Acrodata1, X4=="C15e")

AvB<-rbind(C15a, C15b)
AvB<-droplevels(AvB)
AvC<-rbind(C15a, C15c)
AvC<-droplevels(AvC)
AvD<-rbind(C15a, C15d)
AvD<-droplevels(AvD)
AvE<-rbind(C15a, C15e)
AvE<-droplevels(AvE)
BvC<-rbind(C15b, C15c)
BvC<-droplevels(BvC)
BvD<-rbind(C15b, C15d)
BvD<-droplevels(BvD)
BvE<-rbind(C15b, C15e)
BvE<-droplevels(BvE)
CvD<-rbind(C15c, C15d)
CvD<-droplevels(CvD)
CvE<-rbind(C15c, C15e)
CvE<-droplevels(CvE)
DvE<-rbind(C15d, C15e)
DvE<-droplevels(DvE)

TAvB<-anosim((AvB[,9:21]+1), AvB$X4, permutations = 9999)
TAvC<-anosim((AvC[,9:21]+1), AvC$X4, permutations = 9999)
TAvD<-anosim((AvD[,9:21]+1), AvD$X4, permutations = 9999)
TAvE<-anosim((AvE[,9:21]+1), AvE$X4, permutations = 9999)
TBvC<-anosim((BvC[,9:21]+1), BvC$X4, permutations = 9999)
TBvD<-anosim((BvD[,9:21]+1), BvD$X4, permutations = 9999)
TBvE<-anosim((BvE[,9:21]+1), BvE$X4, permutations = 9999)
TCvD<-anosim((CvD[,9:21]+1), CvD$X4, permutations = 9999)
TCvE<-anosim((CvE[,9:21]+1), CvE$X4, permutations = 9999)
TDvE<-anosim((DvE[,9:21]+1), DvE$X4, permutations = 9999)

A<-data.frame("comp"="AvsB", "stat"=TAvB$stat,"sig"=TAvB$signif)
B<-data.frame("comp"="AvsC", "stat"=TAvC$stat,"sig"=TAvC$signif)
C<-data.frame("comp"="AvsD", "stat"=TAvD$stat,"sig"=TAvD$signif)
D<-data.frame("comp"="AvsE", "stat"=TAvE$stat,"sig"=TAvE$signif)
E<-data.frame("comp"="BvsC", "stat"=TBvC$stat,"sig"=TBvC$signif)
FF<-data.frame("comp"="BvsD", "stat"=TBvD$stat,"sig"=TBvD$signif)
G<-data.frame("comp"="BvsE", "stat"=TBvE$stat,"sig"=TBvE$signif)
H<-data.frame("comp"="CvsD", "stat"=TCvD$stat,"sig"=TCvD$signif)
I<-data.frame("comp"="CvsE", "stat"=TCvE$stat,"sig"=TCvE$signif)
J<-data.frame("comp"="DvsE", "stat"=TDvE$stat,"sig"=TDvE$signif)
rbind(A,B,C,D,E,FF,G,H,I,J)

write.csv(rbind(A,B,C,D,E,FF,G,H,I,J), "anosim.csv")

########FIGURE2##############
########FIGURE2##############
########FIGURE2##############
########FIGURE2##############
########FIGURE2##############
########FIGURE2##############

####### anosim across treatments 

periods<-read.csv("Porites-GCB-21.csv")[,-1]
summary(periods)
dim(periods)
periods<-periods[,1:21]

speci<-c("C15b", "C15a", "C15c", "C15d", "C15e")
ler2<-list()
for(y in 1:5)
{
	Acro1<-droplevels(subset(periods, X4==speci[y]))
	har<-t(Acro1[,9:ncol(Acro1)])
	colnames(har)<-Acro1$local
	stat<-har
	rs<-rownames(stat)
	cls<-colnames(stat)
	Boo<-dim(stat)[[1]]
	datum <- list()
	atum <- list()
	for (i in 1:Boo)
	{		
			if (as.matrix((shapiro.test(stat[i,]))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
				atum[[H]]<-t.test(stat[i,]~Acro1$X1)[[3]]
			}
			else if (as.matrix((shapiro.test(log(stat[i,]+1)))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
				atum[[H]]<-t.test(log(stat[i,]+1)~Acro1$X1)[[3]]
			}
			else if (as.matrix((shapiro.test(sqrt(stat[i,]+1)))[[2]])>0.05)
			{
				H<-rownames(stat)[i]
				datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
				atum[[H]]<-t.test(sqrt(stat[i,]+1)~Acro1$X1)[[3]]
			}
			else
			{	
				H<-rownames(stat)[i]
				datum[[H]]<-((as.numeric(stat[i,])/range(as.numeric(stat[i,]), na.rm=TRUE)[2]))
				atum[[H]]<-wilcox.test(log(stat[i,]+1)~factor(Acro1$X1), p.adjust="none", paired=FALSE)[[3]]
			}
	}
	allerie2<-data.matrix(do.call(rbind, datum))
	ler2[[speci[y]]]<-(data.matrix(do.call(rbind, atum)))
}

perc<-periods
percs<-na.omit(perc)
Nikki<-aggregate(percs[,9:ncol(percs)], by=list(Coral=percs$X4, Temp=percs$X1), FUN=mean)
Nikk<-as.matrix(Nikki[,3:15])
allerie3<-log2(Nikk[6:10,1:13]/Nikk[1:5,1:13])

Stats<-do.call(cbind, ler2)
colnames(Stats)<-c("C15b", "C15a", "C15c", "C15d", "C15e")
Stats

for(L in 1:5)
{
	ler<-data.frame(ler2[L])
	for(I in 1:nrow(ler))
	{
		if(ler[I,1]>=0.05)
		{
			allerie3[L,I]<-NA
		}
	}
}

pdf(file = "heatmap-C15.pdf", width = 9.5, height = 3.5, bg="transparent")
col_breaks = c(seq(-1.25,-0.3,length.out=115), seq(-0.2,0.2,length.out=275), seq(0.3,1.25,length.out=115))
my_palette <- colorRampPalette(c("#d73027","white","#2166ac"))(length(col_breaks)-1)
par(oma=c(0.1,0.1,0.1,0.1))
heatmap.2(allerie3, col=my_palette, breaks=col_breaks, margins = c(1,1), density.info="none", trace="none", dendrogram = "none", symm=FALSE, symkey=FALSE, symbreaks=FALSE, sepwidth=F,na.color = "grey", Colv=F, Rowv=F, scale="none", key=FALSE, key.xlab="Z-score", key.title=NA, cexRow=0.000002, offsetRow=0, srtRow=0, cexCol=0.000000075, srtCol=90, lmat=rbind(c(4,3), c(2,1)), lhei=c(0.1, 3.5), lwid=c(0.5, 8.5))
dev.off()

#######

periods<-read.csv("Porites-GCB-21.csv")[,-1]
summary(periods)

C15a<-droplevels(subset(periods, X4==speci[1]))
mC15a<-aggregate(C15a[,9:ncol(C15a)], by=list(C15a$X1), FUN=mean, na.rm=TRUE)
sC15a<-aggregate(C15a[,9:ncol(C15a)], by=list(C15a$X1), FUN=sd, na.rm=TRUE)
C15b<-droplevels(subset(periods, X4==speci[2]))
mC15b<-aggregate(C15b[,9:ncol(C15b)], by=list(C15b$X1), FUN=mean, na.rm=TRUE)
sC15b<-aggregate(C15b[,9:ncol(C15b)], by=list(C15b$X1), FUN=sd, na.rm=TRUE)
C15c<-droplevels(subset(periods, X4==speci[3]))
mC15c<-aggregate(C15c[,9:ncol(C15c)], by=list(C15c$X1), FUN=mean, na.rm=TRUE)
sC15c<-aggregate(C15c[,9:ncol(C15c)], by=list(C15c$X1), FUN=sd, na.rm=TRUE)
C15d<-droplevels(subset(periods, X4==speci[4]))
mC15d<-aggregate(C15d[,9:ncol(C15d)], by=list(C15d$X1), FUN=mean, na.rm=TRUE)
sC15d<-aggregate(C15d[,9:ncol(C15d)], by=list(C15d$X1), FUN=sd, na.rm=TRUE)
C15e<-droplevels(subset(periods, X4==speci[5]))
mC15e<-aggregate(C15e[,9:ncol(C15e)], by=list(C15e$X1), FUN=mean, na.rm=TRUE)
sC15e<-aggregate(C15e[,9:ncol(C15e)], by=list(C15e$X1), FUN=sd, na.rm=TRUE)

pdf(file = "Density.pdf", width = 7, height = 2, bg="transparent")
par(mfrow=c(1,5), oma = c(0.5, 4.5, 0.75, 0.1), mar = c(1, 0.5, 0.1, 0.1))
colss<-c("#e21f26", "#2078b4", "#f57f20", "#b4d88a", "#34a048")
C15t<-list(mC15b, mC15a, mC15c, mC15d, mC15e)
C15s<-list(sC15b, sC15a, sC15c, sC15d, sC15e)
for(G in 1:5)
{	sC15<-C15s[[G]]
	mC15<-C15t[[G]]
	format(mC15, scientific=T)
	bar<-barplot(mC15$Density, ylim=c(0, 1.75), axes=FALSE, col = colss[G], density=c(1000,15))	
	if(G==1)
	{axis(2, at=, col.axis="black", las=2, cex.axis=2)}
	arrows(bar, (mC15$Density+sC15$Density), bar, (mC15$Density-sC15$Density), lwd=1.0, angle=90, code=3, length=0.025)
	box(col="black")
	par(new=FALSE)
}
dev.off()

#################
#################
#################
#################

PAM <- read.csv("PAM.csv",header=TRUE,sep=",")
summary(PAM)

aPAM<-aggregate(PAM[,6:12], by=list(PAM$Symbiont, PAM$TEMP), FUN=mean)
sPAM<-aggregate(PAM[,6:12], by=list(PAM$Symbiont, PAM$TEMP), FUN=sd)
time<-data.frame("Group.1"="tim1", "Group.2"="tim2", "day1"=1,"day3"=3,"day5"=5,"day7"=7,"day9"=9,"day11"=11,"day13"=13)
aaPAM<-data.frame(rbind(time, aPAM))
speci<-c("C15b", "C15a", "C15c", "C15d", "C15e")
colss<-c("#e21f26", "#2078b4", "#f57f20", "#b4d88a", "#34a048")
sexy<-c(24,24,22,22,22)

pdf(file = "PAM.pdf", width = 2.5, height = 8, bg="transparent")
par(mfrow=c(5,1), oma = c(3, 0.1, 0.1, 0.1), mar = c(1, 3.5, 0.1, 1.1))
for(y in 1:5)
{
	mSingle<-subset(aaPAM, Group.1==speci[y] | Group.2=="tim2")
	sSingle<-subset(sPAM, Group.1==speci[y])
	par(new=FALSE)
	plot(t(mSingle[1,3:9]),t(mSingle[2,3:9]), pch=sexy[y], col=colss[y], bg=colss[y], xlim=c(1,16), ylim=c(0.45,0.7), axes=FALSE, xlab='', ylab='', cex=1.75)
	lines(t(mSingle[1,3:9]),t(mSingle[2,3:9]), col=colss[y], lwd=1)
	epsilon = 0.0
	segments(t(mSingle[1,3:9]), t(mSingle[2,3:9]-sSingle[1,3:9]), t(mSingle[1,3:9]), t(mSingle[2,3:9]+sSingle[1,3:9]), lwd=1, col=colss[y])
	par(new=TRUE)
	plot(t(mSingle[1,3:9]),t(mSingle[3,3:9]), pch=sexy[y], col=colss[y], bg="white", xlim=c(1,16), ylim=c(0.4,0.7), axes=FALSE, xlab='', ylab='', cex=1.75)
	lines(t(mSingle[1,3:9]),t(mSingle[3,3:9]), col=colss[y], lwd=1)
	epsilon = 0.0
	segments(t(mSingle[1,3:9]), t(mSingle[3,3:9]-sSingle[2,3:9]), t(mSingle[1,3:9]), t(mSingle[3,3:9]+sSingle[2,3:9]), lwd=1, col=colss[y])
	axis(2, at=seq(0.45,0.65,0.1), col.axis="black", las=2, cex.axis=1.5)
	#axis(1, at=, col.axis="black", las=1,cex.axis=2)
	grid(col = "grey")
	box(col = "black")
}
axis(1, at=, col.axis="black", las=1,cex.axis=1.5)
dev.off()

##################
##################
##################
##################
##################

YA<-data.frame(na.omit(cbind("Photo"=periods$Photo*100000000, "Density"=periods$Density)))
summary(YA)
YA1<-YA[order(YA$Density),]
Photo<-YA1$Photo
Density<-YA1$Density
try(nls_fit <- nls(Photo ~ a+b*Density^(-c), start = list(a=1, b=1.38, c=0.2),trace = TRUE,control = list(maxiter = 5000), algorithm = "port"), silent = FALSE)
mode(nls_fit)
RR<-cor(YA$Photo,predict(nls_fit))
R<-round(RR, digits=3)
coef <- summary(nls_fit)$coefficients
Aest <- coef[1,1]
Best <- coef[2,1]
Cest <- coef[3,1]
Yest<-Aest+Best*Density^(-Cest)
plot(Density, Photo)
lines(Density, Yest, col = "Black", lty=2, lwd=3)

pdf(file = "Den-Photo.pdf", width = 4, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(3.1, 3.1, 0.1, 0.1), mar = c(0.1, 0.1, 0.1, 0.1))
plot(periods$Density, periods$Photo*100000000, pch=c(1,19)[as.factor(periods$X1)], col=c("#2078b4", "#e21f26", "#f57f20", "#b4d88a", "#34a048")[as.factor(periods$X4)], cex=0.55+(periods$Chla/4), axes=FALSE, lwd=2.5)
lines(Density, Yest, col = "Black", lty=1, lwd=2)
axis(1, at=, col.axis="black", las=1,cex.axis=1.35)
axis(2, at=, col.axis="black", las=1,cex.axis=1.35)
#text(1.15, 0.82, labels=paste("Correlation Coef = ",R,sep=""), cex=1.25)
#text(1.33, 0.76, labels=paste("P-value < 0.0001"), cex=1.25)
box(col="black")
grid(col="grey")
dev.off()


FU<-data.frame(na.omit(cbind("Photo"=periods$Photo*100000000, "NPQ"=periods$NPQ)))
a<-summary(lm(Photo ~ NPQ, data=FU))[[4]][1,1]
b<-summary(lm(Photo ~ NPQ, data=FU))[[4]][2,1]
pdf(file = "Photo-NPQ.pdf", width = 4, height = 4, bg="transparent")
par(mfrow=c(1,1), oma = c(3.1, 3.1, 0.5, 0.5), mar = c(0.1, 0.1, 0.1, 0.1))
plot(periods$Photo*100000000, periods$NPQ, pch=c(1,19)[as.factor(periods$X1)],  col=c("#2078b4", "#e21f26", "#f57f20", "#b4d88a", "#34a048")[as.factor(periods$X4)], cex=0.55+(periods$Chla/4), axes=FALSE, lwd=2.5)
axis(1, at=, col.axis="black", las=1,cex.axis=1.35)
axis(2, at=, col.axis="black", las=1,cex.axis=1.35)
box(col="black")
grid(col="grey")
dev.off()

##################