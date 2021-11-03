#call libraries
library(factoextra)
library(reshape2)
library(cluster)
library(ggplot2)
library(FactoMineR)
library(corrplot)

#choose variables from my original dataset
mydata <- read.csv2(file ="C:\\Users\\DELL\\Downloads\\ESS1.csv")
df <- mydata [1:10]
df <- na.omit(df)

#create dataset grouping by countries and add new variables
a<- aggregate(df[-1], by = list(df$ï..country), FUN =mean)
a$freedom <- c(93, 96, 96, 94, 97, 90, 100, 90, 93, 69, 97, 90, 98, 100, 82, 96, 100, 95)
a$density <- c(106,376,206,232, 135,92, 16,118, 272, 105, 69, 200,421, 14, 123,112,23,102)
a$TAI <- c (0.617, 0.604, 0.813, 0.658, 0.666, 0.534, 0.633, 0.622, 0.546, 0.516, 0.682, 0.507, 0.745, 0.626, 0.626, 0.467, 0.685, 0.556)
dati1<-a[,-1]

#normalize variables
dati1<-scale(dati1[,1:12])
rownames(dati1)<-c("Austria", "Belgium", "Switzerland", "Germany", "Denmark", "Spain", "Finland", "France", "UK", "Hungary", "Ireland", "Italy", "Netherland", "Norway", "Poland", "Portugal", "Sweden", "Slovenia")

#choose the distances
d1 <- dist(dati1, method="euclidean", diag=F, upper=F)
d2<- dist (dati1, method = "manhattan", diag=F, upper= F)

#function to generate the agglomeration program
agglo <- function(hc){
  data.frame(row.names=paste0("Cluster",seq_along(hc$height)),
             height=hc$height,
             components=ifelse(hc$merge<0, 
                               hc$labels[abs(hc$merge)], paste0("Cluster",hc$merge)),
             stringsAsFactors=FALSE) }

#with complete linkage
h1 <- hclust(d1, method="complete"); h1
agglo(h1)
plot(h1, main="complete linkage")
complete <- cutree(h1, k=4)
rect.hclust(h1, 4)
h1cluster <- cutree(h1, k=4)
h1cluster 

#with manhattan distance
h12 <- hclust(d2, method="complete"); h12
agglo(h12)
plot(h12, main="complete linkage")
complete <- cutree(h12, k=4)
rect.hclust(h12, 4)
h12cluster <- cutree(h12, k=4)
h12cluster 

#with avg linkage
h2<-hclust(d1,method="average");h2
agglo(h2)
plot(h2, main="average linkage")
average <- cutree(h2, k=4)
rect.hclust(h2, 4)
h2cluster <- cutree(h2, k=4)
h2cluster 

#with manhattan distance
h22 <- hclust(d2, method="average"); h22
agglo(h22)
plot(h22, main="average linkage")
complete <- cutree(h22, k=4)
rect.hclust(h22, 4)
h22cluster <- cutree(h22, k=4)
h22cluster 

#with single linkage
h3<-hclust(d1,method="single");h3
agglo(h3)
plot(h3, main="single linkage")
single<- cutree(h3, k=4)
rect.hclust(h3, 4)
h3cluster <- cutree(h3, k=4)
h3cluster

#with manhattan distance
h32 <- hclust(d2, method="single"); h32
agglo(h32)
plot(h32, main="single linkage")
complete <- cutree(h32, k=4)
rect.hclust(h32, 4)
h32cluster <- cutree(h32, k=4)
h32cluster 

#with ward linkage
h4<-hclust(d1,method="ward.D");h4
agglo(h4)
plot(h4, main="Ward linkage")
ward<- cutree(h4, k=4)
rect.hclust(h4, 4)
h4cluster <- cutree(h4, k=4)
h4cluster
plot(a[-1], col=h4cluster, main="Ward Methods")

#with manhattan distance
h42 <- hclust(d2, method="ward.D"); h42
agglo(h42)
plot(h42, main="ward.D linkage")
complete <- cutree(h42, k=4)
rect.hclust(h42, 4)
h42cluster <- cutree(h42, k=4)
h42cluster 

#add the cluster to the dataset
h4cluster<- cutree(h4, k=4)
dati1 <- as.data.frame(dati1)
a$clu<-h4cluster
dati1$clu<-h4cluster

#means for variables (non normalized)
medie<-aggregate(a[,2:13], list(h4cluster), mean)
medie 

#calculus of R^2 for each variables
mydata<-dati1
R2 <- rep(NA, (ncol(mydata)-1))
for(i in 1:(ncol(mydata)-1)) 
  R2[i] <- anova(aov(mydata[,i] ~ mydata[,ncol(mydata)]))[1,2]/(anova(aov(mydata[,i] ~ mydata[,ncol(mydata)]))[1,2]+anova(aov(mydata[,i] ~ mydata[,ncol(mydata)]))[2,2])
R2
col<-colnames(mydata[-13])
finali<-cbind(col,round(R2,2))

#plots for cluster interpretation
col<-colnames(mydata)
mydataz<-data.frame(scale(mydata))
mydataz$clu<-h4cluster
dati <- melt(mydataz, measure.vars=col)
ggplot(dati[1:216,], aes(x = variable, y = value, color=variable)) +
  geom_boxplot() +
  facet_wrap(~ mydataz$clu) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

#K-Means
dati1.stand <- dati1[,-13] 

#select number of K (Elbow method)
set.seed(123)
wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}
wssplot(dati1.stand)

#select number of K (Silhouette method)
fviz_nbclust(dati1.stand, kmeans,nstart=15, method = "silhouette")

#k = 2
k.means.fit <- kmeans(dati1.stand, 2) 
str(k.means.fit)
clusplot(dati1.stand, k.means.fit$cluster, 
         main='2D representation of the Cluster solution',
         color=TRUE,
         labels=2, lines=0)

set.seed(1233)
final <- kmeans(dati1.stand, 2, nstart = 15)
print(final)
fviz_cluster(final, data = dati1.stand)

#correlation matrix
corrplot(cor(dati1[-13]))

#PCA
datipuliti<- dati1.stand
n <- nrow(datipuliti)
p <- ncol(datipuliti)

#PCA from correlation
rho <- cor(datipuliti)
eigen(rho)
autoval <- eigen(rho)$values
autovec <- eigen(rho)$vectors

#select components
pvarsp = autoval/p
pvarspcum = cumsum(pvarsp)
pvarsp

#scree Diagram
plot(autoval, type="b", main="Scree Diagram", xlab="Components", ylab="Eigenvalue")
abline(h=1, lwd=3, col="red")

#number of components
world.pca<-prcomp(dati1[-13], center = TRUE,scale. = TRUE)
fviz_eig(world.pca,barcolor = "red",
         barfill = "red",geom = c("bar"),addlabels= TRUE ) +labs(title = "Variances - PCA",
                                                                 x = "Principal Components", y = "% of variances") +theme_bw()+theme(plot.title = element_text(hjust = 0.5))
summary(world.pca)

#interpret the components
eigen(rho)$vectors[,1:2]

#component matrix
comp<-round(cbind(-eigen(rho)$vectors[,1]*sqrt(autoval[1]),-eigen(rho)$vectors[,2]*sqrt(autoval[2])),3)
colnames(comp)<-c("Comp1","Comp2")
comp
rownames

#commonality
comunalita<-comp[,1]^2+comp[,2]^2
comp<-cbind(comp,comunalita)
comp<- as.data.frame(comp)
rownames(comp)<-colnames(datipuliti)

#scores
datipuliti.scale <- scale(datipuliti, T, T)
punteggi <- datipuliti.scale%*%autovec[,1:2]

#standardized scores
punteggiz<-round(cbind(-punteggi[,1]/sqrt(autoval[1]),-punteggi[,2]/sqrt(autoval[2])),2)
plot(punteggiz, main="Score plot",
     xlab="comp1",ylab="comp2")
text(punteggiz, rownames(dati1))
abline(v=0,h=0,col="red")

#loadinngs
plot(comp[,1:2], main="Loadings plot",
     xlab="comp1",ylab="comp2", xlim=range(-1,1))
text(comp, rownames(comp))
abline(v=0,h=0,col="red")

#princomp 
acp<-princomp(datipuliti, cor=T)
summary(princomp(datipuliti, cor=T))

#biplot 
biplot(acp)
fviz_pca_var(world.pca,col.var = "contrib",
             gradient.cols = c("red","orange","blue"),
             repel = TRUE,col.circle = "black",arrowsize = 1,
             labelsize = 0.5,jitter = list(what = "both", width = 1, height = 1) )+theme_bw()+theme(plot.title = element_text(hjust = 0.5))
