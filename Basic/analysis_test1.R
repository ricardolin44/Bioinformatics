library(tidyverse)
library(corrr)
library(factoextra)
library(pheatmap)

ds4B <- read.csv('DATA_FSB_SET_4B.csv', row.names = 1)
if (any(!complete.cases(ds4B)))
{
  ds4B.nona <- ds4B %>% drop_na()
}
ds4B.nona <- as.matrix(ds4B.nona)
View(ds4B.nona)

hist(ds4B.nona, xlab='Normal Expression values')
ds4B.log <- log1p(ds4B.nona)
hist(ds4B.log, xlab='Normal Expression values')


venus.idx <- grep('Venus', colnames(ds4B.log))
earth.idx <- grep('Earth', colnames(ds4B.log))

qqnorm(ds4B.log[2,venus.idx])
qqline(ds4B.log[2,venus.idx], col='steelblue',lwd=2)

shapiro.test(ds4B.log[2,venus.idx])$p.value
shapiro.test(ds4B.log[2,earth.idx])$p.value

p.values <- sapply(1:nrow(ds4B.log),
  function(x) t.test(ds4B.log[x,earth.idx],ds4B.log[x,venus.idx])[c('p.value')]
)
p.values <- as.data.frame(p.values)
table(p.values<0.05)  #We can deny null hypothesis which means that we know that earth and venus are different

#plot heatmap
pval.idx <- p.values<0.05
pheatmap(ds4B.log[pval.idx,])

par(mfrow=c(1,2))
ds4B.log.pca <- prcomp(t(ds4B.log), scale=T, center=T)
plot(ds4B.log.pca,
     xlab='Dimension',
     main='Screenplot')

#Amount of explained variance
cp <- cumsum(ds4B.log.pca$sdev^2/sum(ds4A.log.pca$sdev^2))
plot(cp,
     xlab='PC#',
     ylab='Amount of explained variance',
     main='Cumulative variance plot')
#or
fviz_eig(ds4B.log.pca, addlabels=T, ylim=c(0,70))


#Plotting pca
fviz_pca_biplot(ds4B.log.pca,
                label='var')
                #habillage=ds4B.log$class

fviz_pca_ind(ds4B.log.pca,
             axes=c(3,1),
             geom=c('point'))



