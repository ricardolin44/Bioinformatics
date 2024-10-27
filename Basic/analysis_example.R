library(tidyverse)
library(corrr)
library(factoextra)
library(pheatmap)

#Test dataset4 of individuals living in Earth and in Venus
ds4A <- read.csv('./R/DATA_FSB_SET_4A.csv', row.names = 1)
any(!complete.cases(ds4A))
ds4A <- as.matrix((ds4A))
View(ds4A)

hist(ds4A, xlab='Normal Expression Values')
ds4A.log <- log1p((ds4A))
hist(ds4A.log, xlab='Log scaled Normal Expression Values')  #Need to check whether there's minus or pseudo value

data <- ds4A.log %>%
  as.data.frame %>%
  pivot_longer(names_to = 'patient', values_to = 'expression', cols = 1:ncol(ds4A))
View(data)
data$planet <- 'Venus'
data$planet[grep('Earth', data$patient)] <- 'Earth'  #using patient pattern function(grep) to make a new column named planet

plt <- data %>%
  ggplot()+
  geom_histogram(aes(expression, after_stat(density), fill=planet),
                 binwidth = .5, alpha=.6)+
  xlab('Expression [log1p]')
plt  #not very informative because we want to know how each gene is different to each other

#grabbing the index of Venus and Earth respectively from log ds4A data
venus.idx <- grep('Venus', colnames(ds4A.log))
earth.idx <- grep('Earth', colnames(ds4A.log))

#Outputting a qqplot but only for gene3 in venus for test
qqnorm(ds4A.log[3, venus.idx])
qqline(ds4A.log[3, venus.idx], col='steelblue', lwd=2)

#It's difficult to see if there's normality in the graphic so we performed a shapiro test
shapiro.test(ds4A.log[3, venus.idx])
shapiro.test(ds4A.log[3, earth.idx])

#The normality of the data is confirmed but we would like to check the structure of the data in different dimensions
ds4A.log.pca <- prcomp(t(ds4A.log), scale=F, center=F)
options(repr.plot.height = 7, repr.plot.width = 7)
plot(ds4A.log.pca,
     xlab='Dimension',
     main='Screenplot')

#cumulative explained variability plot (roughly confirm that the 1st variance can explain about 98% of the variant)
cp <- cumsum(ds4A.log.pca$sdev^2/sum(ds4A.log.pca$sdev^2))
plot(cp,
     xlab='PC#',
     ylab='Amount of explained variance',
     main='Cumulative variance plot')

col.by.planet <- rep('Earth', ncol(ds4A.log))
col.by.planet[grep('Venus',colnames(ds4A.log))] <- 'Venus'

#Visualize pca
fviz_pca_ind(ds4A.log.pca,
             axes=c(1,2),
             geom=c('point'),
             col.ind = col.by.planet)

#t test 1 gene to see if there's a difference between means in Venus and Earth in Gene1 expression values
t.test(ds4A.log['Gene1', venus.idx], ds4A.log['Gene1',earth.idx])

#do this to all genes with loop function
p.vals <- sapply(1:nrow(ds4A.log),
                 function(i)
                   t.test(ds4A.log[i, earth.idx], ds4A.log[i,venus.idx])[c('p.value')]
                 )
table(p.vals<0.05)

#plot it into heatmap
de.idx <- p.vals < 0.05
options(repr.plot.height = 20, repr.plot.width = 7)
pheatmap(ds4A.log[de.idx,])

#Investigating the variability in genes
#Means vs variance plot -> show that statistics vary between highly and lowly expressed genes
plot(apply(ds4A.log[,earth.idx],1,mean),apply(ds4A.log[,earth.idx],1,var),
     xlab='Mean expression [log]',
     ylab='Expression variance [log]',
     main='Expression Mean vs Variance for Earth samples',
     pch=19)

apply(ds4A.log[,venus.idx],1,mean)
ds4A.log[,venus.idx]

#Volcano plot
fc.log <- -log10(apply(ds4A.log[,venus.idx],1,mean)/apply(ds4A.log[,earth.idx],1,mean))
col.fc <- rep('black', nrow(ds4A.log))
col.fc[p.vals<0.01 & fc.log <0] <- 'red'
col.fc[p.vals<0.01 & fc.log >0] <- 'green'

plot(fc.log, -log10(unlist(p.vals)),
     main='Volcano plot',
     xlab='mean log expression',
     ylab='sd log expression',
     col=col.fc,
     pch=19)
abline(h=-log10(0.01), v=0)
