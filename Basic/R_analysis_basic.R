names = c('John', 'Smith', 'John')
ages = c(5, 6, 5)
data = data.frame(names,ages)
View(data)

distinct(data)

data()
View(starwars)

library(tidyverse)

starwars %>% 
  filter(height > 150 & mass < 200) %>% 
  mutate(height_in_meter = height/100) %>% 
  select(height_in_meter, mass) %>% 
  arrange(-mass) %>% 
  View()
  
glimpse(msleep)
names(msleep)
unique(msleep)

missing <- !complete.cases(msleep)
msleep[missing,]

starwars %>% 
  select(ends_with('color'))

df <- starwars
df$sex <- as.factor(df$sex)

df <- df %>% 
  mutate(sex = factor(sex,
                      levels= c('male', 'female', 'hermaphroditic', 'none')))
levels(df$sex)

starwars %>% 
  select(sex) %>% 
  mutate(sex = recode(sex,
                      'male'='man',
                      'female'='woman'))
starwars %>% 
  select(name, height) %>% 
  mutate(tallness = 
           if_else(height<100,
                   'short',
                   'tall'))

library(gapminder)
View(gapminder)

data <- select(gapminder, country, year, lifeExp)
wide_data <- data %>% 
  pivot_wider(names_from = year, values_from = lifeExp)

long_data <- wide_data %>% 
  pivot_longer(2:13, names_to = 'year', values_to = 'lifeExp')

msleep %>% 
  drop_na(vore) %>% 
  group_by(vore) %>% 
  summarise(Lower = min(sleep_total),
            Average = mean(sleep_total),
            Max = max(sleep_total),
            Difference = max(sleep_total) - min(sleep_total)) %>% 
  arrange(Average) %>% 
  View()

starwars
ggplot(data=starwars,
       mapping=aes(x=gender))+
  geom_bar()

starwars %>% 
  drop_na(height) %>% 
  ggplot(mapping = aes(x=height))+  #mapping is not essential
  geom_boxplot(fill='steelblue')+
  theme_bw()+
  labs(title='Boxplot of height',
       x='Height of chars')

#Density plots
starwars %>% 
  drop_na(height) %>% 
  filter(sex %in% c('male','female', 'none')) %>% 
  ggplot(aes(height,
             color=sex,
             fill=sex))+
  geom_density(alpha=0.2)+
  theme_bw()

starwars %>% 
  drop_na(height) %>% 
  filter(mass<200) %>% 
  ggplot(aes(height, mass, color=sex))+ 
  geom_point(size=3,alpha=0.8)+
  geom_smooth()+   #add one layer
  facet_wrap(~sex)+  
  theme_bw()+
  labs(title='Boxplot of height',
       x='Height of chars')

library(tidyverse)
flower <- iris %>%
  mutate(Size = cut(Sepal.Length,
                    breaks = 3,
                    labels = c('small', 'medium', 'large'))) %>% 
  select(Species, Size)
View(flower)

flower %>% 
  select(Size) %>% 
  table() %>% 
  chisq.test()

flower %>% 
  table() %>% 
  chisq.test()

model2 <- lm(dist~speed, data=cars)
cars %>% 
  lm(dist~speed, data=.) %>%
  ggplot(aes(speed,dist))+
  geom_point(size=3,alpha=0.5,color='#4073FF')+
  geom_smooth(method = 'lm', se=FALSE, color='#FF1234')+
  theme_bw()+
  labs(title='The relationship between speed and stopping distance',
       x='Speed of car',
       y = 'Distance taken to stop')+
  annotate('Text', x=5, y=110,
           label=paste('Slope:', round(coef(model2)[2]),2),
           color='black', size=5)+
  annotate('Text', x=5, y=100,
            label=paste('Intercept:', round(coef(model2)[1]),2),
            color='black', size=5)


coef(model2)
for (a in coef(summary(model2)))
{
  print(a)
}

coef(summary(model2))[,2][1]

Puromycin %>% 
  ggplot(aes(conc,rate, color=state)) +
  geom_point(size=3, alpha=0.7)+
  geom_smooth(method='loess', se=FALSE)+
  facet_wrap(~state)+
  theme_bw()+
  labs(title='The relationship between medicine concentration and its rate',
       x='Conc',
       y = 'Rate')

ds2 <- read.csv('DATA_FSB_SET_2.csv', row.names=1)
ds2 <- ds2[complete.cases(ds2),]
View(ds2)

#set the plot dimensions here
options(repr.plot.width = 4.5, repr.plot.height = 3)
library(tidyverse)
ggplot(ds2, aes(x=Weight, y=LDL)) +
  geom_point(aes(col=Hospital_Visits), size=3)+
  xlim(c(0,350))+
  ylim(c(0,300))

idx <- grep('Color_House', colnames(ds2))
ds2.onlyNumerical <- ds2[,-idx]

pca <- prcomp(ds2.onlyNumerical, scale=TRUE, center=T)
summary(pca)

#Visualize Eigenvalues
options(repr.plot.height = 6, repr.plot.width = 5)
plot(pca,
     xlab='Dimension',
     main='Screen plot')

cp <- cumsum(pca$sdev^2/ sum(pca$sdev^2))
plot(cp,
     xlab='PC#',
     ylab='Amount of explained variance',
     main='Cumulative variance plot')

#pca$rotation 
#Visualize individual patienst using only the first two components
plot(pca$x[,1], pca$x[,2],
     xlab='PC 1',
     ylab='PC 2',
     main='PCA across Patients')

library('factoextra')
options(repr.plot.height = 10, repr.plot.width = 10)
fviz_pca_var(pca, geom=c('point','text'))
fviz_pca_biplot(pca, geom=c('point','text'))

#constructing a distance matrix between all sample combinations
ds2.onlyNumerical.sc <- scale(ds2.onlyNumerical, scale=T, center=T)
dist.by.variable <- dist(t(ds2.onlyNumerical.sc))
dist.by.patient <- dist(ds2.onlyNumerical.sc)    #use in heatmaps

#Hierarchical clustering
dt.clust.by.var <- hclust(dist.by.variable)
dt.clust.by.patient <- hclust(dist.by.patient)  #for heatmaps

#Plot a dendrogram
options(repr.plot.height = 5, repr.plot.width = 5)
plot(dt.clust.by.var)

#plot heatmap
library(pheatmap)
options(repr.plot.height = 15, repr.plot.width = 6)
pheatmap(ds2.onlyNumerical.sc,
         show_rownames = F,
         show_colnames = T,
         cluster_rows = dt.clust.by.patient,
         cluster_cols = dt.clust.by.var)

ds2 <- read.csv('DATA_FSB_SET_2.csv', row.names=1)
ds2 <- ds2[complete.cases(ds2),]
ds1 <- read.csv('DATA_FSB_SET_1.csv', row.names=1)
ds1 <- ds1[complete.cases(ds1),]
View(ds2)

wt.vs.ldl.lm <- lm(Weight ~ LDL, ds1)
plt <- ds1 %>% 
  ggplot(aes(LDL, Weight))+
  geom_point()+
  geom_smooth(method='lm', se=F)

mean(wt.vs.ldl.lm$residuals)  #need to be near 0
hist(wt.vs.ldl.lm$residuals)
wt.vs.ldl.lm
plot(wt.vs.ldl.lm,1)  # Residuals vs Fitted
plot(wt.vs.ldl.lm,2)  #Q-Q Residual Plot

require(lawstat)
?levene.test

library(tidyverse)
library(corrr)
ds3 <- read.csv('DATA_FSB_SET_3.csv',row.names=1)
summary(ds3)

ds3_corr <- ds3 %>% 
  select(-Planet) %>% 
  select(-Planet2) %>% 
  cor(method='pearson')

boxplot(ds3_corr)

ds3.df <- ds3 %>% 
  mutate(Planet=factor(Planet)) %>% 
  mutate(Planet2=factor(Planet2))
View(ds3.df)

ds3.df %>% 
  group_by(Planet) %>% 
  summarise(
    count_planet=n(),
    mean_LDL = mean(LDL_levels, na.rm=T),
    sd_LDL = sd(LDL_levels, na.rm=T)
  )


#Visualizing density whether it's normalized or not
planet <- 'Venus'
plt <- ds3 %>% 
  filter(Planet==planet) %>% 
  ggplot()+
  geom_histogram(aes(x=LDL_levels,y=after_stat(density), fill=Planet), binwidth=5)+
  xlab('LDL [mg/dl]')+
  stat_function(fun=dnorm,
                args=list(mean= mean(ds3$LDL_levels[ds3$Planet==planet]),
                          sd = sd(ds3$LDL_levels[ds3$Planet==planet])))+
  scale_fill_manual(values=c('#E69F00','#56B4E9','#999999'))
plt

#Assumption 1, all are normally distributed
#divide into 2 
options(repr.plot.width=7, repr.plot.height=4)
par(mfrow=c(1,2), bg='white')

#For Earth
qqnorm(ds3$LDL_levels[ds3$Planet=='Earth'], pch=3, frame=F, main='Earth')
qqline(ds3$LDL_levels[ds3$Planet=='Earth'], col='steelblue', lwd=2)

#For Venus
qqnorm(ds3$LDL_levels[ds3$Planet=='Venus'], pch=3, frame=F, main='Venus')
qqline(ds3$LDL_levels[ds3$Planet=='Venus'], col='steelblue', lwd=2)

#Shapiro Wilk Normality Test -> null hypothesis = normally distributed
shapiro.test(ds3$LDL_levels[ds3$Planet=='Earth'])
shapiro.test(ds3$LDL_levels[ds3$Planet=='Venus'])


#Assumption 2: each group has the same variance
#Bartlett test -> null hypothesis = same variance in each group
ldl_vs_planet.varTest <- bartlett.test(LDL_levels ~ Planet, data =ds3)
ldl_vs_planet.varTest


#Assumption 3: Independence(No bias in our selection)
#Running t-test 
#t-test (need to be caution with Variance is different or indifferent -> different function will be used)
ldl.vs.planet.ttest <- t.test(LDL_levels ~ Planet, data=ds3)    #Var is different
ldl.vs.planet.ttest
signif(ldl.vs.planet.ttest$p.value, digits=3)

ldl.vs.planet.ttest.2 <- t.test(LDL_levels ~ Planet, data=ds3, var.equal=T)    #Var is equal
ldl.vs.planet.ttest.2
signif(ldl.vs.planet.ttest.2$p.value, digits=3)

ds3 %>% 
  group_by(Planet2) %>% 
  summarise(
    count_planet=n(),
    mean_LDL = mean(LDL_levels, na.rm=T),
    sd_LDL = sd(LDL_levels, na.rm=T)
  )

plt <- ggplot(ds3.df, aes(Planet2, LDL_levels, fill=Planet2))+   #theres no diff with ds3.df(factorized) and ds3 in the results
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
plt

#One way anova test -> can only output 1 Pr, can't output relation hypothesis between groups
anova_one_way <- aov(Exercise ~ Planet2, data=ds3)
summary(anova_one_way)

#Tukey's range test -> similar to t-test and used with anova, it can compare many means datas to each other 
TukeyHSD(anova_one_way, ordered=T)
plot(TukeyHSD(anova_one_way, ordered=T))
