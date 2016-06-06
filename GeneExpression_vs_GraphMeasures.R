#======================================================================
# Structural connectome analysis - Allen Institute data
#======================================================================

require(Hmisc)
require(psych)
require(reshape)
require(ggplot2)
require(plyr)

##############################
# Getting descriptives
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data = read.csv('Graph_measures_and_expression.csv')
describe(data)

#====================================================================================
# Group by group comparison
data = read.csv('Graph_measures_and_expression.csv')
measures <- c('NodeDegree','Strength','ClustCoeff','BetweenCent','LocEff')
for (measure in measures){
  variable <- c(paste('patient.',measure,sep=""),paste('control.',measure,sep=""))
  data = read.csv('Graph_measures_and_expression.csv')
  data <- data[,c('ZDHHC9',variable[1],variable[2])]
  data <- melt(data,id=c('ZDHHC9'))
  
  figure <- ggplot(data,aes(x=ZDHHC9,y=value,colour=factor(variable))) + 
    geom_point() + 
    xlab('ZDHHC9 expression') + ylab(measure) +    
    stat_smooth(method='lm') 
  show(figure)
  print(summary(aov(data=data,value ~ variable + variable*ZDHHC9)))
  print(rcorr(data$ZDHHC9,data$value))
}

ggplot(data,aes(x=factor(variable),y=value)) + 
  geom_boxplot()

#====================================================================================
# Group differences
data = read.csv('Graph_measures_and_expression.csv')
measures <- c('NodeDegree','Strength','ClustCoeff','BetweenCent','LocEff')

for (measure in measures){
  variable <- c(paste('patient.',measure,sep=""),paste('control.',measure,sep=""))
  data = read.csv('Graph_measures_and_expression.csv')
  data <- data[,c('ZDHHC9',variable[1],variable[2])]
  data$value <- data[[variable[1]]] - data[[variable[2]]]

  figure <- ggplot(data,aes(x=ZDHHC9,y=value)) + 
    geom_point() + 
    stat_smooth(method='lm')
  show(figure)
  results <- t.test(data[[paste('patient.',measure,sep="")]],data[[paste('control.',measure,sep="")]],paired=TRUE)
  print(results)
}

#====================================================================================
# Mean correlation analysis
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data = read.csv('Graph_measures_and_expression.csv')

measures <- c('Strength','ClustCoeff')

for (measure in measures){
  variable <- c(paste('patient.',measure,sep=""),paste('control.',measure,sep=""))
  x <- data$ZDHHC9
  #y <- data[[variable[2]]] - data[[variable[1]]]
  y <- data[[variable[2]]]
    
  figure <- ggplot(data,aes(x=x,y=y)) + 
    geom_point() + 
    stat_smooth(method='lm')
  
  show(figure)
  results <- rcorr(x,y)
  print(results)
}

#====================================================================================
# Individual Correlations
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data = read.csv('rh_NodeStrength.csv')
data <- melt(data,id=c('X','ZDHHC9','normalizedExpression'))
data <- data[ order(data$ZDHHC9), ]

ggplot(data,aes(x=ZDHHC9, y=value, colour=factor(variable))) +
  geom_point() +
  stat_smooth(method='lm',se=FALSE)

setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data = read.csv('rh_NodeStrength.csv')

participants <- c('z1','z2','z3','z4','z5','z6','z8','c1','c2','c3','c5','c6','c7','c8')

for (participant in participants){
  results <- rcorr(data$normalizedExpression,data[[participant]])
  print(results$P[2])
}


#====================================================================================
# Complete model
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data <- read.csv('NodeStrength.csv')
data <- melt(data,id=c('X','ZDHHC9','normalizedExpression'))
data <- rename(data,c("X"="cortical_region","ZDHHC9"="ZDHHC9_expression","value"="graph_measure","variable"="participant"))

summary(aov(data=data, graph_measure ~ ZDHHC9_expression + factor(cortical_region)))
# There is a significant effect of ZDHHC9 expression and cortical region


#====================================================================================
# Regression model

# Complete model
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data <- read.csv('group_mean_graph_measures.csv')
measures <- data[,c('Total.Node.Strength','Total.Clustering.Coefficient','Total.Local.Efficiency','Total.Betweenness.Centrality')]
measures <- measures[complete.cases(measures),]
cor(measures)

## Groups combined
fit <- lm(ZDHHC9.expression ~ Total.Node.Strength + Total.Clustering.Coefficient + Total.Local.Efficiency + Total.Betweenness.Centrality, data=data)
coefficients(fit)
anova(fit)

## Groups separate
fit <- lm(ZDHHC9.expression ~ ZDHHC9.Node.Strength + ZDHHC9.Clustering.Coefficient + ZDHHC9.Local.Efficiency + ZDHHC9.Betweenness.Centrality, data=data)
coefficients(fit)

#====================================================================================
# Structural equation modelling
setwd("/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/")
data <- read.csv('data_for_SEM.csv')
values <- data[,c('NodeStrength','ClustCoeff','BetweenCent','LocEff')]
values <- values[complete.cases(values),]

# Checking normality
shapiro.test(data$NodeStrength)
shapiro.test(data$ClustCoeff)
shapiro.test(data$BetweenCent)
shapiro.test(data$LocEff)
shapiro.test(data$ZDHHC9)

# Applying z-scores
data$ZDHHC9 <- scale(data$ZDHHC9)
data$GAPDH <- scale(data$GAPDH)
data$FMR1 <- scale(data$FMR1)
data$FOXP2 <- scale(data$FOXP2)
data$GRIN2A <- scale(data$GRIN2A)
data$NodeStrength <- scale(data$NodeStrength)
data$ClustCoeff <- scale(data$ClustCoeff)
data$BetweenCent <- scale(data$BetweenCent)
data$LocEff <- scale(data$LocEff)

# Checking for multi-colinearity
cor(values)

# Regression model with gene expression as predictor
summary(lm(data=data,NodeStrength ~ ZDHHC9 + Group))
summary(lm(data=data,ClustCoeff ~ ZDHHC9 + Group))
summary(lm(data=data,BetweenCent ~ ZDHHC9 + Group))
summary(lm(data=data,LocEff ~ ZDHHC9 + Group))

summary(lm(data=data,NodeStrength ~ GAPDH + Group))
summary(lm(data=data,ClustCoeff ~ GAPDH + Group))
summary(lm(data=data,BetweenCent ~ GAPDH + Group))
summary(lm(data=data,LocEff ~ GAPDH + Group))

summary(lm(data=data,NodeStrength ~ GRIN2A + Group))
summary(lm(data=data,ClustCoeff ~ GRIN2A + Group))
summary(lm(data=data,BetweenCent ~ GRIN2A + Group))
summary(lm(data=data,LocEff ~ GRIN2A + Group))

summary(lm(data=data,NodeStrength ~ FMR1 + Group))
summary(lm(data=data,ClustCoeff ~ FMR1 + Group))
summary(lm(data=data,BetweenCent ~ FMR1 + Group))
summary(lm(data=data,LocEff ~ FMR1 + Group))

summary(lm(data=data,NodeStrength ~ FOXP2 + Group))
summary(lm(data=data,ClustCoeff ~ FOXP2 + Group))
summary(lm(data=data,BetweenCent ~ FOXP2 + Group))
summary(lm(data=data,LocEff ~ FOXP2 + Group))



# Plots
require(reshape)
melted_data <- melt(data,id=c('ZDHHC9','GAPDH','FMR1','FOXP2','GRIN2A','X','Group'))
melted_data <- subset(melted_data,melted_data$variable != 'NodeDegree')
melted_data$Group <- factor(melted_data$Group,levels=c(1,0),labels=c('ZDHHC9','control'))
melted_data$variable <- factor(melted_data$variable,levels=c('NodeStrength','ClustCoeff','BetweenCent','LocEff'),
                               labels=c('Node Strength','Clustering Coefficient','Betweenness Centrality','Local Efficiency'))
library(plyr)
melted_data <- rename(melted_data,c("variable"="graph.metric", "value"="graph.value"))
melted_data <- melt(melted_data,id=c('graph.metric','graph.value','Group','X'))
melted_data <- rename(melted_data,c("variable"="gene.name", "value"="gene.expression"))

require(ggplot2)
melted_data <- subset(melted_data,melted_data$graph.metric == 'Clustering Coefficient')
g <- ggplot(data=melted_data,aes(x=gene.expression,y=graph.value,colour=factor(Group),fill=factor(melted_data$Group))) +
  stat_smooth(method='lm') +
  geom_point(size=1.5,alpha=0.9) +
  facet_wrap(~variable,scales='fixed') +
  xlab('gene expression') +
  ylab('z-transformed graph measure') + 
  theme_minimal() + 
  theme(strip.text.x = element_text(size=9, face='bold'),
        axis.text=element_text(size=8),
        legend.position = 'None') +
    facet_wrap(~gene.name)

ggsave(g,file='/Users/joebathelt1/Documents/Other_Projects/zDHHC9/connectivity/results/regression.png',dpi=300,width=90,height=90,unit='mm')

require(psych)
describeBy(data$ClustCoeff,data$Group,mat=TRUE)


require(ggplot2)
melted_data <- subset(melted_data,melted_data$gene.name == 'ZDHHC9')

svg(paste("/Users/joebathelt1/Desktop/ZDHHC9_vs_connectome.svg",sep=""), 
    width=5.5, 
    height=3.5, 
    pointsize=12)

ggplot(data=melted_data,aes(x=gene.expression,y=graph.value,colour=factor(Group),fill=factor(melted_data$Group))) +
    stat_smooth(method='lm') +
    geom_point(size=1.5,alpha=0.9) +
    facet_wrap(~variable,scales='fixed') +
    xlab(expression(paste('z-transformed ',italic('ZDHHC9'),' expression'))) +
    ylab('z-transformed graph measure') + 
    theme_classic() +
    theme(legend.position = 'right',
          text=element_text(size=12, family="CMU Serif"),
          legend.title = element_text(size=12, face='bold'),
          axis.text = element_text(size=12),
          axis.title = element_text(size=12, face='bold'),
          strip.text = element_text(size=12, face='bold'),
          axis.line.x = element_line(colour="black", size=0.5),
          axis.line.y = element_line(colour="black", size=0.5)) +
    facet_wrap(~graph.metric)

dev.off()

