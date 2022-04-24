library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyr)

#-------------------------------------------------------------------------------------------
# The code below was used to generate figures and run analysis of data for simulations in 3D
#-------------------------------------------------------------------------------------------
MillionSim3D <- read.table(file.choose(), sep =",",header=FALSE)
names(MillionSim3D) <- c("xRange","yRange","zRange","crossings", "num4Cluster", "num5Cluster")


test0 <- MillionSim3D %>% pivot_longer(cols = c("xRange","yRange","zRange"),names_to="type") %>% select(-crossings) %>% select(-num4Cluster) %>% select(-num5Cluster)
test <- MillionSim3D %>% pivot_longer(cols = c("num4Cluster", "num5Cluster"),names_to="type") %>% select(-xRange) %>% select(-yRange) %>% select(-zRange) %>% select(-crossings)

p0 <- test0 %>% ggplot(aes(x=value,fill=type))+geom_histogram(color="#e9ecef",alpha=0.6,position= 'identity',bins=50) +
  scale_fill_manual(values=c("#69b3a2","#f485aa","#404080")) + theme_ipsum() + labs(fill="") + labs(x="Maximum spread in Angstroms",y="Number of simulations")
p0
ggsave("RangesIn3D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

q<-MillionSim3D %>% ggplot(aes(x=crossings))+geom_histogram(color="#e9ecef",fill="#805d0f",alpha=0.6,position= 'identity',bins=80, binwidth=7) + theme_ipsum() + labs(fill="") + labs(x="Number of crossings",y="Number of simulations")
q
ggsave("CrossingsIn3D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

p <- test %>% ggplot(aes(x=value, fill=type))+geom_histogram(color="#e9ecef",alpha=0.6, position='dodge',bins=50, binwidth=5) + 
  scale_fill_manual(values=c("#9a33ff", "#ff9a33"))+theme_ipsum() + labs(fill="") + labs(x="Number of clusters", y="Number of simulations")
p
ggsave("ClusterComparison3D.tiff", plot = last_plot(), width =5, height=5, device='tiff',dpi=700)

num4_3D <- subset(test, type== "num4Cluster", select = "value") %>% ggplot(aes(x=value))+geom_histogram(color="#e9ecef", fill="#9a33ff",alpha=0.6, position='identity', bins=50, binwidth=1)+theme_ipsum()+labs(fill="")+labs(x="Number of clusters", y="Number of simulations")
num4_3D
ggsave("4PointsInCluster3D.tiff", plot = last_plot(), width = 5, height =5, device='tiff',dpi=700)

data1 <- subset(test, type=="num4Cluster", select = "value")
summary(data1$value)
sd(data1$value)
1-sum(data1$value !=0)/length(data1$value)

data2 <- subset(test, type=="num5Cluster", select = "value")
summary(data2$value)
sd(data2$value)
1-sum(data2$value !=0)/length(data2$value)


#-------------------------------------------------------------------------------------------
# The code below was used to generate figures and run analysis of data for simulations in 2D
#-------------------------------------------------------------------------------------------
MillionSim2D <- read.table(file.choose(),sep=",",header=FALSE)
names(MillionSim2D) <- c("xRange","yRange","crossings","crossings_normalized", "num4Cluster", "num5Cluster")

test2 <- MillionSim2D %>% pivot_longer(cols = c("xRange","yRange"),names_to="type") %>% select(-crossings) %>% select(-crossings_normalized) %>% select(-num4Cluster) %>% select(-num5Cluster)
test1 <- MillionSim2D %>% pivot_longer(cols = c("num4Cluster", "num5Cluster"),names_to="type") %>% select(-xRange) %>% select(-yRange)%>% select(-crossings) %>% select(-crossings_normalized)

p2 <- test2 %>% ggplot(aes(x=value,fill=type))+geom_histogram(color="#e9ecef",alpha=0.6,position='identity',bins=50) +
  scale_fill_manual(values=c("#69b3a2","#f485aa")) + theme_ipsum() + labs(fill="") + labs(x="Maximum spread in Angstroms",y="Number of simulations")
p2
ggsave("RangesIn2D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

q2 <- MillionSim2D %>% ggplot(aes(x=crossings))+geom_histogram(color="#e9ecef",fill="#805d0f",alpha=0.6,position='identity',bins=80, binwidth=7) + theme_ipsum() + labs(fill="") + labs(x = "Number of crossings", y= "Number of simulations")
q2
ggsave("CrossingsIn2D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

p1 <- test1 %>% ggplot(aes(x=value, fill=type))+geom_histogram(color="#e9ecef",alpha=0.6, position='dodge',bins=50, binwidth=5) + 
  scale_fill_manual(values=c("#9a33ff", "#ff9a33"))+theme_ipsum() + labs(fill="") + labs(x="Number of clusters", y="Number of simulations")
p1
ggsave("ClusterComparison.tiff", plot = last_plot(), width =5, height=5, device='tiff',dpi=700)

num4 <- subset(test1, type== "num4Cluster", select = "value") %>% ggplot(aes(x=value))+geom_histogram(color="#e9ecef", fill="#9a33ff",alpha=0.6, position='identity', bins=50, binwidth=1)+theme_ipsum()+labs(fill="")+labs(x="Number of clusters", y="Number of simulations")
num4
ggsave("4PointsInCluster.tiff", plot = last_plot(), width = 5, height =5, device='tiff',dpi=700)

data <- subset(test1, type=="num4Cluster", select = "value")
summary(data$value)
sd(data$value)
1-sum(data$value !=0)/length(data$value)

data <- subset(test1, type=="num5Cluster", select = "value")
summary(data$value)
sd(data$value)
1-sum(data$value !=0)/length(data$value)