library(ggplot2)
library(dplyr)
library(hrbrthemes)

MillionSim3D <- read.table(file.choose(), sep =",",header=FALSE)
names(MillionSim3D) <- c("xRange","yRange","zRange","crossings")


test <- MillionSim3D %>% pivot_longer(cols = c("xRange","yRange","zRange"),names_to="type") %>% select(-crossings)

p <- test %>% ggplot(aes(x=value,fill=type))+geom_histogram(color="#e9ecef",alpha=0.6,position= 'identity',bins=50) +
  scale_fill_manual(values=c("#69b3a2","#f485aa","#404080")) + theme_ipsum() + labs(fill="") + labs(x="Maximum spread in Angstroms",y="Number of simulations")
p
ggsave("RangesIn3D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

q<-MillionSim3D %>% ggplot(aes(x=crossings))+geom_histogram(color="#e9ecef",fill="#805d0f",alpha=0.6,position= 'identity',bins=80) + theme_ipsum() + labs(fill="") + labs(x="Number of crossings",y="Number of simulations")
q


MillionSim2D <- read.table(file.choose(),sep=",",header=FALSE)
names(MillionSim2D) <- c("xRange","yRange","crossings","crossings_normalized")

test2 <- MillionSim2D %>% pivot_longer(cols = c("xRange","yRange"),names_to="type") %>% select(-crossings) %>% select(-crossings_normalized)

p2 <- test2 %>% ggplot(aes(x=value,fill=type))+geom_histogram(color="#e9ecef",alpha=0.6,position='identity',bins=50) +
  scale_fill_manual(values=c("#69b3a2","#f485aa")) + theme_ipsum() + labs(fill="") + labs(x="Maximum spread in Angstroms",y="Number of simulations")
p2
ggsave("RangesIn2D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)

q2 <- MillionSim2D %>% ggplot(aes(x=crossings))+geom_histogram(color="#e9ecef",fill="#805d0f",alpha=0.6,position='identity',bins=80) + theme_ipsum() + labs(fill="") + labs(x = "Number of crossings", y= "Number of simulations")
q2
ggsave("CrossingsIn2D.tiff",plot = last_plot(), width = 5, height =5, device = 'tiff', dpi=700)
