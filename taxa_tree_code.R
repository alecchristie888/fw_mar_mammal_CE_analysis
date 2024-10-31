#you may need to install some of these
#
devtools::install_github("GuangchuangYu/ggtree",force=TRUE)
library(ape)
library(cluster) 
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidytree)
library(ggtree)
library(treeio)
library(data.table)

#Create a phylogenetic tree of marine mammals using published trees and taxonomies
#Jefferson, T. A.; Leatherwood, S.; Webber, M. A. (1994). Marine Mammals of the World. Food and Agriculture Department of the United Nations. pp. 1-2. ISBN 978-92-5-103292-3. OCLC 30643250.
#Hrbek T, da Silva VMF, Dutra N, Gravena W, Martin AR, Farias IP (2014) A New Species of River Dolphin from Brazil or: How Little Do We Know Our Biodiversity. PLoS ONE 9(1): e83623. https://doi.org/10.1371/journal.pone.0083623

tt="(
(Dugongidae, Trichechidae),
((Phocidae, Odobenidae, Otariidae),
(((
(Kogiidae, Physeteridae),
(Ziphiidae,
(Platanistidae,
((Iniidae,  Pontoporiidae),
(Delphinidae, Monodontidae, Phocoenidae)))))),
(Balaenidae,
(Balaenopteridae, Eschrichtiidae)
)))
);"

tttree=read.tree(tt,file="")

#check basic tree
plot(tttree)

#get labels of family names
family_labs <- tttree$tip.label

#make a tibble (type of data frame or table) of data from the tree with information on nodes and parent groups
x = as_tibble(tttree)

#get CE dataset to get numbers of studies for each family
master1 = read.csv(choose.files()) #(fwmarinemammal_data.csv)
#subset to only mammals
marinemams = subset(master1,master1$class=="Mammalia")

#get number of studies for each family name from CE dataset
nstuds <- data.table(family=family_labs,studies=sapply(family_labs,function(i){ length(unique(marinemams$pageid[which(marinemams$family==i)]))}))
nstuds

#load in number of species per family 
#Species_per_family.csv
nspecies <- fread(choose.files())
names(nspecies) <- c("family","num")
nspecies[family=="Trichenchidae",family:="Trichechidae"]
nspecies[family=="Plantanistidae",family:="Platanistidae"]
#check which family names are found in the CE dataset
setdiff(nstuds$family,nspecies$family) #all matches
setdiff(nspecies$family,nstuds$family) #all matches

#make tibble of family names, studies, and number of species
d = tibble(merge(nstuds,nspecies,by="family",all=TRUE))
d$studspersp <- d$studies/d$num
d$label <- d$family 
d$studspersp[which(d$studspersp==0)] <- NA
#link the tree data with the data on studies and species
treecon1 = full_join(x, d, by='label')
treecon2 = as.treedata(treecon1)

#plot data as a circular tree (ignore warning message)
#just plotting data on number of studies
ggtree(treecon2, layout="fan", open.angle=45) + 
  geom_tiplab2(hjust=-0.25, size=5.95) + 
  geom_tippoint(aes(fill=studspersp),shape=21,size=15)+
  scale_fill_gradientn(name="Studies per species",breaks=c(0,2,4,6,8),limits=c(0,8),
                       colours=c("white",(brewer.pal(n=9, name="Blues"))[1:9]),
                       na.value="grey50")+
  theme(legend.position=c(0.85,0.45),axis.text = element_text(size=15),
        legend.text = element_text(size=17),
        legend.title = element_text(size=20))+xlim(0,10)

#we downloaded this and then prettified in inkscape with phylopic images
#ggsave("marinemammalstree.svg",device="svg",width=30,height=30,units="cm", dpi=600)


