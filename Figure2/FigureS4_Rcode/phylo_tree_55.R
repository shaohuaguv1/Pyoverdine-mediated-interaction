library(circlize)
library(ggtreeExtra) # Set up overlay packages
library(ggplot2) # 
library(ggtree) # Draw a phylogenetic tree
library(treeio)
library(ggnewscale)


# ------Clear history variables-------------
rm(list = ls())

#---------------------Import data-----------------
phylo = read.csv('./data/phylo_out55.csv',header = T)
hc = read.tree('./data/test_gbk.tre')

#----------- phylogeny tree -----------
p <- ggtree(hc, layout = 'circular',branch.length='none') 
p

#----------- Siderophore production -----------
p1 <- p + new_scale_fill() + 
  geom_fruit(
    data=phylo,
    geom=geom_point,
    pwidth=0.1,
    mapping=aes(y=strains, x = 0.1, color=log10(abs(RFU))),
    #orientation="y",
    stat = 'identity'
  )+scale_color_gradientn(colours=c("grey","blue"))

p1

#----------- strain types -----------
p2 <- p1 + new_scale_fill() + 
  geom_fruit(
    data=phylo,
    geom=geom_bar,
    pwidth=0.1,
    mapping=aes(y=strains, x = 0.1, fill=sptype1),
    orientation="y",
    stat = 'identity'
  ) + 
  scale_fill_manual(
    
    values=c("#0E8A3B", "#E89112", "#B71C25","#FFFFFF"),
    
  ) 
p2



#----------- Synthetase groups -----------
p3 <- p2 + new_scale_fill() + 
  geom_fruit(
    data=phylo,
    geom=geom_bar,
    pwidth=0.1,
    mapping=aes(y=strains, x = 0.1, fill=syng1),
    orientation="y",
    stat = 'identity'
  )  + 
  scale_fill_manual(
    
    values=c( "#696969","#2e8b57", "#800000","#191970","#808000","#ff0000","#ff8c00","#ffd700","#ba55d3","#00fa9a","#00ffff","#0000ff","#f08080","#adff2f","#ff00ff","#1e90ff","#dda0dd","#ff1493","#87cefa","#ffe4b5","#FFFFFF"),
   
  ) 

p3

#----------- Receptor numbers / genome -----------
p4 <- p3 +  new_scale_fill() +
  geom_fruit(
    data=phylo,
    geom=geom_bar,
    pwidth=0.5,
    mapping=aes(y=strains, x=recnum, fill=recnum),
    orientation="y",
    stat = 'identity'
  )+scale_fill_gradientn(colours=c("#1e90ff","white","#ff00ff"))

p4
