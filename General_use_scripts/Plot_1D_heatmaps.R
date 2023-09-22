#The input for this script is a tab-delimited text (.txt) matrix generated using deeptools plotProfile
#Typically, I will take the individual matrices for each replicate, average them in Excel bin-by-bin
#and save that as a new, 2-column txt file: Distance (-2kb, -1980 bp, etc.) and Coverage (bigwig score number)

setwd("~/Box/FACT Paper/Figures (old)//Fig 2 (CUT&RUN)/1D_Heatmaps/processed_matrices/")
#Change this location to match where your raw files are
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(scales)
library(cowplot)
#Load your package dependencies for plotting
list.files("~/Box/FACT Paper/Figures (old)//Fig 2 (CUT&RUN)/1D_Heatmaps/processed_matrices/")
#I just do this to make sure my wd is correct and my file names are correct- totally optional

K56ac_data=read.table("~/Desktop/K56ac.txt",header=TRUE,sep="\t")
#Read in your two-column .txt table
K56ac_data
#Take a look a the data, make sure everything looks right
K56ac_10E=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=X10E_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
#For each set of samples, you want to repeat this step, which will generate one 1D heatmap
#Since I have a lot of samples (7), there will be a lot of commands here
#R won't let you start a variable with a number, hence "X10E" even though it's just "10E" in the file
#Distance = bins/x-axis for plotting
#y = 1 (don't need to change this)
#fill = the value of color intensity at each position (so the y-axis on a metaplot, but the Z-axis on these 1D heatmaps)
#limits = the min and max values for heatmap color intensity
K56ac_11D=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=X11D_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_10E_No1=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=X10E_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_11D_No1=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=X11D_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_TIR1_V5=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=TIR1_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_TIR1_No1=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=TIR1_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_Input=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=Input_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_Spt16_ChIP=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=SPT16_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
K56ac_SSRP1_ChIP=ggplot(K56ac_data, aes(x = Distance, y = 1, fill=SSRP1_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)

K56ac<-(plot_grid(K56ac_TIR1_No1,K56ac_10E_No1,K56ac_11D_No1,K56ac_TIR1_V5,K56ac_10E,K56ac_11D,K56ac_Spt16_ChIP,K56ac_SSRP1_ChIP,K56ac_Input, ncol=1, align="hv"))
K56ac
#This plots all of the 1D heatmaps made above together in a row using the cowplot function plot_grid
#I like ncol=1 to keep everything aligned, but increase the number if you want more columns of heatmaps to save space
#align hv sets the heatmaps to anchor themselves together as the same size - highly recommend keeping this option
ggsave2(file="K56ac.png", plot=K56ac, path = "~/Desktop/", width = 25, height = 25, device='png', dpi=500)
#Saves the plot grid to a png in your designated path If you want tiff, change 'png' to 'tiff' after file and device
#500 dpi should be all you ever need, but change as necessary
#width and height are size in cm. If you go over 50 cm total, you have to add a flag to ignore size (it will tell you what)
