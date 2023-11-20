#The input for this script is a tab-delimited text (.txt) matrix generated using deeptools plotProfile
#Typically, I will take the individual matrices for each replicate, average them in Excel bin-by-bin
#and save that as a new, 2-column txt file: Distance in bp (-2000, -1980, etc.) and Coverage (bigwig score number)

library(ggplot2)
library(RColorBrewer)
library(gplots)
library(scales)
library(cowplot)

#Load your package dependencies for plotting
list.files("INPUT_DIR")
#I just do this to make sure my wd is correct and my file names are as expected- totally optional

10E=read.table("10E_Input_data.txt",header=TRUE,sep="\t")
#Read in your two-column .txt table
10E
#Take a look a the data, make sure everything looks right
10E=ggplot(data, aes(x = Distance, y = 1, fill=X10E_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
#For each set of samples, you want to repeat this step, which will generate one 1D heatmap
#Since I have a lot of samples (7), there will be a lot of commands here
#R won't let you start a variable with a number, hence "X10E" even though it's just "10E" in the file
#Distance = bins/x-axis for plotting
#y = 1 (don't need to change this)
#fill = the value of color intensity at each position (so the y-axis on a metaplot, but the Z-axis on these 1D heatmaps)
#limits = the min and max values for heatmap color intensity
11D=ggplot(data, aes(x = Distance, y = 1, fill=X11D_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
10E_No1=ggplot(data, aes(x = Distance, y = 1, fill=X10E_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
11D_No1=ggplot(data, aes(x = Distance, y = 1, fill=X11D_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
TIR1_V5=ggplot(data, aes(x = Distance, y = 1, fill=TIR1_V5_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
TIR1_No1=ggplot(data, aes(x = Distance, y = 1, fill=TIR1_No1_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
Input=ggplot(data, aes(x = Distance, y = 1, fill=Input_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
Spt16_ChIP=ggplot(data, aes(x = Distance, y = 1, fill=SPT16_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)
SSRP1_ChIP=ggplot(data, aes(x = Distance, y = 1, fill=SSRP1_ChIP_Coverage)) +
  geom_tile() +scale_fill_gradientn(limits = c(2,6), colours =  rich.colors(100),oob=squish)

combined_plot<-(plot_grid(TIR1_No1,10E_No1,11D_No1,TIR1_V5,10E,11D,Spt16_ChIP,SSRP1_ChIP,Input, ncol=1, align="hv"))
combined_plot
#This plots all of the 1D heatmaps made above together in a row using the cowplot function plot_grid
#I like ncol=1 to keep everything aligned, but increase the number if you want more columns of heatmaps to save space
#align hv sets the heatmaps to anchor themselves together as the same size - highly recommend keeping this option
ggsave2(file="cowplot.png", plot=combined_plot, path = "workdir", width = 25, height = 25, device='png', dpi=300)
#Saves the plot grid to a png in your designated path If you want tiff, change 'png' to 'tiff' after file and device
#300 dpi should be all you ever need, but change as necessary
#width and height are size in cm. If you go over 50 cm total, you have to add a flag to ignore size (it will tell you what)
