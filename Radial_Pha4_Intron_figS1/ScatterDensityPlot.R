## Load the necessary libraries
library(ggplot2)
library(ggExtra)

## Set the working directory (if needed)
setwd("/Users/cokiny55/Desktop/Distribution/")
getwd()

# read csv file and save in dat list
list_of_files <- list.files(path = getwd(),
                            recursive = TRUE,
                            pattern = "\\.csv$",
                            full.names = FALSE)

#function for generating the stat_density plot
new.fig <- function(data,name_list) {
  myplot <- ggplot(data, aes(V1, V2))+
    geom_density_2d(aes(color = ..level..), size=1) +
    scale_color_viridis_c()+ # left plot with stat_density_2d
    geom_point(aes(x=V1, y=V2), size = 0.5,color='red',alpha = 1/3)+
    scale_x_continuous(expand = c(0, 0),limits = c(0, 0.7), breaks = seq(0, 1.0, by = 0.1))+
    scale_y_continuous(limits = c(2.5, 4.0), breaks = seq(2.5, 3.75, by = 0.25))+
    theme_classic()+
    theme(text = element_text(size = 18))+
    theme(panel.grid.major.x = element_line(color = "gray",
                                            linewidth = 0.5,
                                            linetype = 2))+
    theme(panel.grid.major.y = element_line(color = "gray",
                                            linewidth = 0.5,
                                            linetype = 2))+
    theme(legend.position = "none")+
    #ggtitle(sub(".csv", "", name_list))+ # get the string before csv
    ggtitle(qdapRegex::ex_between(name_list, "tion", ".csv"))+
    ylab("normalized intensity (log10)")+
    xlab("normalized distance to the periphery")+
    theme(plot.title = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(axis.text.x = element_text(angle = 90))
  
  plot<-ggMarginal(myplot, type="density", fill = "mediumpurple1", xparams = list(bins=400),yparams = list(bins=400))
}

plotlist <- list()
dat <- list()
for(i in 1:length(list_of_files)) {
  df <- read.csv(list_of_files[[i]], header=TRUE, sep = ",")
  #rename the plot
  colnames(df)[1]  <- "V1"
  colnames(df)[2]  <- "V2"
  df$V2 = log10(df$V2)
  dat[[i]] <- as.list(df)
  #generate the plot and save
  plotlist[[i]]<-new.fig(df,list_of_files[[i]])
  ggsave(sub(".csv", ".png", list_of_files[[i]]),plot = plotlist[[i]], device = "png")
}
