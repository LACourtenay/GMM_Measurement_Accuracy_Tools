# Example of ggplot code used to create convex hulls and plot DBSCAN/MS results

# In the python code, extract the data you wish to plot via:
# np.savetxt("data_to_plot.txt",
#            np.stack((raw_data[:,0],raw_data[:,1],labels), axis = 1))
# This code has been included in the python only file

library(ggplot2)
library(tibble)
library(dplyr)
library(tidyverse)

directory<-"C:\\"
data_for_plot<-read.csv(paste(directory, "data_to_plot.txt", sep = "\\"), sep = ",") # Find data to be plot
data_for_plot$label<-as.factor(DBSCAN_all$label)

points<-data_for_plot %>% filter(label != "-1"); hull<-points %>% group_by(label) %>% slice(chull(x,y))
noise<-data_for_plot %>% filter(label == "-1")
centroids<-tibble(x = numeric(), y = numeric())
for (i in 1:length(levels(points$label))) {
  split_label<-split(points, points$label)
  split_label<-split_label[[i]]
  centroids<-add_row(centroids, x = mean(as.matrix(split_label[1])),
                     y = mean(as.matrix(split_label[2])))
  }; centroids<-centroids[2:nrow(centroids),] %>% as.data.frame
ggplot(data = points, aes(x = x, y = y, colour = label)) +
  geom_point(stat = "identity", size = 1.5, shape = 19) +
  geom_polygon(data = hull, aes(x = x, y = y), alpha = 0.2, size = 1) +
  geom_point(data = noise, aes(x = x, y = y), colour = "black",
             stat = "identity", size = 1, shape = 19) +
  geom_point(data = centroids, aes(x = x, y = y), colour = "black", stat = "identity",
             size = 15, shape = "*") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        plot.title = element_text(face = "bold", size = 20),
        plot.subtitle = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none") +
  coord_flip() +
  scale_y_reverse()
