# Libraries ----------------------

library(geomorph)
library(shapes)
library(tibble)
library(dplyr)
library(tidyverse)
library(spatstat)
library(MVN)
library(asbio)
library(car)
library(RVAideMemoire)
library(Evomorph)
library(e1071)
library(cluster)
library(asbio) # If errors occur with this package consult ReadMe.

# obtaining data ---------------

directory<-"C:\\"

a<-read.morphologika(paste(directory, "All_Landmarks.txt", sep = "\\"))
A1<-read.morphologika(paste(directory, "A1.txt", sep = "\\"))
A2<-read.morphologika(paste(directory, "A2.txt", sep = "\\"))
A3<-read.morphologika(paste(directory, "A3.txt", sep = "\\"))
Wolf<-read.morphologika(paste(directory, "Wolf.txt", sep = "\\"))
Dog<-read.morphologika(paste(directory, "Dog.txt", sep = "\\"))

GPA<-gpagen(a$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
            max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
A1_GPA<-gpagen(A1$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
A2_GPA<-gpagen(A2$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                 max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
A3_GPA<-gpagen(A3$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Wolf_GPA<-gpagen(Wolf$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                 max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Dog_GPA<-gpagen(Dog$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                  max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)

raw_mshape<-data.frame(mshape(a$coords))

absolute_mshape<-data.frame(mshape(GPA$coords));names(absolute_mshape)[1]<-"X";names(absolute_mshape)[2]<-"Y";names(absolute_mshape)[3]<-"Z"
A1_mshape<-data.frame(mshape(A1_GPA$coords));names(A1_mshape)[1]<-"X";names(A1_mshape)[2]<-"Y";names(A1_mshape)[3]<-"Z"
A2_mshape<-data.frame(mshape(A2_GPA$coords));names(A2_mshape)[1]<-"X";names(A2_mshape)[2]<-"Y";names(A2_mshape)[3]<-"Z"
A3_mshape<-data.frame(mshape(A3_GPA$coords));names(A3_mshape)[1]<-"X";names(A3_mshape)[2]<-"Y";names(A3_mshape)[3]<-"Z"
Wolf_mshape<-data.frame(mshape(Wolf_GPA$coords));names(Wolf_mshape)[1]<-"X";names(Wolf_mshape)[2]<-"Y";names(Wolf_mshape)[3]<-"Z"
Dog_mshape<-data.frame(mshape(Dog_GPA$coords));names(Dog_mshape)[1]<-"X";names(Dog_mshape)[2]<-"Y";names(Dog_mshape)[3]<-"Z"

# Automate procrustes data collection

data<-tibble(label = character(), Landmark = factor(), Sample = factor(), Animal = factor(),
             x = numeric(), y = numeric(), z = numeric())
for (i in 1:length(a$labels)){
  sep<-GPA$coords[,,i]
  for (i2 in 1:nrow(sep)){
    LM<-sep[i2,]
    data<-add_row(data, label = rownames(a$labels)[i],
                  Landmark = paste("LM", i2, sep = ""),
                  Sample = a$labels[i],
                  Animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                  x = LM[1], y = LM[2], z = LM[3])
  }
}; data<-as.tibble(data) %>%
  mutate(label = if_else(grepl("Wolf", label) == TRUE,
                               substr(label, 1, 11), substr(label, 1, 12))); data$label<-factor(data$label); rm(sep, LM)

Wolf_data<-tibble(label = character(), Landmark = factor(), Sample = factor(), Animal = factor(),
                  x = numeric(), y = numeric(), z = numeric())
for (i in 1:length(Wolf$labels)){
  sep<-Wolf_GPA$coords[,,i]
  for (i2 in 1:nrow(sep)){
    LM<-sep[i2,]
    Wolf_data<-add_row(Wolf_data, label = rownames(Wolf$labels)[i],
                       Landmark = paste("LM", i2, sep = ""),
                       Sample = Wolf$labels[i],
                       Animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                       x = LM[1], y = LM[2], z = LM[3])
  }
}; Wolf_data<-as.tibble(Wolf_data) %>%
  mutate(label = if_else(grepl("Wolf", label) == TRUE,
                         substr(label, 1, 11), substr(label, 1, 12))); Wolf_data$label<-factor(Wolf_data$label); rm(sep, LM)

Dog_data<-tibble(label = character(), Landmark = factor(), Sample = factor(), Animal = factor(),
                 x = numeric(), y = numeric(), z = numeric())
for (i in 1:length(Dog$labels)){
  sep<-Dog_GPA$coords[,,i]
  for (i2 in 1:nrow(sep)){
    LM<-sep[i2,]
    Dog_data<-add_row(Dog_data, label = rownames(Dog$labels)[i],
                      Landmark = paste("LM", i2, sep = ""),
                      Sample = Dog$labels[i],
                      Animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                      x = LM[1], y = LM[2], z = LM[3])
  }
}; Dog_data<-as.tibble(Dog_data) %>%
  mutate(label = if_else(grepl("Wolf", label) == TRUE,
                         substr(label, 1, 11), substr(label, 1, 12))); Dog_data$label<-factor(Dog_data$label); rm(sep, LM)

# Distance Calculations

distances<-tibble(label = character(), Landmark = factor(), Animal = factor(),
                  A1_x = numeric(), A1_y = numeric(), A1_z = numeric(),
                  A3_x = numeric(), A3_y = numeric(), A3_z = numeric(),
                  A2_x = numeric(), A2_y = numeric(), A2_z = numeric())
split_label<-split(data, data$label)
for (i in 1:length(split_label)){
  seperate_data<-split_label[[i]]
  split_landmark<-split(seperate_data, seperate_data$Landmark)
  for (i2 in 1:length(split_landmark)){
    distances<-add_row(distances, label = split_landmark[[i2]]$label[1],
                       Landmark = split_landmark[[i2]]$Landmark[1],
                       Animal = split_landmark[[i2]]$Animal[1],
                       A1_x = split_landmark[[i2]]$x[1], A1_y = split_landmark[[i2]]$y[1], A1_z = split_landmark[[i2]]$z[1],
                       A2_x = split_landmark[[i2]]$x[2], A2_y = split_landmark[[i2]]$y[2], A2_z = split_landmark[[i2]]$z[2],
                       A3_x = split_landmark[[i2]]$x[3], A3_y = split_landmark[[i2]]$y[3], A3_z = split_landmark[[i2]]$z[3])
  }
}; rm(seperate_data, split_landmark, split_label)
distances<-as.tibble(distances) %>%
  mutate(
    A1_A3 = sqrt(((A1_x-A3_x)^2)+((A1_y-A3_y)^2)+((A1_y-A3_y)^2)),
    A1_A2 = sqrt(((A1_x-A2_x)^2)+((A1_y-A2_y)^2)+((A1_y-A2_y)^2)),
    A3_A2 = sqrt(((A3_x-A2_x)^2)+((A3_y-A2_y)^2)+((A3_y-A2_y)^2))
  )

#

# Obtaining interanalyst data - Seperating by Animals ----------------

# Wolf

Wolf_A1<-read.morphologika(paste(directory, "\\Wolf_A1.txt", sep = ""))
Wolf_A2<-read.morphologika(paste(directory, "\\Wolf_A2.txt", sep = ""))
Wolf_A3<-read.morphologika(paste(directory, "\\Wolf_A3.txt", sep = ""))
Wolf_A1_GPA<-gpagen(Wolf_A1$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                     max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Wolf_A2_GPA<-gpagen(Wolf_A2$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Wolf_A3_GPA<-gpagen(Wolf_A3$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                     max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Wolf_A1_mshape<-data.frame(mshape(Wolf_A1_GPA$coords))
Wolf_A2_mshape<-data.frame(mshape(Wolf_A2_GPA$coords))
Wolf_A3_mshape<-data.frame(mshape(Wolf_A3_GPA$coords))

# Dog

Dog_A1<-read.morphologika(paste(directory, "\\Dog_A1.txt", sep = ""))
Dog_A2<-read.morphologika(paste(directory, "\\Dog_A2.txt", sep = ""))
Dog_A3<-read.morphologika(paste(directory, "\\Dog_A3.txt", sep = ""))
Dog_A1_GPA<-gpagen(Dog_A1$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Dog_A2_GPA<-gpagen(Dog_A2$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                       max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Dog_A3_GPA<-gpagen(Dog_A3$coords, curves = NULL, surfaces = NULL, PrinAxes = TRUE,
                      max.iter = NULL, ProcD = TRUE, Proj = TRUE, print.progress = TRUE)
Dog_A1_mshape<-data.frame(mshape(Dog_A1_GPA$coords))
Dog_A2_mshape<-data.frame(mshape(Dog_A2_GPA$coords))
Dog_A3_mshape<-data.frame(mshape(Dog_A3_GPA$coords))

#

# Calculate distributions --------------------

ALM<-split(data, data$Landmark)
for (i in 1:length(ALM)){
  Landmark_n<-paste("LM", i, sep = "")
  LM<-ALM[[Landmark_n]]
  xrange<-c(min(LM$y),max(LM$y))
  yrange<-c(min(LM$z),max(LM$z))
  x<-ppp(LM$x, LM$z, xrange, yrange) # Creating Spatial Window
  print(quadrat.test(x,5,5)$p.value) # p<0.05 - distribuci�n inhomogeneo
}; rm(LM, Landmark_n, xrange, yrange, x)
for (i in 1:length(ALM)){
  Landmark_n<-paste("LM", i, sep = "")
  LM<-ALM[[Landmark_n]]
  xrange<-c(min(LM$y),max(LM$y))
  yrange<-c(min(LM$z),max(LM$z))
  x<-ppp(LM$y, LM$z, xrange, yrange) # Creating Spatial Window
  print(quadrat.test(x,5,5)$statistic[[1]])
}; rm(LM, Landmark_n, xrange, yrange, x)

# For each analyst seperately

ALM<-split(data, data$Animal); ALM<-ALM[[2]]; ALM<-split(ALM, ALM$Landmark)
for (i in 1:length(ALM)){
  Landmark_n <- paste("LM", i, sep = "")
  LM<-ALM[[Landmark_n]]; LM<-LM[-c(2)]
  xrange<-c(min(LM$x),max(LM$x))
  yrange<-c(min(LM$z),max(LM$z))
  x<-ppp(LM$x, LM$z, xrange, yrange) # Creating Spatial Window
  print(quadrat.test(x,5,5)$p.value) # p<0.05 - distribuci�n inhomogeneo
}; rm(LM, Landmark_n, xrange, yrange, x)
for (i in 1:length(ALM)){
  Landmark_n <- paste("LM", i, sep = "")
  LM<-ALM[[Landmark_n]]; LM<-LM[-c(2)]
  xrange<-c(min(LM$x),max(LM$x))
  yrange<-c(min(LM$z),max(LM$z))
  x<-ppp(LM$x, LM$z, xrange, yrange) # Creating Spatial Window
  print(quadrat.test(x,5,5)$statistic[[1]])
}; rm(LM, Landmark_n, xrange, yrange, x)

#

# Absolute Error Calculations and Distributions -------------------
# A1 -----------------------------------

# Error distance calculations

A1_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A1") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((A1_mshape$X[1]-x)^2)+((A1_mshape$Y[1]-y)^2)+((A1_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((A1_mshape$X[2]-x)^2)+((A1_mshape$Y[2]-y)^2)+((A1_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((A1_mshape$X[3]-x)^2)+((A1_mshape$Y[3]-y)^2)+((A1_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((A1_mshape$X[4]-x)^2)+((A1_mshape$Y[4]-y)^2)+((A1_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((A1_mshape$X[5]-x)^2)+((A1_mshape$Y[6]-y)^2)+((A1_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((A1_mshape$X[6]-x)^2)+((A1_mshape$Y[6]-y)^2)+((A1_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((A1_mshape$X[7]-x)^2)+((A1_mshape$Y[7]-y)^2)+((A1_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((A1_mshape$X[8]-x)^2)+((A1_mshape$Y[8]-y)^2)+((A1_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((A1_mshape$X[2]-x)^2)+((A1_mshape$Y[2]-y)^2)+((A1_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((A1_mshape$X[10]-x)^2)+((A1_mshape$Y[10]-y)^2)+((A1_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((A1_mshape$X[11]-x)^2)+((A1_mshape$Y[11]-y)^2)+((A1_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((A1_mshape$X[12]-x)^2)+((A1_mshape$Y[12]-y)^2)+((A1_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((A1_mshape$X[13]-x)^2)+((A1_mshape$Y[13]-y)^2)+((A1_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((A1_mshape$X[14]-x)^2)+((A1_mshape$Y[14]-y)^2)+((A1_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((A1_mshape$X[15]-x)^2)+((A1_mshape$Y[15]-y)^2)+((A1_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((A1_mshape$X[16]-x)^2)+((A1_mshape$Y[16]-y)^2)+((A1_mshape$Z[16]-z)^2)),
                                  sqrt(((A1_mshape$X[17]-x)^2)+((A1_mshape$Y[17]-y)^2)+((A1_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  ))

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(A1_error_distances$distance,
     main = "Histogram of Distribution of Distances from Centroid for A1",
     sub = if (
       shapiro.test(A1_error_distances$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(A1_error_distances$distance)$`p.value` > 0.05 && shapiro.test(A1_error_distances$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(A1_error_distances$distance)$`p.value` > 0.001 && shapiro.test(A1_error_distances$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(A1_error_distances$distance),
                    sd(A1_error_distances$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(A1_error_distances$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(A1_error_distances$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of Distances from Centroid for A1",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(A1_error_distances$distance)

# NMAD and BWMV Calculations

A1_landmark_errors<-tibble(Landmark = factor(), Median = numeric(), NMAD = numeric(), BWMV = numeric())
for (i in 1:length(levels(A1_error_distances$Landmark))){
  A1_landmark_error<-as.tibble(A1_error_distances) %>% filter(Landmark == paste("LM", i, sep = ""))
  A1_landmark_errors<-add_row(A1_landmark_errors,
                               Landmark = paste("LM", i, sep = ""),
                               Median = median(A1_landmark_error$distance),
                               NMAD = mad(A1_landmark_error$distance, constant = 1.4826),
                               BWMV = sqrt(r.bw(A1_landmark_error$distance)$"S.xx"[1]))
}; rm(A1_landmark_error); A1_landmark_errors

#

# A2 -----------------------------------

# Error distance calculations

A2_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A2") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((A2_mshape$X[1]-x)^2)+((A2_mshape$Y[1]-y)^2)+((A2_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((A2_mshape$X[2]-x)^2)+((A2_mshape$Y[2]-y)^2)+((A2_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((A2_mshape$X[3]-x)^2)+((A2_mshape$Y[3]-y)^2)+((A2_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((A2_mshape$X[4]-x)^2)+((A2_mshape$Y[4]-y)^2)+((A2_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((A2_mshape$X[5]-x)^2)+((A2_mshape$Y[6]-y)^2)+((A2_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((A2_mshape$X[6]-x)^2)+((A2_mshape$Y[6]-y)^2)+((A2_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((A2_mshape$X[7]-x)^2)+((A2_mshape$Y[7]-y)^2)+((A2_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((A2_mshape$X[8]-x)^2)+((A2_mshape$Y[8]-y)^2)+((A2_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((A2_mshape$X[2]-x)^2)+((A2_mshape$Y[2]-y)^2)+((A2_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((A2_mshape$X[10]-x)^2)+((A2_mshape$Y[10]-y)^2)+((A2_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((A2_mshape$X[11]-x)^2)+((A2_mshape$Y[11]-y)^2)+((A2_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((A2_mshape$X[12]-x)^2)+((A2_mshape$Y[12]-y)^2)+((A2_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((A2_mshape$X[13]-x)^2)+((A2_mshape$Y[13]-y)^2)+((A2_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((A2_mshape$X[14]-x)^2)+((A2_mshape$Y[14]-y)^2)+((A2_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((A2_mshape$X[15]-x)^2)+((A2_mshape$Y[15]-y)^2)+((A2_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((A2_mshape$X[16]-x)^2)+((A2_mshape$Y[16]-y)^2)+((A2_mshape$Z[16]-z)^2)),
                                  sqrt(((A2_mshape$X[17]-x)^2)+((A2_mshape$Y[17]-y)^2)+((A2_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  ))

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(A2_error_distances$distance,
     main = "Histogram of Distribution of Distances from Centroid for A2",
     sub = if (
       shapiro.test(A2_error_distances$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(A2_error_distances$distance)$`p.value` > 0.05 && shapiro.test(A2_error_distances$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(A2_error_distances$distance)$`p.value` > 0.001 && shapiro.test(A2_error_distances$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(A2_error_distances$distance),
                    sd(A2_error_distances$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(A2_error_distances$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(A2_error_distances$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of Distances from Centroid for A2",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(A2_error_distances$distance)

# NMAD and BWMV Calculations

A2_landmark_errors<-tibble(Landmark = factor(), Median = numeric(), NMAD = numeric(), BWMV = numeric())
for (i in 1:length(levels(A2_error_distances$Landmark))){
  A2_landmark_error<-as.tibble(A2_error_distances) %>% filter(Landmark == paste("LM", i, sep = ""))
  A2_landmark_errors<-add_row(A2_landmark_errors,
                                Landmark = paste("LM", i, sep = ""),
                                Median = median(A2_landmark_error$distance),
                                NMAD = mad(A2_landmark_error$distance, constant = 1.4826),
                                BWMV = sqrt(r.bw(A2_landmark_error$distance)$"S.xx"[1]))
}; rm(A2_landmark_error); A2_landmark_errors

#

# A3 -----------------------------------

# Error distance calculations

A3_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A3") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((A3_mshape$X[1]-x)^2)+((A3_mshape$Y[1]-y)^2)+((A3_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((A3_mshape$X[2]-x)^2)+((A3_mshape$Y[2]-y)^2)+((A3_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((A3_mshape$X[3]-x)^2)+((A3_mshape$Y[3]-y)^2)+((A3_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((A3_mshape$X[4]-x)^2)+((A3_mshape$Y[4]-y)^2)+((A3_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((A3_mshape$X[5]-x)^2)+((A3_mshape$Y[6]-y)^2)+((A3_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((A3_mshape$X[6]-x)^2)+((A3_mshape$Y[6]-y)^2)+((A3_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((A3_mshape$X[7]-x)^2)+((A3_mshape$Y[7]-y)^2)+((A3_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((A3_mshape$X[8]-x)^2)+((A3_mshape$Y[8]-y)^2)+((A3_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((A3_mshape$X[2]-x)^2)+((A3_mshape$Y[2]-y)^2)+((A3_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((A3_mshape$X[10]-x)^2)+((A3_mshape$Y[10]-y)^2)+((A3_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((A3_mshape$X[11]-x)^2)+((A3_mshape$Y[11]-y)^2)+((A3_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((A3_mshape$X[12]-x)^2)+((A3_mshape$Y[12]-y)^2)+((A3_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((A3_mshape$X[13]-x)^2)+((A3_mshape$Y[13]-y)^2)+((A3_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((A3_mshape$X[14]-x)^2)+((A3_mshape$Y[14]-y)^2)+((A3_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((A3_mshape$X[15]-x)^2)+((A3_mshape$Y[15]-y)^2)+((A3_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((A3_mshape$X[16]-x)^2)+((A3_mshape$Y[16]-y)^2)+((A3_mshape$Z[16]-z)^2)),
                                  sqrt(((A3_mshape$X[17]-x)^2)+((A3_mshape$Y[17]-y)^2)+((A3_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  ))

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(A3_error_distances$distance,
     main = "Histogram of Distribution of Distances from Centroid for A3",
     sub = if (
       shapiro.test(A3_error_distances$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(A3_error_distances$distance)$`p.value` > 0.05 && shapiro.test(A3_error_distances$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(A3_error_distances$distance)$`p.value` > 0.001 && shapiro.test(A3_error_distances$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(A3_error_distances$distance),
                    sd(A3_error_distances$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(A3_error_distances$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(A3_error_distances$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of Distances from Centroid for A3",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(A3_error_distances$distance)

# NMAD and BWMV Calculations

A3_landmark_errors<-tibble(Landmark = factor(), Median = numeric(), NMAD = numeric(), BWMV = numeric())
for (i in 1:length(levels(A3_error_distances$Landmark))){
  A3_landmark_error<-as.tibble(A3_error_distances) %>% filter(Landmark == paste("LM", i, sep = ""))
  A3_landmark_errors<-add_row(A3_landmark_errors,
                               Landmark = paste("LM", i, sep = ""),
                               Median = median(A3_landmark_error$distance),
                               NMAD = mad(A3_landmark_error$distance, constant = 1.4826),
                               BWMV = sqrt(r.bw(A3_landmark_error$distance)$"S.xx"[1]))
}; rm(A3_landmark_error); A3_landmark_errors

#

# InterAnalyst Variability ----------------------------

A3_inter_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A3") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((absolute_mshape$X[1]-x)^2)+((absolute_mshape$Y[1]-y)^2)+((absolute_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((absolute_mshape$X[3]-x)^2)+((absolute_mshape$Y[3]-y)^2)+((absolute_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((absolute_mshape$X[4]-x)^2)+((absolute_mshape$Y[4]-y)^2)+((absolute_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((absolute_mshape$X[5]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((absolute_mshape$X[6]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((absolute_mshape$X[7]-x)^2)+((absolute_mshape$Y[7]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((absolute_mshape$X[8]-x)^2)+((absolute_mshape$Y[8]-y)^2)+((absolute_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((absolute_mshape$X[10]-x)^2)+((absolute_mshape$Y[10]-y)^2)+((absolute_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((absolute_mshape$X[11]-x)^2)+((absolute_mshape$Y[11]-y)^2)+((absolute_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((absolute_mshape$X[12]-x)^2)+((absolute_mshape$Y[12]-y)^2)+((absolute_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((absolute_mshape$X[13]-x)^2)+((absolute_mshape$Y[13]-y)^2)+((absolute_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((absolute_mshape$X[14]-x)^2)+((absolute_mshape$Y[14]-y)^2)+((absolute_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((absolute_mshape$X[15]-x)^2)+((absolute_mshape$Y[15]-y)^2)+((absolute_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((absolute_mshape$X[16]-x)^2)+((absolute_mshape$Y[16]-y)^2)+((absolute_mshape$Z[16]-z)^2)),
                                  sqrt(((absolute_mshape$X[17]-x)^2)+((absolute_mshape$Y[17]-y)^2)+((absolute_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A1_inter_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A1") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((absolute_mshape$X[1]-x)^2)+((absolute_mshape$Y[1]-y)^2)+((absolute_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((absolute_mshape$X[3]-x)^2)+((absolute_mshape$Y[3]-y)^2)+((absolute_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((absolute_mshape$X[4]-x)^2)+((absolute_mshape$Y[4]-y)^2)+((absolute_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((absolute_mshape$X[5]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((absolute_mshape$X[6]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((absolute_mshape$X[7]-x)^2)+((absolute_mshape$Y[7]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((absolute_mshape$X[8]-x)^2)+((absolute_mshape$Y[8]-y)^2)+((absolute_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((absolute_mshape$X[10]-x)^2)+((absolute_mshape$Y[10]-y)^2)+((absolute_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((absolute_mshape$X[11]-x)^2)+((absolute_mshape$Y[11]-y)^2)+((absolute_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((absolute_mshape$X[12]-x)^2)+((absolute_mshape$Y[12]-y)^2)+((absolute_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((absolute_mshape$X[13]-x)^2)+((absolute_mshape$Y[13]-y)^2)+((absolute_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((absolute_mshape$X[14]-x)^2)+((absolute_mshape$Y[14]-y)^2)+((absolute_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((absolute_mshape$X[15]-x)^2)+((absolute_mshape$Y[15]-y)^2)+((absolute_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((absolute_mshape$X[16]-x)^2)+((absolute_mshape$Y[16]-y)^2)+((absolute_mshape$Z[16]-z)^2)),
                                  sqrt(((absolute_mshape$X[17]-x)^2)+((absolute_mshape$Y[17]-y)^2)+((absolute_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A2_inter_error_distances<-tibble(Landmark = data$Landmark, Sample = data$Sample, x = data$x, y = data$y, z = data$z) %>%
  filter(Sample == "A2") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((absolute_mshape$X[1]-x)^2)+((absolute_mshape$Y[1]-y)^2)+((absolute_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((absolute_mshape$X[3]-x)^2)+((absolute_mshape$Y[3]-y)^2)+((absolute_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((absolute_mshape$X[4]-x)^2)+((absolute_mshape$Y[4]-y)^2)+((absolute_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((absolute_mshape$X[5]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((absolute_mshape$X[6]-x)^2)+((absolute_mshape$Y[6]-y)^2)+((absolute_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((absolute_mshape$X[7]-x)^2)+((absolute_mshape$Y[7]-y)^2)+((absolute_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((absolute_mshape$X[8]-x)^2)+((absolute_mshape$Y[8]-y)^2)+((absolute_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((absolute_mshape$X[2]-x)^2)+((absolute_mshape$Y[2]-y)^2)+((absolute_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((absolute_mshape$X[10]-x)^2)+((absolute_mshape$Y[10]-y)^2)+((absolute_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((absolute_mshape$X[11]-x)^2)+((absolute_mshape$Y[11]-y)^2)+((absolute_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((absolute_mshape$X[12]-x)^2)+((absolute_mshape$Y[12]-y)^2)+((absolute_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((absolute_mshape$X[13]-x)^2)+((absolute_mshape$Y[13]-y)^2)+((absolute_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((absolute_mshape$X[14]-x)^2)+((absolute_mshape$Y[14]-y)^2)+((absolute_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((absolute_mshape$X[15]-x)^2)+((absolute_mshape$Y[15]-y)^2)+((absolute_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((absolute_mshape$X[16]-x)^2)+((absolute_mshape$Y[16]-y)^2)+((absolute_mshape$Z[16]-z)^2)),
                                  sqrt(((absolute_mshape$X[17]-x)^2)+((absolute_mshape$Y[17]-y)^2)+((absolute_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); interanalyst_errors<-full_join(A2_inter_error_distances,
                                     A3_inter_error_distances); interanalyst_errors<-full_join(interanalyst_errors,
                                                                                                A1_inter_error_distances)
rm(A2_inter_error_distances, A3_inter_error_distances, A1_inter_error_distances)

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(interanalyst_errors$distance,
     main = "Histogram of Distribution of interanalyst Distances from Centroid",
     sub = if (
       shapiro.test(interanalyst_errors$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(interanalyst_errors$distance)$`p.value` > 0.05 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(interanalyst_errors$distance)$`p.value` > 0.001 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(interanalyst_errors$distance),
                    sd(interanalyst_errors$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(interanalyst_errors$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(interanalyst_errors$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of interanalyst Distances from Centroid",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(interanalyst_errors$distance)

# Shapiros
shapiro_values<-tibble(Sample = factor(), w.value = numeric(), p.value = numeric())
for (i in c("A1", "A2", "A3")){
    shapiro_values<-shapiro_values %>% add_row(Sample = i,
                                               p.value = shapiro.test(
                                                 as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                                 )$"p.value",
                                               w.value = shapiro.test(
                                                 as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                                 )$"statistic"[["W"]]
                                               )
    }; shapiro_values

# NMAD and BWMV Calculations

interanalyst_landmark_errors<-tibble(Landmark = factor(),
                                     overall_median = numeric(),
                                     overall_NMAD = numeric(), overall_BWMV = numeric(),
                                     A3_median = numeric(),
                                     A3_NMAD = numeric(), A3_BWMV = numeric(),
                                     A1_median = numeric(),
                                     A1_NMAD = numeric(), A1_BWMV = numeric(),
                                     A2_median = numeric(),
                                     A2_NMAD = numeric(), A2_BWMV = numeric())
for (i in 1:length(levels(interanalyst_errors$Landmark))){
  interanalyst_landmark_error<-as.tibble(interanalyst_errors) %>% filter(Landmark == paste("LM", i, sep = ""))
  interanalyst_landmark_errors<-add_row(interanalyst_landmark_errors,
                                        Landmark = paste("LM", i, sep = ""),
                                        overall_median = median(interanalyst_landmark_error$distance),
                                        overall_NMAD = mad(interanalyst_landmark_error$distance, constant = 1.4826),
                                        overall_BWMV = sqrt(r.bw(interanalyst_landmark_error$distance)$"S.xx"[1]),
                                        A3_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance))),
                                        A3_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A3_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A1_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance))),
                                        A1_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A1_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A2_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance))),
                                        A2_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A2_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)))$"S.xx"[1]))
  }; rm(interanalyst_landmark_error); interanalyst_landmark_errors

# Comparisons According to Animals -----------------------------

# Wolf

A3_inter_error_distances<-tibble(Landmark = Wolf_data$Landmark, Animal = Wolf_data$Animal,
                                  Sample = Wolf_data$Sample,
                                  x = Wolf_data$x, y = Wolf_data$y,z = Wolf_data$z) %>%
  filter(Animal == "Wolf") %>%
  filter(Sample == "A3") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Wolf_mshape$X[1]-x)^2)+((Wolf_mshape$Y[1]-y)^2)+((Wolf_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Wolf_mshape$X[3]-x)^2)+((Wolf_mshape$Y[3]-y)^2)+((Wolf_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Wolf_mshape$X[4]-x)^2)+((Wolf_mshape$Y[4]-y)^2)+((Wolf_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Wolf_mshape$X[5]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Wolf_mshape$X[6]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Wolf_mshape$X[7]-x)^2)+((Wolf_mshape$Y[7]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Wolf_mshape$X[8]-x)^2)+((Wolf_mshape$Y[8]-y)^2)+((Wolf_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Wolf_mshape$X[10]-x)^2)+((Wolf_mshape$Y[10]-y)^2)+((Wolf_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Wolf_mshape$X[11]-x)^2)+((Wolf_mshape$Y[11]-y)^2)+((Wolf_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Wolf_mshape$X[12]-x)^2)+((Wolf_mshape$Y[12]-y)^2)+((Wolf_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Wolf_mshape$X[13]-x)^2)+((Wolf_mshape$Y[13]-y)^2)+((Wolf_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Wolf_mshape$X[14]-x)^2)+((Wolf_mshape$Y[14]-y)^2)+((Wolf_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Wolf_mshape$X[15]-x)^2)+((Wolf_mshape$Y[15]-y)^2)+((Wolf_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Wolf_mshape$X[16]-x)^2)+((Wolf_mshape$Y[16]-y)^2)+((Wolf_mshape$Z[16]-z)^2)),
                                  sqrt(((Wolf_mshape$X[17]-x)^2)+((Wolf_mshape$Y[17]-y)^2)+((Wolf_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A1_inter_error_distances<-tibble(Landmark = Wolf_data$Landmark, Animal = Wolf_data$Animal, Sample = Wolf_data$Sample, x = Wolf_data$x, y = Wolf_data$y, z = Wolf_data$z) %>%
  filter(Animal == "Wolf") %>%
  filter(Sample == "A1") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Wolf_mshape$X[1]-x)^2)+((Wolf_mshape$Y[1]-y)^2)+((Wolf_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Wolf_mshape$X[3]-x)^2)+((Wolf_mshape$Y[3]-y)^2)+((Wolf_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Wolf_mshape$X[4]-x)^2)+((Wolf_mshape$Y[4]-y)^2)+((Wolf_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Wolf_mshape$X[5]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Wolf_mshape$X[6]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Wolf_mshape$X[7]-x)^2)+((Wolf_mshape$Y[7]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Wolf_mshape$X[8]-x)^2)+((Wolf_mshape$Y[8]-y)^2)+((Wolf_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Wolf_mshape$X[10]-x)^2)+((Wolf_mshape$Y[10]-y)^2)+((Wolf_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Wolf_mshape$X[11]-x)^2)+((Wolf_mshape$Y[11]-y)^2)+((Wolf_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Wolf_mshape$X[12]-x)^2)+((Wolf_mshape$Y[12]-y)^2)+((Wolf_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Wolf_mshape$X[13]-x)^2)+((Wolf_mshape$Y[13]-y)^2)+((Wolf_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Wolf_mshape$X[14]-x)^2)+((Wolf_mshape$Y[14]-y)^2)+((Wolf_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Wolf_mshape$X[15]-x)^2)+((Wolf_mshape$Y[15]-y)^2)+((Wolf_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Wolf_mshape$X[16]-x)^2)+((Wolf_mshape$Y[16]-y)^2)+((Wolf_mshape$Z[16]-z)^2)),
                                  sqrt(((Wolf_mshape$X[17]-x)^2)+((Wolf_mshape$Y[17]-y)^2)+((Wolf_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A2_inter_error_distances<-tibble(Landmark = Wolf_data$Landmark, Animal = Wolf_data$Animal, Sample = Wolf_data$Sample, x = Wolf_data$x, y = Wolf_data$y, z = Wolf_data$z) %>%
  filter(Animal == "Wolf") %>%
  filter(Sample == "A2") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Wolf_mshape$X[1]-x)^2)+((Wolf_mshape$Y[1]-y)^2)+((Wolf_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Wolf_mshape$X[3]-x)^2)+((Wolf_mshape$Y[3]-y)^2)+((Wolf_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Wolf_mshape$X[4]-x)^2)+((Wolf_mshape$Y[4]-y)^2)+((Wolf_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Wolf_mshape$X[5]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Wolf_mshape$X[6]-x)^2)+((Wolf_mshape$Y[6]-y)^2)+((Wolf_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Wolf_mshape$X[7]-x)^2)+((Wolf_mshape$Y[7]-y)^2)+((Wolf_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Wolf_mshape$X[8]-x)^2)+((Wolf_mshape$Y[8]-y)^2)+((Wolf_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Wolf_mshape$X[2]-x)^2)+((Wolf_mshape$Y[2]-y)^2)+((Wolf_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Wolf_mshape$X[10]-x)^2)+((Wolf_mshape$Y[10]-y)^2)+((Wolf_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Wolf_mshape$X[11]-x)^2)+((Wolf_mshape$Y[11]-y)^2)+((Wolf_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Wolf_mshape$X[12]-x)^2)+((Wolf_mshape$Y[12]-y)^2)+((Wolf_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Wolf_mshape$X[13]-x)^2)+((Wolf_mshape$Y[13]-y)^2)+((Wolf_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Wolf_mshape$X[14]-x)^2)+((Wolf_mshape$Y[14]-y)^2)+((Wolf_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Wolf_mshape$X[15]-x)^2)+((Wolf_mshape$Y[15]-y)^2)+((Wolf_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Wolf_mshape$X[16]-x)^2)+((Wolf_mshape$Y[16]-y)^2)+((Wolf_mshape$Z[16]-z)^2)),
                                  sqrt(((Wolf_mshape$X[17]-x)^2)+((Wolf_mshape$Y[17]-y)^2)+((Wolf_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); interanalyst_errors<-full_join(A2_inter_error_distances,
                                     A3_inter_error_distances); interanalyst_errors<-full_join(interanalyst_errors,
                                                                                                A1_inter_error_distances)
rm(A2_inter_error_distances, A3_inter_error_distances, A1_inter_error_distances)

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(interanalyst_errors$distance,
     main = "Histogram of Distribution of interanalyst Distances from Centroid",
     sub = if (
       shapiro.test(interanalyst_errors$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(interanalyst_errors$distance)$`p.value` > 0.05 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(interanalyst_errors$distance)$`p.value` > 0.001 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(interanalyst_errors$distance),
                    sd(interanalyst_errors$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(interanalyst_errors$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(interanalyst_errors$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of interanalyst Distances from Centroid",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(interanalyst_errors$distance)

# Shapiros
shapiro_values<-tibble(Sample = factor(), w.value = numeric(), p.value = numeric())
for (i in c("A1", "A2", "A3")){
  shapiro_values<-shapiro_values %>% add_row(Sample = i,
                                             p.value = shapiro.test(
                                               as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                             )$"p.value",
                                             w.value = shapiro.test(
                                               as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                             )$"statistic"[["W"]]
  )
}; shapiro_values

# NMAD and BWMV Calculations

interanalyst_landmark_errors<-tibble(Landmark = factor(),
                                     overall_median = numeric(),
                                     overall_NMAD = numeric(), overall_BWMV = numeric(),
                                     A3_median = numeric(),
                                     A3_NMAD = numeric(), A3_BWMV = numeric(),
                                     A1_median = numeric(),
                                     A1_NMAD = numeric(), A1_BWMV = numeric(),
                                     A2_median = numeric(),
                                     A2_NMAD = numeric(), A2_BWMV = numeric())
for (i in 1:length(levels(interanalyst_errors$Landmark))){
  interanalyst_landmark_error<-as.tibble(interanalyst_errors) %>% filter(Landmark == paste("LM", i, sep = ""))
  interanalyst_landmark_errors<-add_row(interanalyst_landmark_errors,
                                        Landmark = paste("LM", i, sep = ""),
                                        overall_median = median(interanalyst_landmark_error$distance),
                                        overall_NMAD = mad(interanalyst_landmark_error$distance, constant = 1.4826),
                                        overall_BWMV = sqrt(r.bw(interanalyst_landmark_error$distance)$"S.xx"[1]),
                                        A3_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance))),
                                        A3_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A3_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A1_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance))),
                                        A1_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A1_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A2_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance))),
                                        A2_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A2_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)))$"S.xx"[1]))
}; rm(interanalyst_landmark_error); interanalyst_landmark_errors

# Dog

A3_inter_error_distances<-tibble(Landmark = Dog_data$Landmark, Animal = Dog_data$Animal, Sample = Dog_data$Sample, x = Dog_data$x, y = Dog_data$y, z = Dog_data$z) %>%
  filter(Animal == "Dog") %>%
  filter(Sample == "A3") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Dog_mshape$X[1]-x)^2)+((Dog_mshape$Y[1]-y)^2)+((Dog_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Dog_mshape$X[3]-x)^2)+((Dog_mshape$Y[3]-y)^2)+((Dog_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Dog_mshape$X[4]-x)^2)+((Dog_mshape$Y[4]-y)^2)+((Dog_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Dog_mshape$X[5]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Dog_mshape$X[6]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Dog_mshape$X[7]-x)^2)+((Dog_mshape$Y[7]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Dog_mshape$X[8]-x)^2)+((Dog_mshape$Y[8]-y)^2)+((Dog_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Dog_mshape$X[10]-x)^2)+((Dog_mshape$Y[10]-y)^2)+((Dog_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Dog_mshape$X[11]-x)^2)+((Dog_mshape$Y[11]-y)^2)+((Dog_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Dog_mshape$X[12]-x)^2)+((Dog_mshape$Y[12]-y)^2)+((Dog_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Dog_mshape$X[13]-x)^2)+((Dog_mshape$Y[13]-y)^2)+((Dog_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Dog_mshape$X[14]-x)^2)+((Dog_mshape$Y[14]-y)^2)+((Dog_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Dog_mshape$X[15]-x)^2)+((Dog_mshape$Y[15]-y)^2)+((Dog_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Dog_mshape$X[16]-x)^2)+((Dog_mshape$Y[16]-y)^2)+((Dog_mshape$Z[16]-z)^2)),
                                  sqrt(((Dog_mshape$X[17]-x)^2)+((Dog_mshape$Y[17]-y)^2)+((Dog_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A1_inter_error_distances<-tibble(Landmark = Dog_data$Landmark, Animal = Dog_data$Animal, Sample = Dog_data$Sample, x = Dog_data$x, y = Dog_data$y, z = Dog_data$z) %>%
  filter(Animal == "Dog") %>%
  filter(Sample == "A1") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Dog_mshape$X[1]-x)^2)+((Dog_mshape$Y[1]-y)^2)+((Dog_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Dog_mshape$X[3]-x)^2)+((Dog_mshape$Y[3]-y)^2)+((Dog_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Dog_mshape$X[4]-x)^2)+((Dog_mshape$Y[4]-y)^2)+((Dog_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Dog_mshape$X[5]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Dog_mshape$X[6]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Dog_mshape$X[7]-x)^2)+((Dog_mshape$Y[7]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Dog_mshape$X[8]-x)^2)+((Dog_mshape$Y[8]-y)^2)+((Dog_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Dog_mshape$X[10]-x)^2)+((Dog_mshape$Y[10]-y)^2)+((Dog_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Dog_mshape$X[11]-x)^2)+((Dog_mshape$Y[11]-y)^2)+((Dog_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Dog_mshape$X[12]-x)^2)+((Dog_mshape$Y[12]-y)^2)+((Dog_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Dog_mshape$X[13]-x)^2)+((Dog_mshape$Y[13]-y)^2)+((Dog_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Dog_mshape$X[14]-x)^2)+((Dog_mshape$Y[14]-y)^2)+((Dog_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Dog_mshape$X[15]-x)^2)+((Dog_mshape$Y[15]-y)^2)+((Dog_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Dog_mshape$X[16]-x)^2)+((Dog_mshape$Y[16]-y)^2)+((Dog_mshape$Z[16]-z)^2)),
                                  sqrt(((Dog_mshape$X[17]-x)^2)+((Dog_mshape$Y[17]-y)^2)+((Dog_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); A2_inter_error_distances<-tibble(Landmark = Dog_data$Landmark, Animal = Dog_data$Animal, Sample = Dog_data$Sample, x = Dog_data$x, y = Dog_data$y, z = Dog_data$z) %>%
  filter(Animal == "Dog") %>%
  filter(Sample == "A2") %>%
  mutate(distance = if_else(
    Landmark == "LM1",
    sqrt(((Dog_mshape$X[1]-x)^2)+((Dog_mshape$Y[1]-y)^2)+((Dog_mshape$Z[1]-z)^2)),
    if_else(
      Landmark == "LM2",
      sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
      if_else(
        Landmark == "LM3",
        sqrt(((Dog_mshape$X[3]-x)^2)+((Dog_mshape$Y[3]-y)^2)+((Dog_mshape$Z[3]-z)^2)),
        if_else(
          Landmark == "LM4",
          sqrt(((Dog_mshape$X[4]-x)^2)+((Dog_mshape$Y[4]-y)^2)+((Dog_mshape$Z[4]-z)^2)),
          if_else(
            Landmark == "LM5",
            sqrt(((Dog_mshape$X[5]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
            if_else(
              Landmark == "LM6",
              sqrt(((Dog_mshape$X[6]-x)^2)+((Dog_mshape$Y[6]-y)^2)+((Dog_mshape$Z[6]-z)^2)),
              if_else(
                Landmark == "LM7",
                sqrt(((Dog_mshape$X[7]-x)^2)+((Dog_mshape$Y[7]-y)^2)+((Dog_mshape$Z[7]-z)^2)),
                if_else(
                  Landmark == "LM8",
                  sqrt(((Dog_mshape$X[8]-x)^2)+((Dog_mshape$Y[8]-y)^2)+((Dog_mshape$Z[8]-z)^2)),
                  if_else(
                    Landmark == "LM9",
                    sqrt(((Dog_mshape$X[2]-x)^2)+((Dog_mshape$Y[2]-y)^2)+((Dog_mshape$Z[2]-z)^2)),
                    if_else(
                      Landmark == "LM10",
                      sqrt(((Dog_mshape$X[10]-x)^2)+((Dog_mshape$Y[10]-y)^2)+((Dog_mshape$Z[10]-z)^2)),
                      if_else(
                        Landmark == "LM11",
                        sqrt(((Dog_mshape$X[11]-x)^2)+((Dog_mshape$Y[11]-y)^2)+((Dog_mshape$Z[11]-z)^2)),
                        if_else(
                          Landmark == "LM12",
                          sqrt(((Dog_mshape$X[12]-x)^2)+((Dog_mshape$Y[12]-y)^2)+((Dog_mshape$Z[12]-z)^2)),
                          if_else(
                            Landmark == "LM13",
                            sqrt(((Dog_mshape$X[13]-x)^2)+((Dog_mshape$Y[13]-y)^2)+((Dog_mshape$Z[13]-z)^2)),
                            if_else(
                              Landmark == "LM14",
                              sqrt(((Dog_mshape$X[14]-x)^2)+((Dog_mshape$Y[14]-y)^2)+((Dog_mshape$Z[14]-z)^2)),
                              if_else(
                                Landmark == "LM15",
                                sqrt(((Dog_mshape$X[15]-x)^2)+((Dog_mshape$Y[15]-y)^2)+((Dog_mshape$Z[15]-z)^2)),
                                if_else(
                                  Landmark == "LM16",
                                  sqrt(((Dog_mshape$X[16]-x)^2)+((Dog_mshape$Y[16]-y)^2)+((Dog_mshape$Z[16]-z)^2)),
                                  sqrt(((Dog_mshape$X[17]-x)^2)+((Dog_mshape$Y[17]-y)^2)+((Dog_mshape$Z[17]-z)^2))
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )); interanalyst_errors<-full_join(A2_inter_error_distances,
                                     A3_inter_error_distances); interanalyst_errors<-full_join(interanalyst_errors,
                                                                                                A1_inter_error_distances)
rm(A2_inter_error_distances, A3_inter_error_distances, A1_inter_error_distances)

# normality plots

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(interanalyst_errors$distance,
     main = "Histogram of Distribution of interanalyst Distances from Centroid",
     sub = if (
       shapiro.test(interanalyst_errors$distance)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(interanalyst_errors$distance)$`p.value` > 0.05 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(interanalyst_errors$distance)$`p.value` > 0.001 && shapiro.test(interanalyst_errors$distance)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Distribution of Distances from Centroid"
     ); curve(dnorm(x,
                    mean(interanalyst_errors$distance),
                    sd(interanalyst_errors$distance)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(interanalyst_errors$distance),
                                 col = "red",
                                 lwd = 2); qqPlot(interanalyst_errors$distance,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Distribution of interanalyst Distances from Centroid",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(interanalyst_errors$distance)

# Shapiros
shapiro_values<-tibble(Sample = factor(), w.value = numeric(), p.value = numeric())
for (i in c("A1", "A2", "A3")){
  shapiro_values<-shapiro_values %>% add_row(Sample = i,
                                             p.value = shapiro.test(
                                               as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                             )$"p.value",
                                             w.value = shapiro.test(
                                               as.matrix(as.tibble(interanalyst_errors) %>% filter(Sample == i) %>% select(distance))
                                             )$"statistic"[["W"]]
  )
}; shapiro_values

# NMAD and BWMV Calculations

interanalyst_landmark_errors<-tibble(Landmark = factor(),
                                     overall_median = numeric(),
                                     overall_NMAD = numeric(), overall_BWMV = numeric(),
                                     A3_median = numeric(),
                                     A3_NMAD = numeric(), A3_BWMV = numeric(),
                                     A1_median = numeric(),
                                     A1_NMAD = numeric(), A1_BWMV = numeric(),
                                     A2_median = numeric(),
                                     A2_NMAD = numeric(), A2_BWMV = numeric())
for (i in 1:length(levels(interanalyst_errors$Landmark))){
  interanalyst_landmark_error<-as.tibble(interanalyst_errors) %>% filter(Landmark == paste("LM", i, sep = ""))
  interanalyst_landmark_errors<-add_row(interanalyst_landmark_errors,
                                        Landmark = paste("LM", i, sep = ""),
                                        overall_median = median(interanalyst_landmark_error$distance),
                                        overall_NMAD = mad(interanalyst_landmark_error$distance, constant = 1.4826),
                                        overall_BWMV = sqrt(r.bw(interanalyst_landmark_error$distance)$"S.xx"[1]),
                                        A3_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance))),
                                        A3_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A3_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A3") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A1_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance))),
                                        A1_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A1_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A1") %>%
                                            select(distance)))$"S.xx"[1]),
                                        A2_median = median(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance))),
                                        A2_NMAD = mad(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)),
                                          constant = 1.4826),
                                        A2_BWMV = sqrt(r.bw(as.matrix(
                                          as.tibble(interanalyst_landmark_error) %>%
                                            filter(Sample == "A2") %>%
                                            select(distance)))$"S.xx"[1]))
}; rm(interanalyst_landmark_error); interanalyst_landmark_errors

#

# Distribution Analysis ----------------------------------------

# Create plots of distributions

Landmark<-split(data, data$Landmark)
for (i in 1:length(levels(data$Landmark))){
  LM_number<-paste("LM", i, sep ="")
  Landmark_dist<-Landmark[[LM_number]]
  landmark_plot<-ggplot(Landmark_dist, aes(x = x, y = y, colour = Sample)) +
    geom_point(stat = "identity", size = 1.5) +
    stat_ellipse(size = 0.5) +
    scale_color_manual(values = c("A3" = "BA3k", "A1" = "Red", "A2" = "Blue")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = "NA"),
      plot.margin = unit(c(1,1,1,1), "cm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(colour = "NA", fill = "transparent"),
      legend.position = "none"
    ) +
    coord_cartesian(xlim = c(min(data$x)-0.01, max(data$x)+0.01),
                    ylim = c(min(data$y)-0.01, max(data$y)+0.01))
  ggsave(paste(directory, "\\Plots of Distributions\\",
               LM_number, ".svg", sep = ""), plot = landmark_plot,
         device = "svg",
         width = 50,
         height = 50,
         units = "cm",
         dpi = 1000,
         bg = "transparent")
}

# Calculate Significance of Variations in Distribution

LM_number<-"LM17"
Landmark_dist<-Landmark[[LM_number]]; Landmark_dist_coords<-as.matrix(as.tibble(Landmark_dist) %>% select(x, y, z))
pairwise.perm.manova(Landmark_dist_coords, Landmark_dist$Sample,
                     test = c("Wilks"), nperm = 999,
                     progress = TRUE, p.method = "none")

# Animal

# Wolf

Landmark<-split(Wolf_data, Wolf_data$Landmark)
for (i in 1:length(levels(Wolf_data$Landmark))){
  LM_number<-paste("LM", i, sep ="")
  Landmark_dist<-Landmark[[LM_number]]
  landmark_plot<-ggplot(Landmark_dist, aes(x = x, y = y, colour = Sample)) +
    geom_point(stat = "identity", size = 1.5) +
    stat_ellipse(size = 0.5) +
    scale_color_manual(values = c("A3" = "BA3k", "A1" = "Red", "A2" = "Blue")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = "NA"),
      plot.margin = unit(c(1,1,1,1), "cm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(colour = "NA", fill = "transparent"),
      legend.position = "none"
    ) +
    coord_cartesian(xlim = c(min(Wolf_data$x)-0.01, max(Wolf_data$x)+0.01),
                    ylim = c(min(Wolf_data$y)-0.01, max(Wolf_data$y)+0.01))
  ggsave(paste(directory, "\\Plots of Distributions\\",
               LM_number, ".svg", sep = ""), plot = landmark_plot,
         device = "svg",
         width = 50,
         height = 50,
         units = "cm",
         dpi = 1000,
         bg = "transparent")
}

# Multivariate Analyses

LM_number<-"LM17"
Landmark_dist<-Landmark[[LM_number]]; Landmark_dist_coords<-as.matrix(as.tibble(Landmark_dist) %>% select(x, y, z))
pairwise.perm.manova(Landmark_dist_coords, Landmark_dist$Sample,
                     test = c("Wilks"), nperm = 999,
                     progress = TRUE, p.method = "none")

# Dog

Landmark<-split(Dog_data, Dog_data$Landmark)
for (i in 1:length(levels(Dog_data$Landmark))){
  LM_number<-paste("LM", i, sep ="")
  Landmark_dist<-Landmark[[LM_number]]
  landmark_plot<-ggplot(Landmark_dist, aes(x = x, y = y, colour = Sample)) +
    geom_point(stat = "identity", size = 1.5) +
    stat_ellipse(size = 0.5) +
    scale_color_manual(values = c("A3" = "BA3k", "A1" = "Red", "A2" = "Blue")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", colour = "NA"),
      plot.margin = unit(c(1,1,1,1), "cm"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.border = element_rect(colour = "NA", fill = "transparent"),
      legend.position = "none"
    ) +
    coord_cartesian(xlim = c(min(Dog_data$x)-0.01, max(Dog_data$x)+0.01),
                    ylim = c(min(Dog_data$y)-0.01, max(Dog_data$y)+0.01))
  ggsave(paste(directory, "\\Plots of Distributions\\",
               LM_number, ".svg", sep = ""), plot = landmark_plot,
         device = "svg",
         width = 50,
         height = 50,
         units = "cm",
         dpi = 1000,
         bg = "transparent")
}

# Multivariate Analyses

LM_number<-"LM17"
Landmark_dist<-Landmark[[LM_number]]; Landmark_dist_coords<-as.matrix(as.tibble(Landmark_dist) %>% select(x, y, z))
pairwise.perm.manova(Landmark_dist_coords, Landmark_dist$Sample,
                     test = c("Wilks"), nperm = 999,
                     progress = TRUE, p.method = "none")

#

# Procrustes Distances ##########################

# Analyst -----------------------------------

Proc_dist<-ShapeDist(shapes = GPA$coords, reference = mshape(GPA$coords))
Proc_distances<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(Proc_dist)) {
  Proc_distances<-add_row(Proc_distances,
                          label = rownames(a$labels)[i],
                          sample = a$labels[i],
                          animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                          proc_dist = Proc_dist[i])
  }

# Overall normality of procrustes distances

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(Proc_distances$proc_dist,
     main = "Histogram of Procrustes Distances",
     sub = if (
       shapiro.test(Proc_distances$proc_dist)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(Proc_distances$proc_dist)$`p.value` > 0.05 && shapiro.test(Proc_distances$proc_dist)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(Proc_distances$proc_dist)$`p.value` > 0.001 && shapiro.test(Proc_distances$proc_dist)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Procrustes Distances"
     ); curve(dnorm(x,
                    mean(Proc_distances$proc_dist),
                    sd(Proc_distances$proc_dist)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(Proc_distances$proc_dist),
                                 col = "red",
                                 lwd = 2); qqPlot(Proc_distances$proc_dist,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Procrustes Distances",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")

par(op)
shapiro.test(Proc_distances$proc_dist)

# Proc dist normalities according to ... variables

analyst<-split(Proc_distances, Proc_distances$sample)
Proc_A1<-analyst$A1
Proc_A2<-analyst$A2
Proc_A3<-analyst$A3

test<-Proc_A2$proc_dist
shapiro.test(test)
mean(test)
median(test)
sd(test)
skewness(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(test,
     main = "Histogram of Procrustes Distances",
     sub = if (
       shapiro.test(test)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(test)$`p.value` > 0.05 && shapiro.test(test)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(test)$`p.value` > 0.001 && shapiro.test(test)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Procrustes Distances"
     ); curve(dnorm(x,
                    mean(test),
                    sd(test)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(test),
                                 col = "red",
                                 lwd = 2); qqPlot(test,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Procrustes Distances",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles")
par(op)

compare<-rbind(Proc_A2, Proc_A3)
kruskal.test(proc_dist~sample, data = compare)

distances_for_plot<-tibble(label = Proc_distances$label, distance = Proc_distances$proc_dist) %>%
  as.data.frame(); rownames(distances_for_plot)<-distances_for_plot$label; distances_for_plot<-distances_for_plot[-c(1)]
agn<-agnes(distances_for_plot, method = "average")
plot(agn, which.plot = 2, cex = 0.65, main = "UPGMA Tree")

boxplot(proc_dist~sample*animal, data = Proc_distances, col =  c("red","blue","green"))

op<-par(mfrow = c(1,2))
boxplot(proc_dist~animal, data = Proc_distances); boxplot(proc_dist~sample, data = Proc_distances)
par(op)

# Repeatability measure - proportion of variance due to true variation

((sum(as.matrix(tibble(s2a = Proc_A1$proc_dist^2))))/
    (nrow(Proc_distances)-length(levels(Proc_distances$sample))))/
  ((sum(as.matrix(tibble(s2a = Proc_distances$proc_dist^2))))/
     (nrow(Proc_distances)-1))*100

((sum(as.matrix(tibble(s2a = Proc_A2$proc_dist^2))))/
    (nrow(Proc_distances)-length(levels(Proc_distances$sample))))/
  ((sum(as.matrix(tibble(s2a = Proc_distances$proc_dist^2))))/
     (nrow(Proc_distances)-1))*100

((sum(as.matrix(tibble(s2a = Proc_A3$proc_dist^2))))/
    (nrow(Proc_distances)-length(levels(Proc_distances$sample))))/
  ((sum(as.matrix(tibble(s2a = Proc_distances$proc_dist^2))))/
     (nrow(Proc_distances)-1))*100


#

# Only with wolves ---------------------------

Proc_dist_wolf<-ShapeDist(shapes = Wolf_GPA$coords, reference = mshape(Wolf_GPA$coords))
Proc_distances_wolf<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(Proc_dist_wolf)) {
  Proc_distances_wolf<-add_row(Proc_distances_wolf,
                               label = rownames(Wolf$labels)[i],
                               sample = Wolf$labels[i],
                               animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                               proc_dist = Proc_dist_wolf[i])
}

# Overall normality of procrustes distances

shapiro.test(Proc_distances_wolf$proc_dist)

# Proc dist normalities according to ... variables

# Analyst

analyst<-split(Proc_distances_wolf, Proc_distances_wolf$sample)
Proc_A1<-analyst$A1
Proc_A2<-analyst$A2
Proc_A3<-analyst$A3

test<-Proc_A1$proc_dist
shapiro.test(test)
mean(test)
median(test)
sd(test)
skewness(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(test,
     main = "Histogram of Procrustes Distances",
     sub = if (
       shapiro.test(test)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(test)$`p.value` > 0.05 && shapiro.test(test)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(test)$`p.value` > 0.001 && shapiro.test(test)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Procrustes Distances"
     ); curve(dnorm(x,
                    mean(test),
                    sd(test)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(test),
                                 col = "red",
                                 lwd = 2); qqPlot(test,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Procrustes Distances",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles"); par(op)

compare<-rbind(Proc_A3, Proc_A1)
kruskal.test(proc_dist~sample, data = compare)

par(mfrow = c(1,1))
boxplot(proc_dist~sample, data = Proc_distances_wolf)

# Dog -------------------------

Proc_dist_Dog<-ShapeDist(shapes = Dog_GPA$coords, reference = mshape(Dog_GPA$coords))
Proc_distances_Dog<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(Proc_dist_Dog)) {
  Proc_distances_Dog<-add_row(Proc_distances_Dog,
                                label = rownames(Dog$labels)[i],
                                sample = Dog$labels[i],
                                animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                                proc_dist = Proc_dist_Dog[i])
}

# Overall normality of procrustes distances

shapiro.test(Proc_distances_Dog$proc_dist)

# Proc dist normalities according to ... variables

# Analyst

analyst<-split(Proc_distances_Dog, Proc_distances_Dog$sample)
Proc_A1<-analyst$A1
Proc_A2<-analyst$A2
Proc_A3<-analyst$A3

test<-Proc_A3$proc_dist
shapiro.test(test)
mean(test)
median(test)
sd(test)
skewness(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

op<-par(mfrow = c(1,2)) # mar = bottom, left, top, right
hist(test,
     main = "Histogram of Procrustes Distances",
     sub = if (
       shapiro.test(test)$`p.value` > 0.5
       ) {"Shapiro test: Significant Gaussian Distribution"} else if (
         shapiro.test(test)$`p.value` > 0.05 && shapiro.test(test)$`p.value` < 0.5
         ) {"Shapiro test: Guassian Distribution"} else if (
           shapiro.test(test)$`p.value` > 0.001 && shapiro.test(test)$`p.value` < 0.05
           ) {"Shapiro test: Non-Gaussian Distribution"} else {"Shapiro test: V. Significant Non-Gaussian Distribution"},
     col.sub = "red",
     font.sub = 2,
     col = "grey",
     prob = TRUE,
     xlab = "Procrustes Distances"
     ); curve(dnorm(x,
                    mean(test),
                    sd(test)),
              col = "blue",
              lwd = 2,
              add = TRUE); lines(density(test),
                                 col = "red",
                                 lwd = 2); qqPlot(test,
                                                  envelope = 0.95, col.lines = "red",
                                                  pch = 16, cex = 0.5,
                                                  id = FALSE,
                                                  main = "Procrustes Distances",
                                                  ylab = "Actual Quantiles",
                                                  xlab = "Theoretical Quantiles"); par(op)

compare<-rbind(Proc_A1, Proc_A2)
kruskal.test(proc_dist~sample, data = compare)

par(mfrow = c(1,1))
boxplot(proc_dist~sample, data = Proc_distances_Dog)

# Indiviual observer ------------------------

# A3

A3_dist<-ShapeDist(shapes = A3_GPA$coords, reference = mshape(A3_GPA$coords))
A3_distances<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(A3_dist)) {
  A3_distances<-add_row(A3_distances,
                         label = rownames(A3$labels)[i],
                         sample = A3$labels[i],
                         animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                         proc_dist = A3_dist[i])
}

analyst<-split(A3_distances, A3_distances$animal)
Proc_Dog<-analyst$Dog
Proc_Wolf<-analyst$Wolf

test<-Proc_Dog$proc_dist
shapiro.test(test)
mean(test)
sd(test)
skewness(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

# not normal
kruskal.test(proc_dist~animal, data = A3_distances)

boxplot(proc_dist~animal, data = A3_distances)

# A1

A1_dist<-ShapeDist(shapes = A1_GPA$coords, reference = mshape(A1_GPA$coords))
A1_distances<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(A1_dist)) {
  A1_distances<-add_row(A1_distances,
                         label = rownames(A1$labels)[i],
                         sample = A1$labels[i],
                         animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                         proc_dist = A1_dist[i])
}

analyst<-split(A1_distances, A1_distances$animal)
Proc_Dog<-analyst$Dog
Proc_Wolf<-analyst$Wolf

test<-Proc_Wolf$proc_dist
shapiro.test(test)
mean(test)
sd(test)
median(test)
skewness(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

# not normal
kruskal.test(proc_dist~animal, data = A1_distances)

boxplot(proc_dist~animal, data = A1_distances)

# A2

A2_dist<-ShapeDist(shapes = A2_GPA$coords, reference = mshape(A2_GPA$coords))
A2_distances<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(A2_dist)) {
  A2_distances<-add_row(A2_distances,
                          label = rownames(A2$labels)[i],
                          sample = A2$labels[i],
                          animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                          proc_dist = A2_dist[i])
}

analyst<-split(A2_distances, A2_distances$animal)
Proc_Dog<-analyst$Dog
Proc_Wolf<-analyst$Wolf

test<-Proc_Dog$proc_dist
shapiro.test(test)
mean(test)
sd(test)
skewness(test)
median(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

# not normal
kruskal.test(proc_dist~animal, data = A2_distances)

boxplot(proc_dist~animal, data = A2_distances)

