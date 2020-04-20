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

absolute_proc<-procGPA(a$coords, scale = FALSE)
A3_proc<-procGPA(A3$coords, scale = FALSE)
A1_proc<-procGPA(A1$coords, scale = FALSE)
A2_proc<-procGPA(A2$coords, scale = FALSE)
Wolf_proc<-procGPA(Wolf$coords, scale = FALSE)
Dog_proc<-procGPA(Dog$coords, scale = FALSE)

absolute_mshape<-data.frame(mshape(absolute_proc$rotated));names(absolute_mshape)[1]<-"X";names(absolute_mshape)[2]<-"Y";names(absolute_mshape)[3]<-"Z"
A3_mshape<-data.frame(mshape(A3_proc$rotated));names(A3_mshape)[1]<-"X";names(A3_mshape)[2]<-"Y";names(A3_mshape)[3]<-"Z"
A1_mshape<-data.frame(mshape(A1_proc$rotated));names(A1_mshape)[1]<-"X";names(A1_mshape)[2]<-"Y";names(A1_mshape)[3]<-"Z"
A2_mshape<-data.frame(mshape(A2_proc$rotated));names(A2_mshape)[1]<-"X";names(A2_mshape)[2]<-"Y";names(A2_mshape)[3]<-"Z"
Wolf_mshape<-data.frame(mshape(Wolf_proc$rotated));names(Wolf_mshape)[1]<-"X";names(Wolf_mshape)[2]<-"Y";names(Wolf_mshape)[3]<-"Z"
Dog_mshape<-data.frame(mshape(Dog_proc$rotated));names(Dog_mshape)[1]<-"X";names(Dog_mshape)[2]<-"Y";names(Dog_mshape)[3]<-"Z"

# Obtaining interanalyst data - Seperating by Animals ----------------

# Wolf

Wolf_A1<-read.morphologika(paste(directory, "\\Wolf_A1.txt", sep = ""))
Wolf_A2<-read.morphologika(paste(directory, "\\Wolf_A2.txt", sep = ""))
Wolf_A3<-read.morphologika(paste(directory, "\\Wolf_A3.txt", sep = ""))
Wolf_A1_mshape<-data.frame(mshape(Wolf_A1$coords))
Wolf_A2_mshape<-data.frame(mshape(Wolf_A2$coords))
Wolf_A3_mshape<-data.frame(mshape(Wolf_A3$coords))

# Dog

Dog_A1<-read.morphologika(paste(directory, "\\Dog_A1.txt", sep = ""))
Dog_A2<-read.morphologika(paste(directory, "\\Dog_A2.txt", sep = ""))
Dog_A3<-read.morphologika(paste(directory, "\\Dog_A3.txt", sep = ""))
Dog_A1_mshape<-data.frame(mshape(Dog_A1$coords))
Dog_A2_mshape<-data.frame(mshape(Dog_A2$coords))
Dog_A3_mshape<-data.frame(mshape(Dog_A3$coords))

#

# Procrustes Distances ##########################

# Analyst -----------------------------------

Proc_dist<-ShapeDist(shapes = absolute_proc$rotated, reference = mshape(absolute_proc$rotated))
Proc_distances<-tibble(label = character(), sample = factor(), animal = factor(), proc_dist = numeric())
for (i in 1:length(Proc_dist)) {
  Proc_distances<-add_row(Proc_distances,
                          label = rownames(a$labels)[i],
                          sample = a$labels[i],
                          animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                          proc_dist = Proc_dist[i])
  }

# UPGMA Tree

distances_for_plot<-tibble(label = Proc_distances$label, distance = Proc_distances$proc_dist) %>%
  as.data.frame(); rownames(distances_for_plot)<-distances_for_plot$label; distances_for_plot<-distances_for_plot[-c(1)]
agn<-agnes(distances_for_plot, method = "average")
plot(agn, which.plot = 2, cex = 0.65)

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

test<-Proc_A3$proc_dist
shapiro.test(test)
mean(test)
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

compare<-rbind(Proc_A1, Proc_A2)
kruskal.test(proc_dist~sample, data = compare)

distances_for_plot<-tibble(label = Proc_distances$label, distance = Proc_distances$proc_dist) %>%
  as.data.frame(); rownames(distances_for_plot)<-distances_for_plot$label; distances_for_plot<-distances_for_plot[-c(1)]
agn<-agnes(distances_for_plot, method = "average")
plot(agn, which.plot = 2, cex = 0.65)

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

Proc_dist_wolf<-ShapeDist(shapes = Wolf_proc$rotated, reference = mshape(Wolf_proc$rotated))
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

# Analyst

analyst<-split(Proc_distances_wolf, Proc_distances_wolf$sample)
Proc_A1<-analyst$A1
Proc_A2<-analyst$A2
Proc_A3<-analyst$A3

test<-Proc_A3$proc_dist
shapiro.test(test)
mean(test)
sd(test)
median(test)
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
boxplot(proc_dist~sample, data = Proc_distances_wolf)

# Dog -------------------------

Proc_dist_Dog<-ShapeDist(shapes = Dog_proc$rotated, reference = mshape(Dog_proc$rotated))
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
sd(test)
median(test)
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
boxplot(proc_dist~sample, data = Proc_distances_Dog)

# Indiviual observer ------------------------

# A3

A3_dist<-ShapeDist(shapes = A3_proc$rotated, reference = mshape(A3_proc$rotated))
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

test<-Proc_Wolf$proc_dist
shapiro.test(test)
mean(test)
sd(test)
skewness(test)
median(test)
kurtosis(test, type = 1)
mad(test, constant = 1.4826)
sqrt(r.bw(test)$"S.xx"[1])

# not normal
kruskal.test(proc_dist~animal, data = A3_distances)

boxplot(proc_dist~animal, data = A3_distances)

# A1

A1_dist<-ShapeDist(shapes = A1_proc$rotated, reference = mshape(A1_proc$rotated))
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

test<-Proc_Dog$proc_dist
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

A2_dist<-ShapeDist(shapes = A2_proc$rotated, reference = mshape(A2_proc$rotated))
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

test<-Proc_Wolf$proc_dist
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
