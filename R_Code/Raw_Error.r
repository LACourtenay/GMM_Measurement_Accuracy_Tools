# Libraries ----------------------

library(geomorph)
library(shapes)
library(tibble)
library(dplyr)
library(tidyverse)
library(spatstat)
library(MVN)
library(car)
library(RVAideMemoire)
library(Evomorph)
library(e1071)
library(cluster)
library(asbio) # If errors occur with this package consult ReadMe.

# Obtaining Data ---------------

directory<-"C:\\"

a<-read.morphologika(paste(directory, "All_Landmarks.txt", sep = "\\"))
A1<-read.morphologika(paste(directory, "A1.txt", sep = "\\"))
A2<-read.morphologika(paste(directory, "A2.txt", sep = "\\"))
A3<-read.morphologika(paste(directory, "A3.txt", sep = "\\"))
Wolf<-read.morphologika(paste(directory, "Wolf.txt", sep = "\\"))
Dog<-read.morphologika(paste(directory, "Dog.txt", sep = "\\"))

# Extract Data ------------------------------

raw_data<-tibble(label = character(), Landmark = factor(), Sample = factor(), Animal = factor(),
                 x = numeric(), y = numeric(), z = numeric())
for (i in 1:length(a$labels)){
  sep<-a$coords[,,i]
  for (i2 in 1:nrow(sep)){
    LM<-sep[i2,]
    raw_data<-add_row(raw_data, label = rownames(a$labels)[i],
                      Landmark = paste("LM", i2, sep = ""),
                      Sample = a$labels[i],
                      Animal = if_else(grepl("Wolf", label) == TRUE, "Wolf", "Dog"),
                      x = LM[1], y = LM[2], z = LM[3])
  }
}; raw_data<-as.tibble(raw_data) %>%
  mutate(label = if_else(grepl("Wolf", label) == TRUE,
                         substr(label, 1, 11), substr(label, 1, 10))); raw_data$label<-factor(raw_data$label); rm(sep, LM)
distances<-tibble(label = character(), Landmark = factor(), Animal = factor(),
                  A1_x = numeric(), A1_y = numeric(), A1_z = numeric(),
                  A3_x = numeric(), A3_y = numeric(), A3_z = numeric(),
                  A2_x = numeric(), A2_y = numeric(), A2_z = numeric())
split_label<-split(raw_data, raw_data$label)
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

ind_split<-split(distances, distances$label)

for (i in 1:length(ind_split)) {print(shapiro.test(as.matrix(ind_split[[i]][,13:15]))$p.value)}

# Error analyses ====================================

# preprocessing of data --------------------------

# Overall errors per individual studied
ind_error<-tibble(label = character(), Animal = factor(),
                  A1_A3 = numeric (), A1_A2 = numeric(), A3_A2 = numeric())
for (i in 1:length(ind_split)) {
  ind_error<-add_row(ind_error, label = ind_split[[i]]$label[1],
                     Animal = ind_split[[i]]$Animal[1],
                     #A1_A3 = sqrt(r.bw(as.matrix(ind_split[[i]][,13]))$"S.xx"[1]),
                     #A1_A2 = sqrt(r.bw(as.matrix(ind_split[[i]][,14]))$"S.xx"[1]),
                     #A3_A2 = sqrt(r.bw(as.matrix(ind_split[[i]][,15]))$"S.xx"[1]))
                     A1_A3 = median(as.matrix(ind_split[[i]][,13])),
                     A1_A2 = median(as.matrix(ind_split[[i]][,14])),
                     A3_A2 = median(as.matrix(ind_split[[i]][,15]))
                     )
}

# seperate data according to landmarks under observation:

new_distances<-distances %>% filter(Landmark != "LM6" &
                                      Landmark != "LM7" &
                                      Landmark != "LM8" &
                                      Landmark != "LM9" &
                                      Landmark != "LM10" &
                                      Landmark != "LM11" &
                                      Landmark != "LM12" &
                                      Landmark != "LM13" &
                                      Landmark != "LM14" &
                                      Landmark != "LM15" &
                                      Landmark != "LM16" &
                                      Landmark != "LM17")
new_split<-split(new_distances, new_distances$label)
new_error<-tibble(label = character(), Animal = factor(),
                  A1_A3 = numeric (), A1_A2 = numeric(), A3_A2 = numeric())
for (i in 1:length(new_split)) {
  new_error<-add_row(new_error, label = new_split[[i]]$label[1],
                     Animal = new_split[[i]]$Animal[1],
                     A1_A3 = sqrt(r.bw(as.matrix(new_split[[i]][,13]))$"S.xx"[1]),
                     A1_A2 = sqrt(r.bw(as.matrix(new_split[[i]][,14]))$"S.xx"[1]),
                     A3_A2 = sqrt(r.bw(as.matrix(new_split[[i]][,15]))$"S.xx"[1]))
}

# Save this data to reload into R

write.table(new_error, paste(directory, "animal_error.csv", sep = "\\"))

# Load data of interest

animal_error<-read.csv(paste(directory,"\\Animal_Error.csv",sep = ""), sep = ",")

# Descriptive analysis -----------------------

boxplot(Distance~Animal, data = animal_error, outline = FALSE)
kruskal.test(Distance~Animal, data = animal_error)

min(animal_error$Distance)
mean(animal_error$Distance)
sd(animal_error$Distance)
median(animal_error$Distance)
mad(animal_error$Distance, constant = 1.4826)
sqrt(r.bw(animal_error$Distance)$"S.xx"[1])
max(animal_error$Distance)

min(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)
mean(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)
sd(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)
median(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)
mad(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance, constant = 1.4826)
sqrt(r.bw(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)$"S.xx"[1])
max(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)

skewness(split(animal_error, animal_error$Animal)[["Wolf"]]$Distance)

min(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)
mean(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)
sd(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)
median(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)
mad(split(animal_error, animal_error$Animal)[["Dog"]]$Distance, constant = 1.4826)
sqrt(r.bw(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)$"S.xx"[1])
max(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)

skewness(split(animal_error, animal_error$Animal)[["Dog"]]$Distance)

# Statistical analysis ---------------

shapiro.test(animal_error$Distance)
kruskal.test(Distance~Ind, data = animal_error)

shapiro.test(rbind(split(animal_error, animal_error$Ind)[["A1_A3"]],
                   split(animal_error, animal_error$Ind)[["A1_A2"]])[["Distance"]])
kruskal.test(Distance~Ind, data = rbind(split(animal_error, animal_error$Ind)[["A1_A3"]],
                                        split(animal_error, animal_error$Ind)[["A1_A2"]]))

shapiro.test(rbind(split(animal_error, animal_error$Ind)[["A3_A2"]],
                   split(animal_error, animal_error$Ind)[["A1_A2"]])[["Distance"]])
kruskal.test(Distance~Ind, data = rbind(split(animal_error, animal_error$Ind)[["A3_A2"]],
                                        split(animal_error, animal_error$Ind)[["A1_A2"]]))

shapiro.test(rbind(split(animal_error, animal_error$Ind)[["A3_A2"]],
                   split(animal_error, animal_error$Ind)[["A1_A3"]])[["Distance"]])
kruskal.test(Distance~Ind, data = rbind(split(animal_error, animal_error$Ind)[["A3_A2"]],
                                        split(animal_error, animal_error$Ind)[["A1_A3"]]))

boxplot(Distance~Ind, data = animal_error, outline = FALSE)
boxplot(Distance~Ind+Animal, data = animal_error, outline = FALSE)

# Aditional analysis -----------------------------

lm_error_split<-split(distances, distances$Landmark)
lm_error<-tibble(LM = factor(),
                 A1_A3 = numeric(), A1_A2 = numeric(), A3_A2 = numeric())
for (i in 1:length(lm_error_split)) {
  lm_error<-add_row(lm_error, LM = lm_error_split[[i]]$Landmark[1],
                    A1_A3 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,13]))$"S.xx"[1]),
                    A1_A2 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,14]))$"S.xx"[1]),
                    A3_A2 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,15]))$"S.xx"[1]))
}

lm_error_overall<-tibble(LM = factor(), BWMV = numeric(),
                         NMAD = numeric(), Median = numeric())
for (i in 1:length(lm_error_split)) {
  lm_error_overall<-add_row(lm_error_overall, LM = lm_error_split[[i]]$Landmark[1],
                            BWMV = sqrt(r.bw(as.matrix(lm_error_split[[i]][,13:15]))$"S.xx"[1]),
                            NMAD = mad(as.matrix(lm_error_split[[i]][,13:15]), constant = 1.4826),
                            Median = median(as.matrix(lm_error_split[[i]][,13:15])))
}

lm_error_split<-split(distances, distances$Landmark)
lm_error<-tibble(LM = factor(),
                 A1_A3 = numeric(), A1_A2 = numeric(), A3_A2 = numeric())
for (i in 1:length(lm_error_split)) {
  lm_error<-add_row(lm_error, LM = lm_error_split[[i]]$Landmark[1],
                    A1_A3 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,13]))$"S.xx"[1]),
                    A1_A2 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,14]))$"S.xx"[1]),
                    A3_A2 = sqrt(r.bw(as.matrix(lm_error_split[[i]][,15]))$"S.xx"[1]))
}
