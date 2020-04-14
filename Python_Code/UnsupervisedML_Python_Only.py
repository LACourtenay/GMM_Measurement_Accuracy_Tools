import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sklearn
from sklearn.cluster import DBSCAN, KMeans, MeanShift, AgglomerativeClustering, estimate_bandwidth
import scipy
from scipy.spatial.distance import cdist
from scipy import stats
print("Imported Packages:", "\nNumpy: {}.".format(np.__version__),
      "\nMatplotlib: {}.".format(matplotlib.__version__),
      "\nSciPy: {}.".format(scipy.__version__),
      "\nScikit-learn: {}.".format(sklearn.__version__),
      "\nPandas: {}.".format(pd.__version__))

# set internal matplotlib parameters

plt.style.use("ggplot")
plt.rcParams["figure.figsize"] = (16,12)
plt.rcParams["axes.facecolor"] = "None"

# load data

os.chdir("C:/Users/ladc1/Desktop/PhD/Methodological Research")
data = pd.read_csv("proc_coords.txt", sep = "\t")
all_data = data[["Landmark","x","y","z"]]

outline_data = all_data[all_data["Landmark"] != "LM14"]
outline_data = outline_data[outline_data["Landmark"] != "LM15"]
outline_data = outline_data[outline_data["Landmark"] != "LM16"]
outline_data = outline_data[outline_data["Landmark"] != "LM17"]

core_data = outline_data[outline_data["Landmark"] != "LM13"]
core_data = core_data[core_data["Landmark"] != "LM12"]
core_data = core_data[core_data["Landmark"] != "LM11"]
core_data = core_data[core_data["Landmark"] != "LM10"]
core_data = core_data[core_data["Landmark"] != "LM9"]
core_data = core_data[core_data["Landmark"] != "LM8"]
core_data = core_data[core_data["Landmark"] != "LM7"]
core_data = core_data[core_data["Landmark"] != "LM6"]

# DBSCAN

# First sample

raw_data = pd.DataFrame.to_numpy(all_data[["x","y","z"]], dtype = "float64")

db = DBSCAN(eps = 0.04, min_samples = 30).fit(raw_data)

core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
  if k == -1:
  # Black used for noise.
    col = [0, 0, 0, 1]
  class_member_mask = (labels == k)
  xy = raw_data[class_member_mask & core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=5)
  xy = raw_data[class_member_mask & ~core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=1)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()
stats.find_repeats(labels)

# Second sample

raw_data = pd.DataFrame.to_numpy(outline_data[["x","y","z"]], dtype = "float64")

db = DBSCAN(eps = 0.04, min_samples = 30).fit(raw_data)

core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
  if k == -1:
  # Black used for noise.
    col = [0, 0, 0, 1]
  class_member_mask = (labels == k)
  xy = raw_data[class_member_mask & core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=5)
  xy = raw_data[class_member_mask & ~core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=1)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()
stats.find_repeats(labels)

# Third sample

raw_data = pd.DataFrame.to_numpy(core_data[["x","y","z"]], dtype = "float64")

db = DBSCAN(eps = 0.04, min_samples = 30).fit(raw_data)

core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
unique_labels = set(labels)
colors = [plt.cm.Spectral(each)
          for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
  if k == -1:
  # Black used for noise.
    col = [0, 0, 0, 1]
  class_member_mask = (labels == k)
  xy = raw_data[class_member_mask & core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=5)
  xy = raw_data[class_member_mask & ~core_samples_mask]
  plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
           markeredgecolor='k', markersize=1)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()
stats.find_repeats(labels)

# Mean Shift

# First Sample

raw_data = pd.DataFrame.to_numpy(all_data[["x","y","z"]], dtype = "float64")

ms = MeanShift(cluster_all = False, bandwidth = estimate_bandwidth(raw_data, quantile = .05, n_samples = raw_data.shape[0]))
y_pred = ms.fit_predict(raw_data)

labels = ms.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
plt.scatter(raw_data[:, 0], raw_data[:, 1], c = y_pred)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()

# Second Sample

raw_data = pd.DataFrame.to_numpy(outline_data[["x","y","z"]], dtype = "float64")

ms = MeanShift(cluster_all = False, bandwidth = estimate_bandwidth(raw_data, quantile = .05, n_samples = raw_data.shape[0]))
y_pred = ms.fit_predict(raw_data)

labels = ms.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
plt.scatter(raw_data[:, 0], raw_data[:, 1], c = y_pred)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()

# Third Sample

raw_data = pd.DataFrame.to_numpy(core_data[["x","y","z"]], dtype = "float64")

ms = MeanShift(cluster_all = False, bandwidth = estimate_bandwidth(raw_data, quantile = .1, n_samples = raw_data.shape[0]))
y_pred = ms.fit_predict(raw_data)

labels = ms.labels_
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)
plt.scatter(raw_data[:, 0], raw_data[:, 1], c = y_pred)
plt.title('Estimated number of clusters ({}) and noise points ({})'.format(n_clusters_, n_noise_))
plt.axis("off")
plt.show()

# Data you wish to plot in r can be extracted via:

np.savetxt("data_to_plot.txt",
           np.stack((raw_data[:,0],raw_data[:,1],labels), axis = 1))
