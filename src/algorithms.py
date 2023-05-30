"""
This file contains the algorithms used to generate the clusters
"""
#for the run_algorithm function
import tools
import os

#needed in kmeans and KDE
import matplotlib.pyplot as plt

#needed in kmeans
from sklearn.cluster import KMeans
import numpy as np


#for k-means warnings
#remove warnings
import warnings
warnings.filterwarnings("ignore")

#needed in KDE
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from scipy.signal import argrelextrema


########################################################################
########################################################################
########################################################################
"""
Algorithm 1: Decide a number of clusters from the total number of cells
"""
def make_clusters(freq_list,total):
    """
    total: total number of cells
    combined_list: list of the frequency of each mutation in each cluster

    Takes the total number of cells and generates a cluster with the frequency assigned to according to the total number of cells
    So for example if we have 5 cells, 5 clusters can be generated with a frequency of 0.2, 0.4, 0.6, 0.8 and 1.0, since those are all the possible frequencies that can be generated with 5 cells (without counting the noise) (although that does not mean that all these clusters will be generated)
    """
    clusters = {}
    for i in range(total):
        #a dictionary is generated with an element that will be i/total
        clusters.update({(i+1)/total:[]})
    
    count = 0
    for i in freq_list:
        #each mutation is added to the cluster with the closest frequency
        clusters[min(clusters.keys(), key=lambda k: abs(k - i))].append(count)
        count += 1
    
    #empty clusters are removed
    for i in list(clusters):
        if clusters[i] == []:
            del clusters[i]

    return clusters



########################################################################
########################################################################
########################################################################
"""
Algorithm 2: Using K-means
"""

#elbow method can be used to choose the number of clusters

from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt

def elbow_method(frecuencias,total):

    #frequency_list of mutations
    l = frecuencias
    if total > len(l):
        total = len(l)

    #transform the list into a numpy array
    X = np.array(l).reshape(-1, 1)

    #define an empty list to store the values of the sums of the squares within each cluster
    sum_of_squares = []

    #define the maximum number of clusters to try
    max_clusters = total

    # Try different values of k (number of clusters) and store the values of the sums of the squares within each cluster
    for k in range(1, max_clusters):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(X)
        sum_of_squares.append(kmeans.inertia_)

    # Plot an elbow graph to determine the optimal number of clusters
    fig, ax = plt.subplots()
    ax.plot(range(1, max_clusters), sum_of_squares, 'bx-')
    ax.set_xlabel('Number of clusters (k)')
    ax.set_ylabel('Sum of squared distances to the cluster center')
    ax.set_title('Elbow method to determine the optimal K')
    #plt.show()
    return fig


def cluster_kmeans(frequencies, max_num_clusters):
    """
        Generates clusters with kmeans trying from 1 to max_num_clusters clusters
    """
    if max_num_clusters > len(frequencies):
        max_num_clusters = len(frequencies)
    #transform the list into a numpy array
    X = np.array(frequencies).reshape(-1, 1)
    #define an empty list to store the values of the sums of the squares within each cluster
    sum_of_squares = []
    #define an inertia value of -1 to start with an impossible value of inertia
    inertia = -1
    #creates a dictionary to store the clusters
    clusters = dict()
    # Try different values of k (number of clusters) and store the values of the sums of the squares within each cluster
    for k in range(1, max_num_clusters):
        kmeans = KMeans(n_clusters=k)
        kmeans.fit(X)
        labels = kmeans.labels_
        centroids = kmeans.cluster_centers_
        print("k:", k)
        print("Etiquetas de clusters:", labels)
        print("Centroides:", centroids)
        sum_of_squares.append(kmeans.inertia_)
        #if the sum of the squares within the cluster does not improve, it stops, because it probably can't be improved anymore
        if inertia != -1 and kmeans.inertia_ >= inertia:
            break

        #add the clusters to the list

        #makes a list of dictionaries with the clusters and its labels
        l = []
        for i in range(len(centroids)):
            l.append({str(centroids[i][0]):[index for index, element in enumerate(labels) if element == i]})
        clusters.update({str(k):{"k":k, "centroids":l, "inertia":kmeans.inertia_}})
        #clusters.update({str(k):{"k":k,"labels":labels, "centroids":centroids, "inertia":kmeans.inertia_}})

        inertia = kmeans.inertia_
    return clusters


########################################################################
########################################################################
########################################################################
"""
Algorithm 3: Detecting the peaks of the histogram using KDE
"""


def kde_cluster(freq_list,bandwidth=0.02,search_hyperparameters=False):

    #f to ndarray
    freq_list = np.array(freq_list)
    f_d = np.linspace(0, 1, 2000)

    ###############################################
    #cross validation to find the best hyperparameters for the kernel and the bandwidth
    if search_hyperparameters:
        param_grid = {'kernel': ['gaussian', 'epanechnikov', 'exponential', 'linear'],
              'bandwidth' : np.linspace(0.005, 1, 10)
             }

        grid = GridSearchCV(
                estimator  = KernelDensity(),
                param_grid = param_grid,
                n_jobs     = -1,
                cv         = 10, 
                verbose    = 0
            )

        _ = grid.fit(X = freq_list.reshape((-1,1)))

        print(grid.best_params_, ":", grid.best_score_, grid.scoring)

        kde = grid.best_estimator_

    ###############################################

    else:
        kde = KernelDensity(bandwidth=bandwidth, kernel="gaussian")
        kde.fit(freq_list[:, None])


    # score_samples returns the log of the probability density
    logprob = kde.score_samples(f_d[:, None])
    fig, ax = plt.subplots()
    ax.fill_between(f_d, np.exp(logprob), alpha=0.5)
    ax.plot(freq_list, np.full_like(freq_list, -0.01), '|k', markeredgewidth=1)
    #plt.ylim(-0.02, 0.22)
    #plt.show()



    mi, ma = argrelextrema(logprob, np.less)[0], argrelextrema(logprob, np.greater)[0]

    print("Minima:", f_d[mi])
    print("Maxima:", f_d[ma])
    print("********************************")
    #print(f[f < f_d[mi][0]], f[(f >= f_d[mi][0]) * (f <= f_d[mi][1])], f[f >= f_d[mi][1]])
    #print(f[f < f_d[mi][0]], f[(f >= f_d[mi][0]) * (f < f_d[mi][1])], f[f >= f_d[mi][1]])
    print(freq_list)


    f_d[mi]
    #a dictionary that stores clusters of mutations with the same frequency, where the key is the maximum frequency of the cluster
    clusters = {}
    for i in f_d[mi]:
        clusters.update({i:[]})
    clusters.update({1:[]})

    #the frequencies are added to the clusters, skipping the ones with frequency 0
    for i in range(len(freq_list)):
        if freq_list[i] != 0:
            for j in clusters.keys():
                if freq_list[i] <= j:
                    clusters[j].append(i)
                    break

    return clusters, fig


########################################################################
########################################################################
########################################################################


def run_algorithm(PARAMS,algorithm_number):

        #cluster will store the results of the cluster algorithms in a dictionary to be stored later in a json file
        cluster = dict()
        ALG_PARAMS = PARAMS["algorithm"]
        CELL_DATA = tools.read_json(ALG_PARAMS["data_file"])

        if 1 in algorithm_number:
            print("Running the first algorithm...")
            #cluster the mutations
            clusters = make_clusters(CELL_DATA[ALG_PARAMS["frequency_list_name"]],CELL_DATA["total"])
            #store the clusters in the dictionary
            cluster.update({"Algorithm_1":clusters})
        
        if 2 in algorithm_number:
            print("Running the second algorithm(Kmeans)...")
            #cluster the mutations
            fig_name = "elbow_method.png"
            fig_path = os.path.join(PARAMS["out_dir"],fig_name)
            tools.store_fig(elbow_method(CELL_DATA[ALG_PARAMS["frequency_list_name"]],CELL_DATA["total"]),fig_path)
            clusters = cluster_kmeans(CELL_DATA[ALG_PARAMS["frequency_list_name"]],CELL_DATA["total"])
            #store the clusters in the dictionary
            cluster.update({"Algorithm_2_kmeans":clusters, "Algorithm_2_elbow_method":fig_path})
        
        if 3 in algorithm_number:
            print("Running the third algorithm...")
            fig_name = "KDE_histogram.png"
            fig_path = os.path.join(PARAMS["out_dir"],fig_name)
            alg_in = {"freq_list" : CELL_DATA[ALG_PARAMS["frequency_list_name"]]}
            if "KDE" in ALG_PARAMS:
                if "bandwidth" in ALG_PARAMS["KDE"]:
                    alg_in["bandwidth"] = ALG_PARAMS["KDE"]["bandwidth"]
                if "search_hyperparameters" in ALG_PARAMS["KDE"]:
                    alg_in["search_hyperparameters"] = ALG_PARAMS["KDE"]["search_hyperparameters"]
            #cluster the mutations
            clusters,fig = kde_cluster(**alg_in)
            tools.store_fig(fig,fig_path)
            #store the clusters in the dictionary
            cluster.update({"Algorithm_3":clusters, "Algorithm_3_KDE_histogram":fig_path})
        
        return cluster