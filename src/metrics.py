from sklearn.metrics.cluster import v_measure_score,adjusted_rand_score, pair_confusion_matrix
import tools, treebuild

#count the number of clusters in the inferred minus the number of clusters in the ground truth
def clone_number_error(inferred,ground_truth):
    return len(inferred) - len(ground_truth)


def v_measure_metric(inferred,ground_truth,total):
    #v measure
    #print(ground_truth)
    ground_truth_labels = label_list(ground_truth,total)
    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    #print(inferred)
    inferred_labels = label_list(inferred,total)
    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    #print(inferred_labels)
    v_measure = v_measure_score(ground_truth_labels, inferred_labels)
    #print("v_measure: "+str(v_measure))
    return v_measure

def adjusted_rand_index_metric(inferred,ground_truth,total):
    #adjusted rand index
    ground_truth_labels = label_list(ground_truth,total)
    inferred_labels = label_list(inferred,total)
    adjusted_rand_index = adjusted_rand_score(ground_truth_labels, inferred_labels)
    #print("adjusted_rand_index: "+str(adjusted_rand_index))
    return adjusted_rand_index

def pair_confusion_matrix_metric(inferred,ground_truth,total):
    #adjusted rand index
    ground_truth_labels = label_list(ground_truth,total)
    inferred_labels = label_list(inferred,total)
    pair_confusion = pair_confusion_matrix(ground_truth_labels, inferred_labels)
    #print("adjusted_rand_index: "+str(adjusted_rand_index))
    return pair_confusion


def label_list(dictionary,size):
    label_list = ["0"]*size
    for k, v in dictionary.items():
        for i in v:
            label_list[i-1] = k
    return label_list




"""
from stereodist import CASet, DISC
#calculate the Disc and Caset metrics
def tree_distance(tree,g_newick):
    trees = [tree,g_newick]
    print(trees)
    #print(g_newick)
    distance_measure = CASet.caset_intersection
    distance_matrix = [[None for x in range(len(trees))] for y in range(len(trees))]

    for i, t1 in enumerate(trees):
        for j, t2 in enumerate(trees):
            distance_matrix[i][j] = distance_measure(t1, t2)
            distance_matrix[j][i] = distance_matrix[i][j]

    output = CASet.get_matrix_string(distance_matrix)
    print(output)
    return output
"""



def calculate_metrics(inferred,ground_truth,ground_truth_tree = None,remove_root_0=True):
    #get the total number of mutations in the ground truth
    total = ground_truth["total"]
    #an adapted ground truth is created to be able to compare the results of the algorithms
    adapted_ground_truth = ground_truth.copy()
    for key in adapted_ground_truth["cluster_list"].keys():
        print("***************")
        print(adapted_ground_truth["cluster_list"][key])
        adapted_ground_truth["cluster_list"][key] = adapted_ground_truth["cluster_list"][key][1]
        print(adapted_ground_truth["cluster_list"][key])
        print("***************")
    inferred_alg1 = inferred.get("Algorithm_1")
    inferred_alg2 = inferred.get("Algorithm_2_kmeans")
    if inferred_alg2 is not None:
        inferred_alg2_dict = {}
        for i in range(len(inferred_alg2)):
            inferred_alg2_dict.setdefault("Algorithm_2_kmeans_"+(str(i+1)), inferred_alg2[str(i+1)]["centroids"])
            #inferred_alg2_list[i] = cluster_dict(inferred_alg2[i])
        #it takes the elements of inferred_alg2_dict (which are dictionaries of lists) and transform them into one dictionary of lists
        #print("*******************************")
        #print(inferred_alg2_dict)
        #print("*******************************")
        for key in inferred_alg2_dict.keys():
            d = {}
            for dict in inferred_alg2_dict[key]:
                #for element in dict.keys():
                #    
                d.update(dict)
            inferred_alg2_dict[key] = d
        #print("*******************************")
        #print(inferred_alg2_dict)
        #print("*******************************")
    
    inferred_alg3 = inferred.get("Algorithm_3")
    inferred_dict = {"Algorithm_1":inferred_alg1,"Algorithm_3":inferred_alg3}
    inferred_dict.update(inferred_alg2_dict)

    print(inferred_dict)

    #creates the newick tree of the ground truth
    #if ground_truth_tree:
    #    g_newick = tools.tree_to_newick(ground_truth_tree, ground_truth_tree.get_node(ground_truth_tree.root))

    #create dict to store benchmark results
    benchmark = {}
    for key in inferred_dict.keys():
        #the key is added to the benchmark dict
        benchmark.setdefault(key,{})
        #clone number error
        #if "Algorithm_2_kmeans" in key:
            #print("kmeans")
            #in kmeans the clustering includes the 0 frequency cluster, so it is used also in the ground truth
            #benchmark[key].setdefault("clone_number_error",clone_number_error(inferred_dict[key],ground_truth["frequency_clusters"]))
            #benchmark[key].setdefault("clone_number_error_different_frequencies",clone_number_error(inferred_dict[key],ground_truth_no_0["frequency_clusters"]))
        #else:
        benchmark[key].setdefault("clone_number_error",clone_number_error(inferred_dict[key],ground_truth["cluster_list"]))
        #clone number error considering the real frequencies (because may be 2 real clusters with the same frequency but different mutations, and this is impossible to distinghuish without more information)
        benchmark[key].setdefault("clone_number_error_different_frequencies",clone_number_error(inferred_dict[key],ground_truth["frequency_clusters"]))
        #genotype error
        print("genotype error")
        #v_measure
        #if "Algorithm_2_kmeans" in key:
            #in kmeans the clustering includes the 0 frequency cluster, so it is used also in the ground truth
        #    benchmark[key].setdefault("genotype_error",v_measure_metric(inferred_dict[key],adapted_ground_truth["cluster_list"],total))
        #    benchmark[key].setdefault("genotype_error_different_frequencies",v_measure_metric(inferred_dict[key],adapted_ground_truth["frequency_clusters"],total))
        #else:
        benchmark[key].setdefault("genotype_error",v_measure_metric(inferred_dict[key],ground_truth["cluster_list"],total))
        benchmark[key].setdefault("genotype_error_different_frequencies",v_measure_metric(inferred_dict[key],ground_truth["frequency_clusters"],total))
        #adjusted rand index
        #if "Algorithm_2_kmeans" in key:
            #in kmeans the clustering includes the 0 frequency cluster, so it is used also in the ground truth
        #    benchmark[key].setdefault("adjusted_rand_index",adjusted_rand_index_metric(inferred_dict[key],adapted_ground_truth["cluster_list"],total))
        #    benchmark[key].setdefault("adjusted_rand_index_different_frequencies",adjusted_rand_index_metric(inferred_dict[key],adapted_ground_truth["frequency_clusters"],total))
        #else:
        benchmark[key].setdefault("adjusted_rand_index",adjusted_rand_index_metric(inferred_dict[key],ground_truth["cluster_list"],total))
        benchmark[key].setdefault("adjusted_rand_index_different_frequencies",adjusted_rand_index_metric(inferred_dict[key],ground_truth["frequency_clusters"],total))
        #confusion matrix
        #if "Algorithm_2_kmeans" in key:
            #in kmeans the clustering includes the 0 frequency cluster, so it is used also in the ground truth
        #    benchmark[key].setdefault("pair_confusion_matrix",pair_confusion_matrix_metric(inferred_dict[key],adapted_ground_truth["cluster_list"],total))
        #    benchmark[key].setdefault("pair_confusion_matrix_different_frequencies",pair_confusion_matrix_metric(inferred_dict[key],adapted_ground_truth["frequency_clusters"],total))
        #else:
        benchmark[key].setdefault("pair_confusion_matrix",pair_confusion_matrix_metric(inferred_dict[key],ground_truth["cluster_list"],total))
        benchmark[key].setdefault("pair_confusion_matrix_different_frequencies",pair_confusion_matrix_metric(inferred_dict[key],ground_truth["frequency_clusters"],total))
    return benchmark


"""
        #disc and caset
        #it needs to generate the tree and pass it to newick format
        if ground_truth_tree:
            trees = treebuild.build_cluster_tree(inferred_dict[key])
            treelist = []
            for tree in trees:
                tree.show()
                tr = tools.tree_to_newick(tree, tree.get_node(tree.root))
                if remove_root_0:
                    if tr[-1] == "0":
                        tr = tr[1:-2]
                tr += "\n" + tr  + ";"
                treelist.append(tr)
            for tree in treelist:
                print(tree)
                #print(g_newick)
                dist_list = tree_distance(tree,g_newick)
""" 
        #print(tree)
        #print(inferred_dict[key])
        #print(tools.create_tree(inferred_dict[key]))
    

    












