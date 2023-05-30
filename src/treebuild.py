import treelib as tl

def build_genetic_trees(freqs):
    #a tree list is created
    trees = []
    #the initial tree is created with a frequency of 1
    arbol = tl.Tree()
    arbol.create_node(tag="0",identifier="0",data={"frecuencia": 1.0})

    #the initial tree is added to the tree list
    trees.append(arbol)

    #for each frequency in the frequency list the following is done
    while True:
        #the most frequent node is taken
        max_freq_idx = freqs.index(max(freqs))+1
        #max_freq_idx = freqs[np.where(freqs == max(freqs))][0]
        max_freq = freqs[max_freq_idx-1]
        #is removed from the frequency list
        freqs[max_freq_idx-1] = 0.0 
        if max_freq == 0.0:
            break
        #an auxiliary L1 list is created that will later replace the list of trees
        L1 = []
        #for each tree in the tree list the following is done
        for tree in trees:
            #print(tree)
            #for each node in the tree the following is done
            for node in tree.all_nodes():
                #if the node has a frequency greater than zero
                if node.data["frecuencia"] >= max_freq:
                    #print("node freq", node.data["frecuencia"], "max_freq", max_freq)
                    #a tree copy is created
                    tree_copy = tl.Tree(tree.subtree(tree.root), deep=True)
                    
                    #a node is created with the frequency of the most frequent node taken before
                    tree_copy.create_node(tag=str(max_freq_idx)+"_"+str(max_freq),identifier=str(max_freq_idx)+"_"+str(max_freq),data={"frecuencia": max_freq},parent=node.identifier)
                    #print("tree_copy", tree_copy)
                    #print("tree", tree)
                    #the frequency of the node is updated
                    new_node = tree_copy.get_node(node.identifier)
                    new_node.data["frecuencia"] = new_node.data["frecuencia"] - max_freq
                    #print("new_node", new_node)
                    #the tree is added to the auxiliary list
                    L1.append(tree_copy)
        #the tree list is replaced by the auxiliary list
        #for tree in L1:
        #    print(tree)
        #    for node in tree.all_nodes():
                #print(node)
        #print("***********************")
        trees = L1
        
    return trees

def modify_tree_identifier(tree,mutations_list):
    for node in tree.all_nodes():
        node.update_node( (str(mutations_list[int(node.identifier.split("_")[0])])+"_"+str(node.identifier)) )

#builds a genetic tree from a cluster frequency list
#so with the tree built it modifies the name of the nodes to store the mutations that have occurred in the cluster
def build_cluster_tree(cluster_dict):
    freq_list = []
    mutations_list = []
    for key in cluster_dict.keys():
        freq_list.append(float(key))
        mutations_list.append(cluster_dict[key])
    
    mutations_list = [','.join(map(str, sublist)) for sublist in mutations_list]


    trees = build_genetic_trees(freq_list)
    #print("***********************aaaaaaaaaaaaaaaaa")
    #print(mutations_list)
    #print("**********************aaaaaaaaaaaaaaaaaaa")
    for tree in trees:
        #if remove_root_0:
            #tree.show()
            #tree.subtree(tree.get_node(tree.root).fpointer[0])
        for node in tree.all_nodes():
            #print(node.identifier.split("_")[0])
            if int(node.identifier.split("_")[0])-1 >= 0:
                tree.update_node(node.identifier, identifier=(str(mutations_list[int(node.identifier.split("_")[0])-1 ])+"_"+str(node.identifier)),
                             tag =((mutations_list[int(node.identifier.split("_")[0])-1 ])+"_"+str(node.identifier)) )
    return trees


#for a dictionary of clusters it builds a dictionary of lists of trees
def build_cluster_trees_dict(inferred):
    inferred_alg1 = inferred.get("Algorithm_1")
    inferred_alg2 = inferred.get("Algorithm_2_kmeans")
    inferred_alg2_dict = {}
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


    clusters_trees_dict = {}
    for key in inferred_dict.keys():
        clusters_trees_dict[key] = build_cluster_tree(inferred_dict[key])
    return clusters_trees_dict