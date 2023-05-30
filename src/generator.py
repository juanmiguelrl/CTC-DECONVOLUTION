import treelib as tl
import numpy as np
import random
import scipy.stats as stats

def generate_tree(nodes,mutations,mutation_ratio,initial_mutation_ratio,base,min_minimum_mutation_range=0,max_minimum_mutation_range = 1,min_population_multiplier_range=2,max_population_multiplier_range=5):
    """"
    Generates a tree of mutation clusters with random nodes and mutations taking as input:

    nodes: Number of nodes in the tree (the number of clusters of the genetic tree)
    mutations: Number of mutations in the tree (the number of mutations in the genetic tree)
    mutation_ratio: A mutation ratio (to indicate how  possible is a mutation to happen in the initial root mutated parent cluster)
    initial_mutation_ratio: probability of mutation in the initial cluster
    base: number of cells in the initial cluster
    minimum_mutation_range: range of minimum mutations in a cluster

    And returns:
    
    main_tree: the CTC cluster lineages tree generated
    frequency_list: the mutation frequencies
    frequency_list_noisy: the mutation frequencies with noise
    cluster_list: the list of clusters with their corresponding mutations
    total: the total number of cells in the tree 
    
    """

    #cluster_list mutations utilizadas (se usa para mantener el teorema de los sitios infinitos) (ISA)))
    cluster_list_mutations = ""
    #total number of cells
    total = 0

    def mutate(parent,cluster_list_mutations):
        """"
        Creates a cluster_list of mutations and a name for a child cluster, taking the cluster of the parent node and generating new mutations, generates the first mutation and from there there is a possibility that another mutation will be generated from the mutation_ratio parameter, and for each new mutation this possibility is reduced more
        To avoid repeating mutations, the cluster_list of mutations that have already occurred is used to select a mutation at random that has not yet occurred
        The minimum_mutation_range parameter is used to indicate the minimum number of mutations that a cluster must have
        """
        p = mutation_ratio * 2
        min = np.random.randint(min_minimum_mutation_range, max_minimum_mutation_range)
        name = parent.identifier
        while True:
            #print("cluster_list_mutations", cluster_list_mutations)
            """generates a random mutation that has not yet happened"""
            n = np.random.choice([i for i, x in enumerate(cluster_list_mutations) if x == "0"])+1
            #print("n", n)
            if len(name) >= n:
                name = name[:n-1] + "1" + name[n:]
            else:
                name = name[:n-1] + "1"
            #print("name", name)
            #cluster_list_mutations = [ max(x,y) for x in name for y in cluster_list_mutations ]
            cluster_list_mutations = [ max(x, y) for (x,y) in zip(name,cluster_list_mutations) ]
            #print("cluster_list_mutations_modificada", cluster_list_mutations)

            min = min - 1
            #print("min", min)
            if min < 0:
                p = p/2
                if not np.random.choice([False,True], p=[1-p, p]):
                    break
        #cluster_list_mutations = ["1" if x == "1" else "0" for x,y in name,cluster_list_mutations]
        #name = "".join([np.random.choice(["0","1"], p=[mutation_ratio, 1-mutation_ratio]) if x == "0" else "1" for x in parent.identifier])
        #print("name", name)
        #print("cluster_list_mutations", cluster_list_mutations)
        #print("name", name)
        return name,cluster_list_mutations

    def get_population(base,level):
        return np.random.randint(base+level*min_population_multiplier_range, base+1+level*max_population_multiplier_range)


    """" The initial tree is created with the number of nodes """
    main_tree = tl.Tree()
    name = "".join(np.random.choice([0,1], p=[1-initial_mutation_ratio, initial_mutation_ratio], size=mutations).astype(str))
    cluster_list_mutations = name

    """The mutated parent node from which the other nodes will derive is generated"""
    n = np.random.choice([i for i, x in enumerate(cluster_list_mutations) if x == "0"])+1
    #print("n", n)
    if len(name) >= n:
        name = name[:n-1] + "1" + name[n:]
    else:
        name = name[:n-1] + "1"
    #print("name", name)
    #cluster_list_mutations = [ max(x,y) for x in name for y in cluster_list_mutations ]
    """The cluster_list of mutations is a variable that will be used to keep track of which mutations of the possible ones have happened in each cluster, in this way the infinite sites assumtion (ISA) can be maintained"""
    cluster_list_mutations = name

    #the number of cells of the initial mutated parent node is generated
    number = get_population(base,0)
    total = number
    main_tree.create_node(name+ "_" + str(number), name)


    def create_mutated_child(parent,mutation_ratio,cluster_list_mutations):
        """"
        Creates a mutated child node from a parent node and adds it to the tree, returns the number of cells of the child node and the cluster_list of mutations
        """
        #en un string de 0s y 1s, cambiar un numero de mutations al azar a 1s
        #name = "".join([np.random.choice(["0","1"], p=[mutation_ratio, 1-mutation_ratio]) if x == "0" else "1" for x in parent.identifier])
        name,cluster_list_mutations = mutate(parent,cluster_list_mutations)
        number = get_population(base,main_tree.level(parent.identifier )  + 1)
        #add the child to the tree
        #this check now should not be necessary as in the function mutate a correct valid name is always returned
        if main_tree.get_node(name) is None:
            main_tree.create_node(name+ "_" + str(number), name, parent=parent.identifier)
            return number,cluster_list_mutations
        else:
            return 0,cluster_list_mutations


    #choose a random node to mutate
    """
        A parent node from which the mutated child node will derive is chosen randomly, the process is repeated as many times as indicated in the nodes variable
    """
    for i in range(nodes):
        node = np.random.choice(main_tree.all_nodes())
        #print(node.tag)
        n,cluster_list_mutations = create_mutated_child(node,mutation_ratio,cluster_list_mutations)
        total += n

    """Finally generates a cluster_list with the mutations that have occurred in each cluster"""
    #main_tree.show()
    cluster_list = ([x.identifier for x in main_tree.all_nodes()])
    #print(cluster_list)
    cluster_list_imprimir = ([x.tag.split("_")[0] for x in main_tree.all_nodes()])
    numero_celulas = ([x.tag.split("_")[1] for x in main_tree.all_nodes()])
    #make a list of tuples with the mutations and the number of cells
    mix_list = [list(x) for x in zip(cluster_list_imprimir, numero_celulas)]
    #print(combined_list)
    #print(cluster_list_imprimir)
    #print(numero_celulas)

    def get_mutations_rates(tree_list):
        """
            Calculates the frequency of each mutation using the number of cells in each cluster and the total number of cells
        """
        number_clusters = len(tree_list)
        number_mutations = len(tree_list[0][0])
        mutation_rates = np.zeros(number_mutations)
        for tree in tree_list:
            mutation_rates += (np.array([int(x) for x in tree[0]])*int(tree[1]))
            #print(mutation_rates)
        mutation_rates = mutation_rates/total#len(tree_list)
        #print(mutation_rates)
        return mutation_rates

    def apply_noise(frequencies,bound,sigma=100,mu=0):
        """
            Adds noise to the mutation frequencies of each mutation following a truncated normal distribution
        """
        noise = stats.truncnorm( (-bound - mu) / sigma, (bound - mu) / sigma, loc=mu, scale=sigma).rvs(len(frequencies))

        #return frequencies + noise
        return np.where(frequencies>0,frequencies + noise,frequencies)

    #print(get_mutations_rates(cluster_list_imprimir))
    #print("total", total)
    #freq_list = get_mutations_rates(mix_list)
    #combined_list = [list(x) for x in zip(cluster_list_imprimir, freq_list)]
    #return main_tree,freq_list,cluster_list,total,combined_list
    frequency_list = get_mutations_rates(mix_list)
    frequency_list_noisy = apply_noise(frequency_list,(1/total/2))

    #make 0 all the negative frequencies using list comprehension and 1 all the frequencies that are bigger than 1
    frequency_list_noisy = [0 if x < 0 else x for x in frequency_list_noisy]
    frequency_list_noisy = [1 if x > 1 else x for x in frequency_list_noisy]

    #remove the 0s
    frequency_list_no_zeros = [x for x in frequency_list if x != 0]
    frequency_list_noisy_no_zeros = [x for x in frequency_list_noisy if x != 0]

    return main_tree,frequency_list,frequency_list_noisy,frequency_list_no_zeros,frequency_list_noisy_no_zeros,cluster_list,total







