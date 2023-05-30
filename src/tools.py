import json
import numpy as np
import os
import treelib as tl
#import nni
import pandas as pd

def remove_same_seed(df):
    #removes the rows with the same seed and same parameters
    df = df.drop_duplicates(subset=['nodes', 'total_number_of_mutations','mutation_ration','initial_mutation_ratio','base','total_number_of_cells',
                                    'seed','frequency_list','algorithm','k_means_number_of_clusters','KDE_bandwidth','KDE_search_hyperparameters',
                                    'clusters_ground_truth'],
                                    keep='first')
    return df

def join_csvs(csvs):
    #for each csv in the list, it reads the simulation column and adds it to the next csv
    df = pd.read_csv(csvs[0], sep=',', header=0)
    for i in range(1,len(csvs)):
        df2 = pd.read_csv(csvs[i], sep=',', header=0)
        #it takes the maximum value of the simulation column in the previous csv and adds it to the next csv
        df2['simulation'] = df2['simulation'] + df['simulation'].max()
        df = df.append(df2, ignore_index=True)
    return df
    

def csv_to_pandas(csv_file):
    df = pd.read_csv(csv_file, sep=',', header=0)
    return df
"""
def nni_parameterss(parameters):
    print("updating params with nni")
    updated_params = nni.get_next_parameter()
    parameters.update(updated_params)
    print("Config after update with NNI config:")
    print(json.dumps(parameters, indent=2))
    return parameters
"""
    
#adapt the output of generate_tree stored in a dictionary to the format of the tree to be able to store it in a json file
def adapt_output(datain):
    data = datain.copy()
    if type(data) == dict:
        for key in data:
            #if parameter is a tree, convert it to a dictionary
            if type(data[key]) == tl.Tree:
                data[key] = data[key].to_dict()
            #if parameter is a numpy array, convert it to a list
            elif type(data[key]) == np.ndarray:
                data[key] = data[key].tolist()
            #if parameter is another dictionary, call the function again
            elif type(data[key]) == dict:
                data[key] = adapt_output(data[key])
            #if parameter is a list, call the function again
            elif type(data[key]) == list:
                data[key] = adapt_output(data[key])
    
    #this part is not neccessary because the output of generate_tree is always a dictionary
    #but it is left here in case it is needed in the future
    elif type(data) == list:
        for i in range(len(data)):
            #if parameter is a tree, convert it to a dictionary
            if type(data[i]) == tl.Tree:
                data[i] = data[i].to_dict()
            #if parameter is a numpy array, convert it to a list
            elif type(data[i]) == np.ndarray:
                data[i] = data[i].tolist()
            #if parameter is another dictionary, call the function again
            elif type(data[i]) == dict:
                data[i] = adapt_output(data[i])
            #if parameter is a list, call the function again
            elif type(data[i]) == list:
                data[i] = adapt_output(data[i])
    return data

#reads the parameters stored at a json file
def read_json(json_file):
    with open(json_file) as f:
        data = json.load(f)
    return data

def to_tuple(data):
    for key in data:
        if isinstance(data[key], list):
            data[key] = tuple(data[key])
    return data

#check for each parameter if it is a tuple and if it is, it uses np.random.randint
# to generate a random number between the two values of the tuple
def check_tuple(data):
    for key in data:
        if isinstance(data[key], tuple):
            #if both values are ints it generates a random int
            if isinstance(data[key][0], int) and isinstance(data[key][1], int):
                data[key] = np.random.randint(data[key][0], data[key][1])
            #else it generates a random float
            else:
                data[key] = np.random.uniform(data[key][0], data[key][1])
    return data

#check if the parameter set_seed is True and if it is, it sets the seed
def set_seed(data):
    if data["set_seed"]:
        np.random.seed(data["seed"])
    print("The seed is: {}".format(np.random.get_state()[1][0]))
    return data

#stores figures into a directory
def store_fig(fig, name):
    # if the directory doesn't exist it creates it
    dir = os.path.dirname(name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    fig.savefig(name)

#stores figures into a directory
def store_fig_ax(fig, name):
    # if the directory doesn't exist it creates it
    dir = os.path.dirname(name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    fig.figure.savefig(name)

#stores a dictionary into a json file
def store_dict(data, name):
    # if the directory doesn't exist it creates it
    dir = os.path.dirname(name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    #if the file doesn't exist it creates it
    with open(name, 'w') as fp:
        json.dump(data, fp,indent=4)

def store_string(data, name):
    # if the directory doesn't exist it creates it
    dir = os.path.dirname(name)
    if not os.path.exists(dir):
        os.makedirs(dir)
    #if the file doesn't exist it creates it
    with open(name, 'w') as fp:
        fp.write(data)


# Función recursiva para agregar nodos y ramas al árbol
def add_node(tree, parent, node_dict):
    for element in node_dict["children"]:
        #if the element is a string it means that it is a leaf
        if isinstance(element, str):
            tree.create_node(element, element, parent=parent)
        #if the element is a dictionary it means that it is a node
        elif isinstance(element, dict):
            name = list(element.keys())[0]
            tree.create_node(name, name, parent=parent)
            new_node = tree.get_node(name)
            add_node(tree, new_node, element[name])



def create_tree(data):
    tree = tl.Tree()
    #take the first element of the dictionary as the root (it should be only one)
    name = list(data.keys())[0]
    tree.create_node(name, name)
    #root = Node("main_tree")
    #tree.add_node(root)
    add_node(tree, tree.get_node(name), data[name])
    return tree

    

#it uses the trees transformed into strings to compare them
def compare_trees(t1,t2):
    print(t1.to_json())
    print(t2.to_json())
    return t1.to_json() == t2.to_json()



#generate a cluster list with the mutations that have occurred in each cluster
def cluster_treelist(tree):
    cluster = dict()
    #take the root of the tree
    root = tree.get_node(tree.root)
    number_mutations = len(root.tag.split("_")[0])
    #mutations_added = np.zeros(number_mutations)
    mutations_added = ["0"] * number_mutations
    total = 0
    #for each node in the tree
    for node in tree.all_nodes():
        #gets the number of cells in the cluster
        #print(node.identifier)
        total+= int(node.tag.split("_")[1])
    
    mutations_added = get_cluster_mutations(root,total,cluster,mutations_added,tree)
    #print(cluster)
    return cluster

def get_cluster_mutations(node,total,cluster,mutations_added,tree):
    name,number = node.tag.split("_")
    mutations = [name[i] if name[i] != mutations_added[i] else '0' for i in range(len(name))]
    #print(len(mutations_added))
    #print(len(name))
    #returns the positions of the mutations that have occurred in the cluster
    mutations_positions = [i for i, x in enumerate(mutations) if x == "1"]
    #cluster.update( {name: {(int(number)/total) :mutations_positions} })
    cluster.update( {name: [(int(number)/total) ,mutations_positions] })
    mutations_added = [ max(x, y) for (x,y) in zip(name,mutations_added) ]
    #for each child of the node
    for child in node.fpointer:
        mutations_added = get_cluster_mutations(tree.get_node(child),total,cluster,mutations_added,tree)

    return mutations_added

#makes a cluster dictionary from a frequency list
def cluster_dict(frequency_list):
    cluster = {}
    for i in range(len(frequency_list)):
        cluster.setdefault(frequency_list[i], []).append(i)
    return cluster


#pass tree to newick format
def tree_to_newick(tree,node):
    newick = ""
    if len(node.identifier.split("_")[0].split(",")) > 1:
        name = "{" + ",".join(node.identifier.split("_")[0].split(",")) + "}"
    else:
        name = node.identifier.split("_")[0]
    if node.is_leaf():
        newick += name
    else:
        newick += "(" + ",".join([tree_to_newick(tree,tree.get_node(n)) for n in node.fpointer]) + ")" + name
    return newick

"""

#pass tree to newick format
def tree_to_newick(tree):
    newick = ""
    for node in tree.all_nodes():
        if node.is_leaf():
            newick += node.identifier
        else:
            newick += "(" + ",".join([tree.get_node(n).identifier for n in node.fpointer]) + ")"
    return newick + ";"

#pass tree to newick format
def tree_to_newick(tree):
    newick = ""
    for node in tree.all_nodes():
        if len(node.identifier.split("_")[0].split(",")) > 1:
            name = "{" + ",".join(node.identifier.split("_")[0].split(",")) + "}"
        else:
            name = node.identifier.split("_")[0]
        if node.is_leaf():
            newick += name
        else:
            newick += "(" + ",".join([tree.get_node(n).identifier for n in node.fpointer]) + ")"
    return newick + ";"

def tree_to_newick(tree):
    newick = ""
    root = tree.root
    newick += newick_node(tree,tree.get_node(root)) + ";"
    return newick
"""

def newick_node(tree,node):
    if node.is_leaf():
        return node.identifier
    else:
        newick = "("
        for child in node.fpointer:
            newick += newick_node(tree,tree.get_node(child)) + ","
        return newick[:-1] + ")" + node.identifier
    
def clusters_to_newick(clusters,remove_root_0):
    out = ""
    for key in clusters:
        out += "\n" + key + ":" 
        for element in clusters[key]:
            t = create_tree(element)
            tr = tree_to_newick(t, t.get_node(t.root))
            if remove_root_0:
                if tr[-1] == "0":
                    tr = tr[1:-2]
            out += "\n" + tr  + ";"
        out += "\n"
    return out


"""

#pass tree to newick format
def tree_to_newick(tree,node,remove_root_0):

    newick = ""
    if len(node.identifier.split("_")[0].split(",")) > 1:
        name = "{" + ",".join(node.identifier.split("_")[0].split(",")) + "}"
    else:
        name = node.identifier.split("_")[0]
    if node.is_leaf():
        newick += name
    else:
        newick += "(" + ",".join([tree_to_newick(tree,tree.get_node(n),remove_root_0) for n in node.fpointer]) + ")" + name
    return newick


#pass tree to newick format
def tree_to_newick(tree):
    newick = ""
    for node in tree.all_nodes():
        if node.is_leaf():
            newick += node.identifier
        else:
            newick += "(" + ",".join([tree.get_node(n).identifier for n in node.fpointer]) + ")"
    return newick + ";"

#pass tree to newick format
def tree_to_newick(tree):
    newick = ""
    for node in tree.all_nodes():
        if len(node.identifier.split("_")[0].split(",")) > 1:
            name = "{" + ",".join(node.identifier.split("_")[0].split(",")) + "}"
        else:
            name = node.identifier.split("_")[0]
        if node.is_leaf():
            newick += name
        else:
            newick += "(" + ",".join([tree.get_node(n).identifier for n in node.fpointer]) + ")"
    return newick + ";"

def tree_to_newick(tree):
    newick = ""
    root = tree.root
    newick += newick_node(tree,tree.get_node(root)) + ";"
    return newick

def newick_node(tree,node):
    if node.is_leaf():
        return node.identifier
    else:
        newick = "("
        for child in node.fpointer:
            newick += newick_node(tree,tree.get_node(child)) + ","
        return newick[:-1] + ")" + node.identifier
    
def clusters_to_newick(clusters,remove_root_0):
    out = ""
    for key in clusters:
        out += "\n" + key + ":" 
        for element in clusters[key]:
            t = create_tree(element)
            out += "\n" + tree_to_newick(t, t.get_node(t.root),remove_root_0) + ";"
        out += "\n"
    return out

"""


"""


tree_to_newick(t)


    newick = newick + "," + root.identifier + ")"
    for child in node.fpointer:
        mutations_added = get_cluster_mutations(tree.get_node(child),total,cluster,mutations_added,tree)
    #newick += ",".join([tree.get_node(n).identifier for n in root.fpointer]) + ")"



#create a tree from a dictionary
def create_tree(data):
    tree = tl.Tree()
    tree.create_node(data["name"], data["name"], data=data)
    for child in data["children"]:
        tree.create_node(child["name"], child["name"], parent=data["name"], data=child)
        create_tree(child)
    return tree


#create a tree from a dictionary
def create_tree(data):
    added = set()
    tree = Tree()
    while data:

        for key, value in data.items():
            if value['parent'] in added:
                tree.create_node(key, key, parent=value['parent'])
                added.add(key)
                data.pop(key)
                break
            elif value['parent'] is None:
                tree.create_node(key, key)
                added.add(key)
                data.pop(key)
                break



import json
from treelib import Node, Tree

# Función recursiva para agregar nodos y ramas al árbol
def add_node(tree, parent, node_dict):
    for key, value in node_dict.items():
        if key != "children":
            # Si la clave no es "children", agregamos un nodo
            node = Node(key)
            tree.add_node(node, parent)
            if isinstance(value, dict):
                # Si el valor es un diccionario, llamamos a la función recursivamente
                add_node(tree, node.identifier, value)
            else:
                # Si el valor no es un diccionario, lo agregamos como una hoja del árbol
                leaf = Node(f"{key}: {value}")
                tree.add_node(leaf, node.identifier)
        else:
            # Si la clave es "children", agregamos los nodos hijos
            for child in value:
                if isinstance(child, dict):
                    add_node(tree, parent, child)
                else:
                    leaf = Node(child)
                    tree.add_node(leaf, parent)


def compare_trees(tree1, tree2): 
    # Comprobar si las raíces son iguales
    if tree1.root != tree2.root:
        return False
    
    # Comprobar si los hijos de la raíz son iguales
    if tree1.children(tree1.root) != tree2.children(tree2.root):
        return False
    
    # Recorrer los hijos de la raíz y comprobar si son iguales
    for child in tree1.children(tree1.root):
        if not are_trees_equal(tree1.subtree(child.identifier), tree2.subtree(child.identifier)):
            return False
    
    # Si todo ha ido bien, los árboles son iguales
    return True
"""



