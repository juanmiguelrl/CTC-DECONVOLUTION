import generator, tools,algorithms, metrics
import os
import pandas as pd
import time

def run_benchmark(PARAMS):

    df = pd.DataFrame(columns=['simulation','nodes', 'total_number_of_mutations','mutation_ration','initial_mutation_ratio','base','total_number_of_cells',
                               'frequency_list','algorithm','k_means_number_of_clusters',
                               'clusters_ground_truth','seed'
                               'KDE_bandwidth','KDE_search_hyperparameters',
                               'clone_number_error' ,'v_measure','adjusted_rand_score','pair_confusion_matrix','time'])

    """ First it generates the data to test"""
    #for i in range(len(PARAMS["generate_list"]["names"])):
    simulation = 0
    for data_generator in PARAMS["generate_list"]:
        count = 0
        repeat = True
        #if generate element has the option repeat with a number, it repeats the generation of the dataset that number of times
        while repeat:
            generate = data_generator["params"].copy()
            #to store the number of simulation
            simulation += 1

            if data_generator.get("repeat"):
                if count < data_generator["repeat"]:
                    count += 1
                else:
                    repeat = False
                    continue
            #if the parameters are lists, it transforms them into tuples
            tools.to_tuple(generate)
            #if the parameters are tuples, it generates a random number between the two values
            #else it keeps the value
            generate = tools.check_tuple(generate)
            
            print("Generating the dataset...")
            resultado = dict(zip(['main_tree', 'frequency_list', 'frequency_list_noisy','frequency_list_no_zeros','frequency_list_noisy_no_zeros', 'cluster_list', 'total'], generator.generate_tree(**generate)))
            print("The dataset has been generated")
            resultado_out = tools.adapt_output(resultado)
            if PARAMS["out_file"]:
                tools.store_dict(resultado_out, os.path.join(PARAMS["out_dir"],"nodes="+str(generate["nodes"])+"_mutations="+str(generate["mutations"])+"_mutation_ratio="+str(generate["mutation_ratio"])+"_initial_mutation_ratio="+str(generate["initial_mutation_ratio"])+"_base="+str(generate["base"])+"repeat="+str(count)+".json"))

        
        #Then it generates the ground truth clusters to be used in the benchmark

            total = resultado["total"]
            for nombre_lista_frecuencias, lista_frecuencias in resultado.items():
                cluster_ground_truth = tools.cluster_treelist(resultado["main_tree"])
                total_mutations = generate["mutations"]
                if nombre_lista_frecuencias.startswith("freq"):
                    clusters_freq = tools.cluster_dict(lista_frecuencias)
                    #for list in clusters_freq.values():
                        #total += len(list)
                    if PARAMS["out_file"]:
                        dict_to_store = {"cluster_list":cluster_ground_truth,"frequency_clusters":clusters_freq,"total":total_mutations}
                        tools.store_dict(dict_to_store, os.path.join(PARAMS["out_dir"],nombre_lista_frecuencias+"_ground_truth"+"_nodes="+str(generate["nodes"])+"_mutations="+str(generate["mutations"])+"_mutation_ratio="+str(generate["mutation_ratio"])+"_initial_mutation_ratio="+str(generate["initial_mutation_ratio"])+"_base="+str(generate["base"])+"repeat="+str(count)+".json"))


                    #Then it runs the algorithms
                    cluster = dict()    

                    #alg1
                    start = time.time()
                    cluster.update({"Algorithm_1":algorithms.make_clusters(lista_frecuencias,total)})
                    end = time.time()
                    t1 = end - start
                    #alg2
                    if PARAMS["out_file"]:
                        fig_name = "elbow_method_"+nombre_lista_frecuencias+"_ground_truth"+"_nodes="+str(generate["nodes"])+"_mutations="+str(generate["mutations"])+"_mutation_ratio="+str(generate["mutation_ratio"])+"_initial_mutation_ratio="+str(generate["initial_mutation_ratio"])+"_base="+str(generate["base"])+"repeat="+str(count)+".png"
                        fig_path = os.path.join(PARAMS["out_dir"],fig_name)
                        tools.store_fig(algorithms.elbow_method(lista_frecuencias,total),fig_path)
                    start = time.time()
                    cluster.update({"Algorithm_2_kmeans":algorithms.cluster_kmeans(lista_frecuencias,total), "Algorithm_2_elbow_method":fig_path})
                    end = time.time()
                    t2 = end - start

                    alg_in = {"freq_list" : lista_frecuencias}
                    if "KDE" in PARAMS:
                        for element in PARAMS["KDE"]:
                            if "bandwidth" in element:
                                alg_in["bandwidth"] = element["bandwidth"]
                            if "search_hyperparameters" in element:
                                alg_in["search_hyperparameters"] = element["search_hyperparameters"]


                            start = time.time()
                            kde,fig = algorithms.kde_cluster(**alg_in)
                            end = time.time()
                            t3 = end - start
                            if PARAMS["out_file"]:
                                fig_name = "KDE_histogram"+"_bandwith="+str(alg_in["bandwidth"])+"_search_hyperparameters"+str(alg_in["search_hyperparameters"]) + "_"+nombre_lista_frecuencias+"_ground_truth"+"_nodes="+str(generate["nodes"])+"_mutations="+str(generate["mutations"])+"_mutation_ratio="+str(generate["mutation_ratio"])+"_initial_mutation_ratio="+str(generate["initial_mutation_ratio"])+"_base="+str(generate["base"])+"repeat="+str(count)+".png"
                                fig_path = os.path.join(PARAMS["out_dir"],fig_name)
                                tools.store_fig(fig,fig_path)
                            cluster.update({
                                "Algorithm_3"+"_bandwith="+str(alg_in["bandwidth"])+"_search-hyperparameters="+str(alg_in["search_hyperparameters"]):kde,
                                "Algorithm_3_KDE_histogram"+"_bandwith="+str(alg_in["bandwidth"])+"_search_hyperparameters"+str(alg_in["search_hyperparameters"]):fig_path})
                    
                    if PARAMS["out_file"]:
                        tools.store_dict(tools.adapt_output(cluster), os.path.join(PARAMS["out_dir"],nombre_lista_frecuencias+"_inferred_clusters"+"_nodes="+str(generate["nodes"])+"_mutations="+str(generate["mutations"])+"_mutation_ratio="+str(generate["mutation_ratio"])+"_initial_mutation_ratio="+str(generate["initial_mutation_ratio"])+"_base="+str(generate["base"])+"repeat="+str(count)+".json"))



                    #the metrics are calculated
                    inferred_alg2 = cluster.get("Algorithm_2_kmeans")
                    inferred_alg2_dict = {}
                    for i in range(len(inferred_alg2)):
                        inferred_alg2_dict.setdefault("Algorithm_2_kmeans_"+(str(i+1)), inferred_alg2[str(i+1)]["centroids"])
                    for key in inferred_alg2_dict.keys():
                        d = {}
                        for dictionar in inferred_alg2_dict[key]:
                            #for element in dict.keys():
                            #    
                            d.update(dictionar)
                        inferred_alg2_dict[key] = d
                    cluster.pop("Algorithm_2_kmeans",None)
                    cluster.update(inferred_alg2_dict)
                    #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
                    #print(cluster)
                    for key in cluster.keys():
                        extra = {}
                        if key.startswith("Algorithm_1"):
                            extra = {"name" : "Algorithm_1", "time" : t1}
                        if key.startswith("Algorithm_2"):
                            if key.startswith("Algorithm_2_elbow_method"):
                                continue
                            #Algorithm_2_kmeans_68
                            extra = {"k_means_number_of_clusters" : key.split("_")[3],"name" : "Algorithm_2_kmeans", "time" : t2}
                        if key.startswith("Algorithm_3"):
                            if key.startswith("Algorithm_3_KDE_histogram"):
                                #print("histogram")
                                continue
                                #Algorithm_3_bandwith=0.1_search_hyperparametersTrue
                                #print(key)
                            extra = {"KDE_bandwidth" : key.split("_")[2].split("=")[1],"KDE_search_hyperparameters" : key.split("_")[3].split("=")[1],"name" : "Algorithm_3_KDE", "time" : t3}
                        
                        #print(total_mutations)

                        #it modifies the cluster ground truth to adapt it to the format of the inferred clusters (removing the first element of the list, which is the frequency of the cluster)
                        #(and leaves only the mutations which compose the cluster)
                        adapted_ground_truth = cluster_ground_truth.copy()
                        for el in adapted_ground_truth.keys():
                            #print("***************")
                            adapted_ground_truth[el] = adapted_ground_truth[el][1]
                        #print(cluster[key])
                        df = df.append({'simulation':simulation,'nodes': generate["nodes"], 'total_number_of_mutations' : total_mutations,'mutation_ration':generate["mutation_ratio"]
                                        ,'initial_mutation_ratio':generate["initial_mutation_ratio"],'base':generate["base"], 'total_number_of_cells' : total,
                                    'frequency_list': nombre_lista_frecuencias,'algorithm': key,'k_means_number_of_clusters' : extra.get("k_means_number_of_clusters"),
                                    'KDE_bandwidth' : extra.get("KDE_bandwidth") ,'KDE_search_hyperparameters' : extra.get("KDE_search_hyperparameters"),

                                    'clusters_ground_truth' : 'real_different_frequencies','seed' : PARAMS["seed"],

                                    'clone_number_error' : metrics.clone_number_error(cluster[key],clusters_freq),
                                    'v_measure' : metrics.v_measure_metric(cluster[key],clusters_freq,total_mutations),
                                    'adjusted_rand_score' : metrics.adjusted_rand_index_metric(cluster[key],clusters_freq,total_mutations),
                                    'pair_confusion_matrix' : metrics.pair_confusion_matrix_metric(cluster[key],clusters_freq,total_mutations),
                                    'time' : extra.get("time")
                                    },
                                        ignore_index=True)
                        
                        df = df.append({'simulation':simulation,'nodes': generate["nodes"], 'total_number_of_mutations' : total_mutations,'mutation_ration':generate["mutation_ratio"]
                                        ,'initial_mutation_ratio':generate["initial_mutation_ratio"],'base':generate["base"], 'total_number_of_cells' : total,
                                    'frequency_list': nombre_lista_frecuencias,'algorithm': key,'k_means_number_of_clusters' : extra.get("k_means_number_of_clusters"),
                                    'KDE_bandwidth' : extra.get("KDE_bandwidth") ,'KDE_search_hyperparameters' : extra.get("KDE_search_hyperparameters"),

                                    'clusters_ground_truth' : 'real_clusters','seed' : PARAMS["seed"],

                                    'clone_number_error' : metrics.clone_number_error(cluster[key],adapted_ground_truth),
                                    'v_measure' : metrics.v_measure_metric(cluster[key],adapted_ground_truth,total_mutations),
                                    'adjusted_rand_score' : metrics.adjusted_rand_index_metric(cluster[key],adapted_ground_truth,total_mutations),
                                    'pair_confusion_matrix' : metrics.pair_confusion_matrix_metric(cluster[key],adapted_ground_truth,total_mutations),
                                    'time' : extra.get("time")
                                    },
                                        ignore_index=True)
                        
                        #store the file during each iteration so that if the program is interrupted, the results are not lost or if the results are wanted to be seen during the execution
                        #the data is directly overwritten in the csv file as this is faster than appending the data to the csv file
                        if (PARAMS["out_file"] & PARAMS["update_csv_during_run"]):
                            df.to_csv(os.path.join(PARAMS["out_dir"],'benchmark_results.csv'), index=False)
            #changes the seed for not repeating the same dataset in the next loop
            if PARAMS.get("seed"):
                PARAMS["seed"] += 1
                tools.set_seed(PARAMS)



    if PARAMS["out_file"]:
        df.to_csv(os.path.join(PARAMS["out_dir"],'benchmark_results.csv'), index=False)

