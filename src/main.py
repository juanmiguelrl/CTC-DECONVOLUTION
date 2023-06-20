import argparse
import tools,generator,algorithms,vcf,treebuild,metrics, benchmark, graphs
import os



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--generate", help="Generates a cell lineage with the parameters from the json file", required=False, default=False, action="store_true")
    parser.add_argument("--j", help="json file with the variables", required=True)
    parser.add_argument("--alg", help="cluster algorithm to run", required=False,nargs='+', type=int, default=[])
    #parser.add_argument("--all", help="execute all cluster algorithms", required=False, default=False, action="store_true")
    parser.add_argument("--vcf", help="reads the mutation frequencies from a VCF file", required=False, default=False, action="store_true")
    parser.add_argument("--build", help="bulds a genetic tree from the frequencies given", required=False, default=False, action="store_true")
    parser.add_argument("--tree_to_cluster", help="store the tree nodes into clusters in a json file", required=False, default=False, action="store_true")
    parser.add_argument("--metrics", help="calculate multiple metrics for benchmarking of the clusters calculated comparing them with a ground truth", required=False, default=False, action="store_true")
    parser.add_argument("--build_cluster", help="builds a genetic tree from the frequencies and mutations given from a dictionary of different clusterizations", required=False, default=False, action="store_true")
    parser.add_argument("--tree_to_newick", help="pass the trees given to newick format", required=False, default=False, action="store_true")
    #parser.add_argument("--n","--nni", help="to use nni or not during the training", required=False, default=False, action="store_true")
    parser.add_argument("--benchmark", help="benchmark the different algorithms with the different parameters provided", required=False, default=False, action="store_true")
    parser.add_argument("--graphs", help="generate graphs from the results of the benchmark", required=False, default=False, action="store_true")

    args = parser.parse_args()
    print(args)

    PARAMS = tools.read_json(args.j)

    #if args.n:
    #    PARAMS = tools.nni_params(PARAMS)

    tools.set_seed(PARAMS)

    if args.generate:
        gen_params = PARAMS["generate"]
        if gen_params["fixed"]:  
            #if the parameters are lists, it transforms them into tuples  
            tools.to_tuple(gen_params["fixed_params"])
        else:
            #if the parameters are lists, it transforms them into tuples
            tools.to_tuple(gen_params["random_params"])
            #if the parameters are tuples, it generates a random number between the two values
            #else it keeps the value
            gen_params = tools.check_tuple(gen_params["random_params"])
        resultado = generator.run_generator(gen_params)
        print("The cell lineage has been generated")
        print("The tree is:")
        print(resultado["main_tree"])
        print("The frequency list is:")
        print(resultado["frequency_list"])
        print("The noisy frequency list is:")
        print(resultado["frequency_list_noisy"])
        print("The cluster list is:")
        print(resultado["cluster_list"])
        print("The total number of cells is:")
        print(resultado["total"])
        print("Storing the tree...")
        #t2 = resultado["main_tree"]
        resultado = tools.adapt_output(resultado)
        if PARAMS["out_file"]:
            tools.store_dict(resultado, os.path.join(PARAMS["out_dir"],gen_params["out_file"]))
        else:
            print(resultado)
        #t = tools.create_tree(resultado["main_tree"])
        #t.show()
        #print(tools.compare_trees(t,t2))

    if args.vcf:
        data = vcf.read_vcf(PARAMS["VCF"]["vcf_file"])
        print(data.head())
        if not PARAMS["VCF"].get("frequency_column"):
            data = vcf.info_to_columns(data)
        if PARAMS["VCF"].get("separate_column"):
            for i in range (len(PARAMS["VCF"]["separate_column"])):
                data = vcf.separate_column(data,PARAMS["VCF"]["separate_column"][i],PARAMS["VCF"]["separator"][i])
        #apply functions to the dataframe
        if PARAMS["VCF"]["apply_functions"]:
            data = vcf.apply_functions(data,PARAMS["VCF"]["functions"])
        print(data.head())
        store = vcf.vcf_frequencies_to_dict(data,PARAMS["VCF"]["total"],PARAMS["VCF"].get("frequency_column"))
        tools.store_dict(store, os.path.join(PARAMS["out_dir"],PARAMS["VCF"]["out_file"]))


    if args.alg:
        print("Running the cluster algorithm(s)...")
        cluster = algorithms.run_algorithm(PARAMS,args.alg)
        cluster = tools.adapt_output(cluster)
        print("The cluster algorithm(s) has been run")
        print("The cluster list is:")
        print(cluster)
        if PARAMS["out_file"]:
            print("Storing the cluster list...")
            tools.store_dict(cluster, os.path.join(PARAMS["out_dir"],PARAMS["algorithm"]["out_file"]))

    if args.tree_to_cluster:
        dict_tree = tools.read_json(PARAMS["tree_to_cluster"]["in_file"])["main_tree"]
        frequency_list = tools.read_json(PARAMS["tree_to_cluster"]["in_file"])["frequency_list"]
        clusters_freq = tools.cluster_dict(frequency_list)
        total = 0
        for list in clusters_freq.values():
            total += len(list)
        if PARAMS["tree_to_cluster"]["remove_zeros_from_clusters_freq"]:
            #removes "0.0" from the dictionary
            if 0.0 in clusters_freq:
                del clusters_freq[0.0]
        tree = tools.create_tree(dict_tree)
        cluster = tools.cluster_treelist(tree)
        dict_to_store = {"cluster_list":cluster,"frequency_clusters":clusters_freq,"total":total}
        tools.store_dict(dict_to_store, os.path.join(PARAMS["out_dir"],PARAMS["tree_to_cluster"]["out_file"]))

    if args.build:
        tree = tools.read_json(PARAMS["build_tree"]["in_file"])
        trees = treebuild.build_genetic_trees(tree[PARAMS["build_tree"]["frequency_list_name"]])
        trees = tools.adapt_output(trees)
        print(type(trees))
        tools.store_dict(trees, os.path.join(PARAMS["out_dir"],PARAMS["build_tree"]["out_file"]))

    if args.metrics:
        input = tools.read_json(PARAMS["metrics"]["in_file"])
        ground_truth = tools.read_json(PARAMS["metrics"]["ground_truth"])

        dict_tree = tools.read_json(PARAMS["tree_to_cluster"]["in_file"])["main_tree"]
        #ground_truth_tree = tools.create_tree(dict_tree)

        metrics = metrics.calculate_metrics(input,ground_truth)
        #print(metrics)
        metrics = tools.adapt_output(metrics)
        tools.store_dict(metrics, os.path.join(PARAMS["out_dir"],PARAMS["metrics"]["out_file"]))

    if args.build_cluster:
        input = tools.read_json(PARAMS["build_cluster"]["in_file"])
        dict_trees = treebuild.build_cluster_trees_dict(input)
        dict_trees = tools.adapt_output(dict_trees)
        tools.store_dict(dict_trees, os.path.join(PARAMS["out_dir"],PARAMS["build_cluster"]["out_file"]))

    if args.tree_to_newick:
        input = tools.read_json(PARAMS["tree_to_newick"]["in_file"])
        newick = tools.clusters_to_newick(input,PARAMS["tree_to_newick"]["remove_root_0"])
        print(newick)
        tools.store_string(newick, os.path.join(PARAMS["out_dir"],PARAMS["tree_to_newick"]["out_file"]))

    if args.benchmark:
        print("Running the benchmark...")
        benchmark.run_benchmark(PARAMS["benchmark"])
        print("The benchmark has been run")
        
    if args.graphs:
        print("Generating graphs...")
        csv = PARAMS["graphs"]["in_file"]
        #if csv is a list it will join all the csvs in the list
        if isinstance(csv,list):
            df = tools.join_csvs(csv)
            if PARAMS["graphs"]["remove_same_seed"]:
                df = tools.remove_same_seed(df)
        else:
            df = tools.csv_to_pandas(csv)
        graphs.graphs_from_Dataframe(PARAMS["graphs"]["out_dir"],df)
        print("The graphs have been generated")