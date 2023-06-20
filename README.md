# CTC-DECONVOLUTION

This tool allows to carry out multiple CTC-deconvolution related functionalities.\
Mainly it allows to make mutation clustering by using the frequencies and the cell number.\
The functionalities are activated via command line and the different parameters are passed in a JSON file.\
The functionalities and its parameters are explained below, including an example JSON file.\

## Installation
The program requires python and multiple libraries which can be installed from the requirements.txt file.\

## Program options:
--j : followed by the json file path with the program configuration (needed for the different options of the program to work correctly)\
(Options described, multiple of these options requires the output produced by a previous option)\
&emsp;--g,--generate : generates clone trees and mutations frequencies with the parameters passed in the JSON file.\
&emsp;--alg : runs the clustering algorithms selected with the parameters passed in the JSON file. The algorithms are indicated with the numbers 1, 2 and 3.\
&emsp;&emsp; 1: It execute the "closest method" clustering algorithm.\
&emsp;&emsp; 2: It execute the "k-means method" clustering algorithm.\
&emsp;&emsp; 3: It execute the "KDE-smoothing method" clustering algorithm.\
&emsp;--vcf : reads the mutation frequencies from a VCF file with the parameters indicated in the parameters in the JSON file.\
&emsp;--build : builds all the statistically possible trees from the frequency list indicated in the parameters in the JSON file.\
&emsp;--build_cluster : builds all the statistically possible trees from the clusterings frequencies from the file indicated in the parameters in the JSON file.\
&emsp;--tree_to_cluster : store in a JSON file the clones generated, with its mutations and frequencies with the --generate option to use as ground truth to calculate the metrics.\
&emsp;--metrics: calculate multiple metrics for benchmarking of the clusterings calculated comparing them with the ground truth generated with the --tree_to_cluster option.\
&emsp;--tree_to_newick : pass the trees given to newick format.\
&emsp;--benchmark : benchmark the different algorithms with the different parameters provided in the JSON file.\
&emsp;--graphs : generate graphs from the results of the benchmark.\

For example, to run the program using a VCF file the input in the command line would be:\
```console
python main.py --j "example.json" --vcf --alg 1 2 3
```

## json configuration file
Here the program functionalities are explained along with its options to configure in the json file
List with the json element necessary for each option of the program:\
Options:\
"set_seed" : boolean, if true it sets the seed for the random functions to the value indicated in the "seed" element.\
"seed" : integer, the seed to set for the random functions.\
"out_file" : boolean, if true it stores the output of the --generate and --alg options.\
"out_dir" : string, the path where the output files will be stored.\


### generates clone trees (--g,--generate) (in a json element called "generate"):
There are 2 different generation algorithms to use:\
"out_file" : the path where the output will be stored\
"fixed" : Boolean, true executes the 2º "fixed" generation function, false executes the 1º generation function.\
1º generation function:\
It generates a tree with using probabilities for the mutations to happen in each clone without more fixed values than the number of nodes or the maximum quantity of mutations.\
&emsp;"nodes" : indicates the number of clones which will be generated (each clone is a node, containing the number of cells, and its mutations)\
&emsp;"mutations" : Indicates the number of mutations which will be generated (some mutations may have a frequency of 0)\
&emsp;"mutation_ratio" : Indicates the probability of the mutations to happen in a genetic clone\
&emsp;"initial_mutation_ratio" : Indicates how possible is a mutation to happen in the initial root mutated parent clone\
&emsp;"base" : Indicates a base number of cells which will have at least each genetic clone\
&emsp;"min_minimum_mutation_range": The minimum number of mutations which will have the minimum mutation range\
&emsp;"max_minimum_mutation_range": The maximum number of mutations which will have the minimum mutation range\
&emsp;"min_population_multiplier_range": The minimum number of cells which will have the population multiplier mutation range\
&emsp;"max_population_multiplier_range": The maximum number of cells which will have the population multiplier mutation range\
2º (fixed) generation function:\
It generates a tree using a fixed quantity of clones (nodes of the tree), cells and mutations.\
&emsp;"fixed_params" : a dictionary with the parameters to use in the fixed generation function:\
&emsp;&emsp;"clone_number" : indicates the number of clones which will be generated (each clone is a node, containing the number of cells, and its mutations)\
&emsp;&emsp;"cell_number" : Indicates the number of cells which will be generated\
&emsp;&emsp;"mutation_number" : Indicates the number of mutations which will be generated\
&emsp;&emsp;"sequencing_depth" : Indicates the sequencing depth simulated (used to simulate the observational noise)\
&emsp;&emsp;"variance" : Indicates the variance used for generating the observational noise\

### runs the clustering algorithms (--alg)(in a json element called "algorithm"):
"data_file" : the path where the input file is stored\
"out_file" : the name of the output file\
"frequency_list_name" : the name of the frequency list in the input file to use\
"KDE": dictionary with the parameters to use in the KDE-smoothing clustering algorithm:\
&emsp;"bandwidth" : the bandwidth to use in the KDE-smoothing clustering algorithm\
&emsp;"search_hyperparameters" : boolean, if true it searches for the best bandwidth for the KDE-smoothing clustering algorithm\
"chromosone_n_sets" : the number of sets of chromosomes to use in the closest clustering algorithm\


### reads the mutation frequencies from a VCF file (--vcf)(in a json element called "VCF"):
"vcf_file" : the path where the input file is stored\
"apply_functions" : boolean, if true it applies the functions indicated in the "functions" element to the mutation frequencies\
"functions" : dictionary with the functions to apply to the dataframe obtained from the vcf file, the key is the name lambda function and the value is the parameters and the axis of the dataframe to apply the function to, for example:\
&emsp;"lambda n,df: df.head(n[0])": {"n":[20],"axis":0}\
"total" : integer, the total number of cells in the sample\
"out_file" : the name of the output file\


### builds all the statistically possible trees from the frequency list indicated (--build)(in a json element called "build_tree"):
"out_file" : the name of the output file\
"in_file" : the path where the input file is stored\
"frequency_list_name" : the name of the frequency list in the input file to use\


### builds all the statistically possible trees from the clusterings frequencies (--build_cluster)(in a json element called "build_cluster"):
"out_file" : the name of the output file\
"in_file" : the path where the input file is stored\

### store in a JSON file the clones generated, with its mutations and frequencies (--tree_to_cluster)(in a json element called "tree_to_cluster"):\
"out_file" : the name of the output file\
"in_file" : the path where the input file is stored\
"frequency_list_name" : the name of the frequency list in the input file to use\
"remove_zeros_from_clusters_freq" : boolean, if true it removes the mutations with frequency 0 from the frequency list if there are\


### calculate multiple metrics for benchmarking of the clusterings calculated (--metrics)(in a json element called "metrics"):
"out_file" : the name of the output file\
"in_file" : the path where the input file is stored\
"ground_truth" : the path where the ground truth file is stored (generated with --tree_to_cluster)\


### pass the trees given to newick format (--tree_to_newick)(in a json element called "tree_to_newick"):
"out_file" : the name of the output file\
"in_file" : the path where the input file is stored\
"remove_root_0" : boolean, if true it removes the root node if its name is "0"\


### benchmark the different algorithms (--benchmark)(in a json element called "benchmark"):
"set_seed" : boolean, if true it sets the seed indicated in the "seed" element\
"seed" : integer, the seed to use if "set_seed" is true\
"out_figures" : boolean, if true it generates and store the figures of the benchmarking\
"store_simulations_and_clustering" : boolean, if true it stores the simulations and the clusterings files generated\
"out_dir" : the path where the output files will be stored
"update_csv_during_run" : boolean, if true it updates the csv file with the results during the run
"update_each_n_runs" : integer, the number of runs to do before updating the csv file with the results
"generate_list": a list of dictionaries where each dictionary contains 2 elements:\
&emsp;"params" : are the parameters to use in the generation of the trees and are the same as in the --generate option\
&emsp;"repeat" : optional parameter, it is an integer which indicates the number of times to repeat the generation of the trees with the same parameters
"frequency_lists" : a list of names of the frequency lists to use in the benchmarking
"KDE": dictionary with the parameters to use in the KDE-smoothing clustering algorithm (the same as in the --alg option)\
"chromosone_n_sets" : the number of sets of chromosomes to use in the closest clustering algorithm (the same as in the --alg option)


### generate graphs from the results of the benchmark (--graphs)(in a json element called "graphs"):
"in_file" : a list of paths of the input files from which to generate the graphs 
"out_dir" : the path where the output files will be stored
"remove_same_seed" : boolean, if true it removes the results with the same seed

## example JSON file
```json
{
    "generate" : {
        "nodes" : [3,15],
        "mutations" : [150,300],
        "mutation_ratio" : [0.2,0.4],
        "initial_mutation_ratio" : [0.2,0.6],
        "base" : [1,6],
        "min_minimum_mutation_range" : 1,
        "max_minimum_mutation_range" : 3,
        "min_population_multiplier_range" : 1,
        "max_population_multiplier_range" : 1,
        "out_file" : "tree.json",
        "fixed" : true,
        "fixed_params" : {
            "clone_number" : [1,5],
            "cell_number" : [2,20],
            "mutation_number" : 50,
            "sequencing_depth" : 100,
            "variance" : 0.001
        }
    },
    "set_seed" : true,
    "seed" : 1234,
    "out_file" : true,
    "out_dir" : "out",

    "algorithm" : {
        "data_file" : "out\\tree.json",
        "out_file" : "cluster.json",
        "frequency_list_name" : "frequency_list",
        "KDE" : {
            "bandwidth" : 0.02,
            "search_hyperparameters" : false
        },
        "chromosone_n_sets" : 2
    },

    "VCF" : {
        "vcf_file" : "vcf\\example.vcf",
        "apply_functions" : true,
        "functions" : {
            "lambda n,df: df.head(n[0])": {"n":[20],"axis":0}
        },
        "total" : 10,
        "out_file" : "vcf.json"

    },

    "build_tree" : {
        "out_file" : "lineage.json",
        "in_file" : "out\\tree.json",
        "frequency_list_name" : "frequency_list"
    },

    "tree_to_cluster" : {
        "out_file" : "ground_truth.json",
        "in_file" : "out\\tree.json",
        "frequency_list_name" : "frequency_list_noisy_no_zeros",
        "remove_zeros_from_clusters_freq" : true
    },

    "metrics" : {
        "out_file" : "performance_metrics.json",
        "in_file" : "out\\cluster.json",
        "ground_truth" : "out\\ground_truth.json"
    },

    "build_cluster" : {
        "out_file" : "cluster_lineages.json",
        "in_file" : "out\\cluster.json"
    },

    "tree_to_newick" : {
        "out_file" : "tree.nwk",
        "in_file" : "out\\cluster_lineages.json",
        "remove_root_0" : true
    },

    "benchmark" : {
        "set_seed" : true,
        "seed" : 1234,
        "out_figures" : false,
        "store_simulations_and_clustering" : false,
        "out_dir" : "benchmark",
        "update_csv_during_run" : true,
        "update_each_n_runs" : 50,
        "generate_list" :
            [
                {
                    "params" : {
                        "nodes" : [3,15],
                        "mutations" : [150,300],
                        "mutation_ratio" : [0.2,0.4],
                        "initial_mutation_ratio" : [0.2,0.6],
                        "base" : [1,6],
                        "min_minimum_mutation_range" : 1,
                        "max_minimum_mutation_range" : 3,
                        "min_population_multiplier_range" : 1,
                        "max_population_multiplier_range" : 1,
                        "out_file" : "tree.json",
                        "fixed" : true,
                        "fixed_params" : {
                            "clone_number" : [1,5],
                            "cell_number" : [2,20],
                            "mutation_number" : 50,
                            "sequencing_depth" : 100,
                            "variance" : 0.001
                        }
                    },
                    "repeat" : 100
                }
            ],
            "frequency_lists" : ["frequency_list","frequency_list_noisy"], 
        "KDE" : [
            {
                "bandwidth" : 0.02,
                "search_hyperparameters" : false
            }
        ],
        "chromosone_n_sets" : 2
    },
    "graphs" : {
        "in_file" : ["benchmark\\benchmark_results.csv"],
        "out_dir" : "benchmark\\graphs",
        "remove_same_seed" : false
    }
}
```
