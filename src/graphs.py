import tools
import matplotlib.pyplot as plt
import pandas as pd
import os


def graphs_from_Dataframe(dir,df):

    unique_simulation = df.drop_duplicates(subset='simulation')


    # Define a color palette for the 3 algorithms
    color_palette = ['orange', 'green', 'blue']
    number_color = "red"

    num_figures = 11
    figs = []
    for _ in range(num_figures):
        fig, axs = plt.subplots()
        figs.append(axs)

    #generate a histogram of number of clusters per simulation
    #figs[0].hist(unique_simulation["nodes"],len(unique_simulation["nodes"].unique()),edgecolor='black')
    figs[0].bar(unique_simulation['nodes'].unique(),unique_simulation['nodes'].value_counts(), width=1,edgecolor='black')
    # Set labels and title
    figs[0].set_xlabel('Number of Clones')
    figs[0].set_ylabel('Number of Samples')
    figs[0].set_title('Histogram of Clones')
    figs[0].set_xticks(unique_simulation["nodes"])
    figs[0].set_xticklabels(unique_simulation["nodes"])

    #generate a histogram of number of mutations per simulation
    #unique_simulation["total_number_of_mutations"].hist()
    #figs[1].hist(unique_simulation["total_number_of_mutations"],len(unique_simulation["total_number_of_mutations"].unique())
    #                ,width=0.8)
    figs[1].bar(unique_simulation['total_number_of_mutations'].unique(),unique_simulation['total_number_of_mutations'].value_counts(), width=1,edgecolor='black')
    # Set labels and title
    figs[1].set_xlabel('Number of mutations')
    figs[1].set_ylabel('Number of Samples')
    figs[1].set_title('Histogram of mutations')
    figs[1].set_xticks(unique_simulation["total_number_of_mutations"])
    figs[1].set_xticklabels(unique_simulation["total_number_of_mutations"],rotation=60)

    #generate of the number of cells per simulation
    #unique_simulation["total_number_of_cells"].hist()
    #figs[2].hist(unique_simulation["total_number_of_cells"],len(unique_simulation["total_number_of_cells"].unique())
    #                ,width=0.8)
    #figs[2].bar(unique_simulation["total_number_of_cells"].groupby(unique_simulation["total_number_of_cells"]).count(),unique_simulation['total_number_of_cells'].value_counts().sort_index(), width=2,edgecolor='black')
    figs[2].bar(unique_simulation['total_number_of_cells'].unique(),unique_simulation['total_number_of_cells'].value_counts(), width=1,edgecolor='black')
    # Set labels and title
    figs[2].set_xlabel('Number of cells')
    figs[2].set_ylabel('Number of Samples')
    figs[2].set_title('Histogram of cells')
    figs[2].set_xticks(unique_simulation["total_number_of_cells"])
    figs[2].set_xticklabels(unique_simulation["total_number_of_cells"],rotation=60)

    """ Without noise neither zeros"""
    #without noise and for the different clusters by real frequency (we take into account the clusters by the frequency,
    # as if the frequency is the same for 2 clusters in this case are indistinguishable without more information)
    frequency_list_no_zeros = df[(df['frequency_list'] == 'frequency_list') & (df['clusters_ground_truth'] == 'real_clusters')]

    #add a column with the clustering algorithm
    frequency_list_no_zeros['clustering_algorithm'] = frequency_list_no_zeros['algorithm'].apply(lambda x: int(x.split("_")[1]))
    # Define the mapping dictionary
    mapping = {1: 'Closest', 2: 'k-means', 3: 'KDE-smoothing'}

    # Map the integer values to names
    frequency_list_no_zeros['alg_name'] = frequency_list_no_zeros['clustering_algorithm'].map(mapping)
    #take from k means the best result of each simulation for each metric
    #for that it takes the best value of the rows with the same simulation and same frequency_list
    #taking the minimum for clone number error
    #so it takes the row with the minimum value of clone number error for each simulation and each frequency_list of the ones which start with "Algorithm_2_kmeans" in its algorithm column


    # Assuming you have a DataFrame named 'df'
    df_alg2 = frequency_list_no_zeros[frequency_list_no_zeros['clustering_algorithm'] == 2]  # Select rows where 'algorithm' is 'alg2'
    #create a column with the absolute values of the column "clone" inside each group
    df_alg2['abs_clone_number_error'] = df_alg2.groupby('simulation')['clone_number_error'].transform(lambda x: x.abs())
    #filters the dataframe to get the rows with the minimum absolute value in each group
    filtered_df = df_alg2[df_alg2['v_measure'] == df_alg2.groupby('simulation')['v_measure'].transform('max')]
    #for each group it takes the row with the maximum value of the sum of the columns v_measure and adjusted_rand_score
    filtered_df["sum_v_adjusted_rand_score"] = filtered_df["v_measure"] + filtered_df["adjusted_rand_score"]
    filtered_df = filtered_df[filtered_df['sum_v_adjusted_rand_score'] == filtered_df.groupby('simulation')['sum_v_adjusted_rand_score'].transform('max')]
    filtered_df = filtered_df[filtered_df['abs_clone_number_error'] == filtered_df.groupby('simulation')['abs_clone_number_error'].transform('min')]
    #for each group it takes the row with the minimum value of the column cluster_number
    filtered_df = filtered_df.loc[filtered_df.groupby('simulation')['k_means_number_of_clusters'].idxmin()]
    #remove the column with the absolute values of the column "clone" inside each group and with the sum of v_measure and adjusted_rand_score
    filtered_df = filtered_df.drop('abs_clone_number_error', axis=1)
    filtered_df = filtered_df.drop('sum_v_adjusted_rand_score', axis=1)
    # Sort the DataFrame by the absolute difference in ascending order
    #df_sorted = df_alg2.sort_values('difference')
    # Select the row with the closest value to zero among 'alg2' rows
    #closest_to_zero = df_sorted.iloc[0]
    # Filter the original DataFrame to keep rows with 'alg1' and 'alg3' algorithms
    df_filtered = frequency_list_no_zeros[frequency_list_no_zeros['clustering_algorithm'].isin([1, 3])]
    # Concatenate the rows with 'alg2' closest to zero and the filtered rows with 'alg1' and 'alg3'
    frequency_list_no_zeros_kmeans_best_clone_number_error = pd.concat([filtered_df, df_filtered])


    #frequency_list_no_zeros_kmeans_best_clone_number_error = frequency_list_no_zeros[frequency_list_no_zeros['algorithm'] == 2].groupby(['simulation','frequency_list']).min()['clone_number_error']


    #frequency_list_no_zeros.groupby(['simulation','frequency_list']).min().groupby('clustering_algorithm').mean().plot.bar(ax=figs[3])
    #figs[3].bar(frequency_list_no_zeros.groupby('clustering_algorithm').mean().index,frequency_list_no_zeros.groupby('clustering_algorithm').mean()['clone_number_error'])
    #figs[3].bar(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['clone_number_error'],edgecolor='black', color=color_palette)
    # Crear el boxplot
    labels = frequency_list_no_zeros_kmeans_best_clone_number_error['alg_name'].unique()

    figs[3].violinplot(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['clone_number_error'].apply(list).values.tolist())
    frequency_list_no_zeros_kmeans_best_clone_number_error.boxplot(column="clone_number_error",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[3])
    figs[3].get_figure().suptitle('')

    figs[3].set_title('Clone_number_error of the metrics of each clustering algorithm\n without noise')
    figs[3].set_ylabel('Clone_number_error')
    figs[3].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['clone_number_error'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[3].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[3].containers:
    #    figs[3].bar_label(bars, label_type='center', color=number_color)


    #for v_measure_score
    #figs[4].bar(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['v_measure'],edgecolor='black', color=color_palette)
    figs[4].violinplot(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['v_measure'].apply(list).values.tolist())
    frequency_list_no_zeros_kmeans_best_clone_number_error.boxplot(column="v_measure",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[4])
    figs[4].get_figure().suptitle('')
    figs[4].set_title('V_measure of the metrics of each clustering algorithm\n without noise')
    figs[4].set_ylabel('V_measure')
    figs[4].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['v_measure'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[4].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[4].containers:
    #    figs[4].bar_label(bars, label_type='center', color=number_color)
    #for adjusted_rand_score
    #figs[5].bar(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['adjusted_rand_score'],edgecolor='black', color=color_palette)
    figs[5].violinplot(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['adjusted_rand_score'].apply(list).values.tolist())
    frequency_list_no_zeros_kmeans_best_clone_number_error.boxplot(column="adjusted_rand_score",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[5])
    figs[5].get_figure().suptitle('')
    figs[5].set_title('Adjusted_rand_score of the metrics of each clustering algorithm\n without noise')
    figs[5].set_ylabel('Adjusted_rand_score')
    figs[5].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['adjusted_rand_score'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[5].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[5].containers:
    #    figs[5].bar_label(bars, label_type='center', color=number_color)
    #for time
    #figs[6].bar(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['time'],edgecolor='black', color=color_palette)
    figs[6].violinplot(frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['time'].apply(list).values.tolist())
    frequency_list_no_zeros_kmeans_best_clone_number_error.boxplot(column="time",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[6])
    figs[6].get_figure().suptitle('')
    figs[6].set_title('Time used of each clustering algorithm\n without noise')
    figs[6].set_ylabel('Time used')
    figs[6].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['time'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[6].text(i+1, mean, f'mean={mean:.8f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[6].containers:
    #    figs[6].bar_label(bars, label_type='center', color=number_color)
    #barplot of the mean of the metrics of each clustering algorithm\n
    #frequency_list_no_zeros.groupby('clustering_algorithm').mean().plot.bar(ax=figs[3])




    """ With noise and without zeros"""

    #same as before but with frequency_list_noisy_no_zeros
    #here real_clusters are used, as due to the noise probably 2 frequencies are not gonna repeat
    frequency_list_noisy_no_zeros = df[(df['frequency_list'] == 'frequency_list_noisy') & (df['clusters_ground_truth'] == 'real_clusters')]
    #add a column with the clustering algorithm
    frequency_list_noisy_no_zeros['clustering_algorithm'] = frequency_list_noisy_no_zeros['algorithm'].apply(lambda x: int(x.split("_")[1]))
    # Define the mapping dictionary
    mapping = {1: 'Closest', 2: 'k-means', 3: 'KDE-smoothing'}

    # Map the integer values to names
    frequency_list_noisy_no_zeros['alg_name'] = frequency_list_noisy_no_zeros['clustering_algorithm'].map(mapping)
    # Assuming you have a DataFrame named 'df'
    df_alg2 = frequency_list_noisy_no_zeros[frequency_list_noisy_no_zeros['clustering_algorithm'] == 2]  # Select rows where 'algorithm' is 'alg2'
    #create a column with the absolute values of the column "clone" inside each group
    df_alg2['abs_clone_number_error'] = df_alg2.groupby('simulation')['clone_number_error'].transform(lambda x: x.abs())
    #filters the dataframe to get the rows with the minimum absolute value in each group
    #filtered_df = df_alg2[df_alg2['abs_clone_number_error'] == df_alg2.groupby('simulation')['abs_clone_number_error'].transform('min')]
    filtered_df = df_alg2[df_alg2['v_measure'] == df_alg2.groupby('simulation')['v_measure'].transform('max')]
    #for each group it takes the row with the maximum value of the sum of the columns v_measure and adjusted_rand_score
    #filtered_df["sum_v_adjusted_rand_score"] = filtered_df["v_measure"] + filtered_df["adjusted_rand_score"]
    #filtered_df = filtered_df[filtered_df['sum_v_adjusted_rand_score'] == filtered_df.groupby('simulation')['sum_v_adjusted_rand_score'].transform('max')]
    filtered_df["sum_v_adjusted_rand_score"] = filtered_df["v_measure"] + filtered_df["adjusted_rand_score"]
    filtered_df = filtered_df[filtered_df['sum_v_adjusted_rand_score'] == filtered_df.groupby('simulation')['sum_v_adjusted_rand_score'].transform('max')]
    filtered_df = filtered_df[filtered_df['abs_clone_number_error'] == filtered_df.groupby('simulation')['abs_clone_number_error'].transform('min')]
    #for each group it takes the row with the minimum value of the column cluster_number
    filtered_df = filtered_df.loc[filtered_df.groupby('simulation')['k_means_number_of_clusters'].idxmin()]
    #remove the column with the absolute values of the column "clone" inside each group and with the sum of v_measure and adjusted_rand_score
    filtered_df = filtered_df.drop('abs_clone_number_error', axis=1)
    filtered_df = filtered_df.drop('sum_v_adjusted_rand_score', axis=1)
    # Sort the DataFrame by the absolute difference in ascending order
    #df_sorted = df_alg2.sort_values('difference')
    # Select the row with the closest value to zero among 'alg2' rows
    #closest_to_zero = df_sorted.iloc[0]
    # Filter the original DataFrame to keep rows with 'alg1' and 'alg3' algorithms
    df_filtered = frequency_list_noisy_no_zeros[frequency_list_noisy_no_zeros['clustering_algorithm'].isin([1, 3])]
    # Concatenate the rows with 'alg2' closest to zero and the filtered rows with 'alg1' and 'alg3'
    frequency_list_noisy_no_zeros_kmeans_best_clone_number_error = pd.concat([df_filtered, filtered_df])

    #figs[7].bar(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['clone_number_error'],edgecolor='black', color=color_palette)
    labels = frequency_list_noisy_no_zeros_kmeans_best_clone_number_error['alg_name'].unique()
    
    figs[7].violinplot(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['clone_number_error'].apply(list).values.tolist())
    frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.boxplot(column="clone_number_error",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[7])
    figs[7].get_figure().suptitle('')
    figs[7].set_title('Clone_number_error of the metrics of each clustering algorithm\n with noise')
    figs[7].set_ylabel('Clone_number_error')
    figs[7].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['clone_number_error'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[7].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[7].containers:
    #    figs[7].bar_label(bars, label_type='center', color=number_color)

    #for v_measure_score
    #figs[8].bar(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['v_measure'],edgecolor='black', color=color_palette)
    figs[8].violinplot(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['v_measure'].apply(list).values.tolist())
    frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.boxplot(column="v_measure",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[8])
    figs[8].get_figure().suptitle('')
    figs[8].set_title('V_measure of the metrics of each clustering algorithm\n with noise')
    figs[8].set_ylabel('V_measure')
    figs[8].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['v_measure'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[8].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[8].containers:
    #    figs[8].bar_label(bars, label_type='center', color=number_color)

    #for adjusted_rand_score
    #figs[9].bar(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['adjusted_rand_score'],edgecolor='black', color=color_palette)
    figs[9].violinplot(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['adjusted_rand_score'].apply(list).values.tolist())
    frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.boxplot(column="adjusted_rand_score",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[9])
    figs[9].get_figure().suptitle('')
    figs[9].set_title('Adjusted_rand_score of the metrics of each clustering algorithm\n with noise')
    figs[9].set_ylabel('Adjusted_rand_score')
    figs[9].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['adjusted_rand_score'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[9].text(i+1, mean, f'mean={mean:.3f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    #for bars in figs[9].containers:
    #    figs[9].bar_label(bars, label_type='center', color=number_color)
    #for time
    #figs[10].bar(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean().index,frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name').mean()['time'],edgecolor='black', color=color_palette)
    figs[10].violinplot(frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['time'].apply(list).values.tolist())
    frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.boxplot(column="time",
                     by='alg_name',patch_artist=True,boxprops={'alpha': 0.4},ax=figs[10])
    figs[10].get_figure().suptitle('')
    figs[10].set_title('Time used of each clustering algorithm\n with noise')
    figs[10].set_ylabel('Time used')
    figs[10].set_xlabel('Clustering algorithm')
    # Calculate the mean values for each group
    mean_values = frequency_list_noisy_no_zeros_kmeans_best_clone_number_error.groupby('alg_name')['time'].mean()
    # Add mean value annotations to the plot
    for i, mean in enumerate(mean_values):
        figs[10].text(i+1, mean, f'mean={mean:.8f}', ha='center', va='bottom',color=number_color,fontweight='bold')
    for bars in figs[10].containers:
        figs[10].bar_label(bars, label_type='center', color=number_color)


    #store the figures in the directory dir
    for i, fig in enumerate(figs):
        #the name is the title of the figure
        fig_name = str(i)+"_"+figs[i].get_title() + ".png"
        fig_name = fig_name.replace("\n","")
        fig_path = os.path.join(dir,fig_name)
        fig.figure.tight_layout()
        tools.store_fig_ax(fig,fig_path)

    plt.show()

