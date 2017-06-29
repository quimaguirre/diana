###################################

README ANALYSIS OF RESULTS OF DIANA

###################################


1. create_table_results_diana.py
--------------------------------

# Exec

$> python create_table_results_diana.py

# Use

Parse the results from diana and obtain a pandas table with all the results.

# Inputs:

- results_dir 			# Results directory
- data_frame_file 		# Output table of results
- analysis_file 		# Crossings file



2. analyze_features_crossval.py
-------------------------------

# Exec

$> python analyze_features_crossval.py

# Use

Perform an analysis of the relevance of the features. 
Method: Train SVC models using cross-validation, obtain the mean of the AUCs distribution.

# Input parameters

- repetitions = 25 			Number of repetitions
- n_fold = 10     			Number of folds
- num_targets = [3] 		Minimum number of targets for drug
- type_analysis = 'all'		Method of analysis. In this case, we use all kinds of features
- classifier = 'SVC'		Classifier

# Output

- '/tables/feature_selection_crossval.tsv'			Ranking of features 
- '/tables/feature_selection_crossval_groups.tsv'	Ranking of features by group (GUILD-node, GUILD-edge, GUILD-func, Linkers, Seeds)

We obtain a ranking of features by group. How do we do the selection?
From each group, I get the highest ranked feature using spearman and the highest using dot product.
In this way, I obtain a varied group of features.
From Linkers and Seeds, I pick up all the features.



3. choose_classifier.py
-----------------------

# Exec

$> python choose_classifier.py

# Use

Choose the better classifier for our data from a list of different classifiers. 
Method: Train models using cross-validation, compare the distribution of AUCs.

# Input parameters

- repetitions = 25 			Number of repetitions
- n_fold = 10     			Number of folds
- num_targets = [3] 		Minimum number of targets for drug
- type_analysis = 'all'		Method of analysis. In this case, we use all kinds of features
- classifiers: 
        'KNeighbors' : KNeighborsClassifier(3),
        'SVC' : SVC(),
        'SVC linear' : SVC(kernel="linear", C=0.025),
        'SVC rbf' : SVC(gamma=2, C=1),
        'DecisionTree' : DecisionTreeClassifier(max_depth=5),
        'RandomForest' : RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLP' : MLPClassifier(alpha=1),
        'AdaBoost' : AdaBoostClassifier(),
        'GaussianNB' : GaussianNB(),
        'QuadraticDiscr.' : QuadraticDiscriminantAnalysis()


# Output

- '/figures/classifiers_selection_before_tuning.eps'	Distribution of AUCs vs. classifier

We end up choosing the classifier which is more stable (smaller boxplot) and with high median.
In this case we choose SVC.



4. tune_algorithm.py
--------------------

# Exec

$> python tune_algorithm.py

# Use

Choose the best hyperparameters of the chosen algorithm for our data.

# Input parameters

- repetition_analysis = 2 	Number of times that we repeat the analysis
- repetitions = 25 			Number of repetitions of the cross-validation
- n_fold = 10     			Number of folds
- num_targets = [3] 		Minimum number of targets for drug
- type_analysis = 'all'		Method of analysis. In this case, we use all kinds of features
- classifier = 'SVC'		Classifier that we will optimize

# Output

- '/tables/tuning_results.tsv'

We end up choosing the combination of parameters that has been reported as the best performing more times.
In our case, we have chosen the 4 best performing and we will compare them again with 'choose_classifier.py':
- SVC best 1 	{'clf__gamma': 1.0, 'clf__C': 1.0, 'clf__kernel': 'rbf'}		14
- SVC best 2	{'clf__gamma': 0.1, 'clf__C': 10.0, 'clf__kernel': 'rbf'}		9
- SVC best 3 	{'clf__gamma': 0.01, 'clf__C': 1000.0, 'clf__kernel': 'rbf'}	7
- SVC best 4 	{'clf__gamma': 0.1, 'clf__C': 1.0, 'clf__kernel': 'rbf'}		6



5. choose_classifier.py
-----------------------

We will run again choose_classifier.py, but now adding the best combination of SVC parameters to compare them.
The output will be:

- '/figures/classifiers_selection_after_tuning.eps'		Distribution of AUCs vs. classifier

I will choose the 'SVC best 1', which is slightly better than the others



6. analyze_by_targets_and_method_reduced.py
-------------------------------------------

# Exec

$> python analyze_by_targets_and_method_reduced.py

# Use

Analyze the performance of the different methods (guild, linkers, seeds, structure, all) by no. of minimum/maximum targets.
Now we will do it for number of targets equal or greater.

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- num_targets = [3,4,5,6,7,8,9] 											Minimum/Maximum number of targets for drug
- greater_or_smaller = 'greater' 											Get the equal or greater number of targets for drug
- type_analysis = ['all', 'guild', 'linkers', 'seeds', 'struc', 'random']	Method of analysis
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/targets_and_method_greater.eps'

We obtain a plot of the performance of the method depending on the minimum number of targets per drug and the method used.



7. analyze_by_targets_and_method_reduced.py
-------------------------------------------

We will run again analyze_by_targets_and_method_reduced.py, but analyzing by number of targets equal or smaller.
The output will be:

- '/figures/targets_and_method_smaller.eps'



8. analyze_by_targets_and_type_reduced.py
-----------------------------------------

# Exec

$> python analyze_by_targets_and_type_reduced.py

# Use

Analyze the performance of the different types of data (nodes, edges, functions, structure) by no. of minimum/maximum targets.
Now we will do it for number of targets equal or greater.

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- num_targets = [3,4,5,6,7,8,9] 											Minimum/Maximum number of targets for drug
- greater_or_smaller = 'greater' 											Get the equal or greater number of targets for drug
- type_analysis = ['node', 'edge', 'function', 'struc', 'random']			Method of analysis
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/targets_and_type_greater.eps'

We obtain a plot of the performance of the method depending on the minimum number of targets per drug and the method used.



9. analyze_by_targets_and_type_reduced.py
-----------------------------------------

We will run again analyze_by_targets_and_type_reduced.py, but analyzing by number of targets equal or smaller.
The output will be:

- '/figures/targets_and_type_smaller.eps'



10. analyze_by_classification_and_method.py
-------------------------------------------

# Exec

$> python analyze_by_classification_and_method.py

# Use

Analyze the performance of the different methods (guild, linkers, seeds, structure, all) by classification of the drug-drug interactions (how the drugs are related by their targets).

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- classifications = ['Different targets in different biological processes', 
                     'Different targets in related biological processes',
                     'Different targets in same biological process',
                     'Same target']
- type_analysis = ['all', 'guild', 'linkers', 'seeds', 'struc', 'random']	Method of analysis
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/classification_and_method.eps'
- '/tables/classification_and_method_mannwhitney.txt'

We obtain a plot of the performance of the method depending on the classification of the drug-drug interactions.

We also obtain a text file containing the comparisons of All/GUILD/Linkers vs. Seeds for different and related pathways.
To do the comparison we use the Mann-Whitney U test.



11. analyze_by_classification_and_type.py
-------------------------------------------

# Exec

$> python analyze_by_classification_and_type.py

# Use

Analyze the performance of the different types of data (nodes, edges, functions, structure) by classification of the drug-drug interactions (how the drugs are related by their targets).

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- classifications = ['Different targets in different biological processes', 
                     'Different targets in related biological processes',
                     'Different targets in same biological process',
                     'Same target']
- type_analysis = ['node', 'edge', 'function', 'struc', 'random']			Method of analysis
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/classification_and_type.eps'

We obtain a plot of the performance of the type of data depending on the classification of the drug-drug interactions.



12. analyze_by_targets_and_method_and_classification.py
-------------------------------------------------------

# Exec

$> python analyze_by_targets_and_method_and_classification.py

# Use

Analyze the performance of the different methods (guild, linkers, seeds, structure, all) by number of targets of the drug-drug interactions of a certain classification.
In this case we will analyze the targets equal or greater, and the classification by different/related pathways.

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- type_analysis = ['all', 'guild', 'linkers', 'seeds', 'struc', 'random']	Method of analysis
- classification1 = 'Different targets in different biological processes'
- classification2 = 'Different targets in related biological processes'
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/targets_and_method_greater_diffrelpathways.eps'

We obtain a plot of the performance of the method depending on the classification of the drug-drug interactions.

We also obtain a text file containing the comparisons of All/GUILD/Linkers vs. Seeds for different and related pathways.
To do the comparison we use the Mann-Whitney U test.



13. analyze_by_targets_and_type_and_classification.py
-------------------------------------------------------

# Exec

$> python analyze_by_targets_and_type_and_classification.py

# Use

Analyze the performance of the different types of data (nodes, edges, functions, structure) by number of targets of the drug-drug interactions of a certain classification.
In this case we will analyze the targets equal or greater, and the classification by different/related pathways.

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- type_analysis = ['node', 'edge', 'function', 'struc', 'random']			Method of analysis
- classification1 = 'Different targets in different biological processes'
- classification2 = 'Different targets in related biological processes'
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/targets_and_type_greater_diffrelpathways.eps'

We obtain a plot of the performance of the type of data depending on the classification of the drug-drug interactions.

We also obtain a text file containing the comparisons of All/GUILD/Linkers vs. Seeds for different and related pathways.
To do the comparison we use the Mann-Whitney U test.



14. analyze_by_two_ATC.py
-------------------------

# Exec

$> python analyze_by_two_ATC.py

# Use

Analyze the performance of the method 'all' in predicting drug combinations included in specific pairs of ATC with respect to the rest of pairs of ATC.

# Input parameters

- repetitions = 25		 													Number of repetitions of the cross-validation
- n_fold = 2     															Number of folds
- type_analysis = 'all'														Method of analysis
- classifier = 'SVC best 1'													Classifier that we will optimize

# Output

- '/figures/atc_pairs.eps'
- '/tables/ATC_pairs.csv'

We obtain a plot of the performance of the method depending on the ATC pair.

We also obtain a text file containing the comparisons of ATC pairs with respect to the rest.



15. analyze_structural_similarity.py
------------------------------------

# Exec

$> python analyze_structural_similarity.py

# Use

Analyze the structural similarity scores of the drug combinations.
Obtain the pairs of 'me-too' drugs.
Obtain the pairs of 'me-too' drug combinations.

# Input parameters

-

# Output

We obtain a plot of the frequency of drug combinations with certain structural similarity scores:

- '/figures/structural_similarity_in_drug_combinations.eps'

We obtain a pickle file containing a set of me-too drugs:

- me_too_drugs.pcl

We obtain a pickle file containing a set of me-too drug combinations:

- me_too_drug_comb.pcl

We obtain a table containing all the me-too drug pairs

- me_too_drug_pairs.txt

We obtain a table containing all the me-too drug combination pairs

- me_too_drug_comb_pairs.txt



16. analyze_example.py
----------------------

# Exec

$> python analyze_example.py

# Use

Analyze a concrete example of drug pair.
Obtain Venn diagrams of the overlap between targets/nodes/edges/functions.

# Input parameters

The names of the DCDB drugs:

    dcdb_drug1 = 'DCC0073'
    dcdb_drug2 = 'DCC0290'

The threshold and type of threshold to analyze:

    threshold = 10
    type_threshold = 'percentage' # percentage or item

The folders of the results:

    current_dir = '/home/quim/project/diana_results'
    results_dir = current_dir + '/3_targets_analysis'

# Output

Venn diagrams of the overlap between targets/nodes/edges/functions:

    targets_venn = results_dir + '/venn_targets_{}.png'.format(comb_name)
    nodes_venn = results_dir + '/venn_nodes_{}.png'.format(comb_name)
    edges_venn = results_dir + '/venn_edges_{}.png'.format(comb_name)
    functions_venn = results_dir + '/venn_functions_{}.png'.format(comb_name)


