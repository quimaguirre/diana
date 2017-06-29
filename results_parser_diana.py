import diana_vbiana as diana
import matplotlib.pyplot as plt
import numpy as np
import re
import time
from os import listdir
from os.path import isfile, isdir, join
from itertools import cycle

from sklearn import svm
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from scipy import interp

analysis_file = "100_analyses/diana_100_analyses_27_10_16.sh"
main_directory = "100_analyses/results"

# Parse the analysis file in order to obtain the names of the DDI pairs and non-DDI pairs
fa = open(analysis_file,"r")
ddi = False
non_ddi = False
ddi_pairs = {}
non_ddi_pairs = {}
pair_regex = re.compile('# "*(.+?)"* vs\. "*(.+?)"*\n')

for line in fa:
    if line == '# DDI EXPERIMENTS\n':
        ddi = True
        non_ddi = False
    elif line == '# NO DDI EXPERIMENTS\n':
        ddi = False
        non_ddi = True

    if ddi == True:
        m = pair_regex.search(line)
        if m != None:
            drug1 = diana.check_special_drug_name(m.group(1))
            drug2 = diana.check_special_drug_name(m.group(2))
            drug_pair = main_directory+'/results_'+drug1+'_'+drug2
            ddi_pairs[drug_pair] = {}
    elif non_ddi == True:
        m = pair_regex.search(line)
        if m != None:
            drug1 = diana.check_special_drug_name(m.group(1))
            drug2 = diana.check_special_drug_name(m.group(2))
            drug_pair = main_directory+'/results_'+drug1+'_'+drug2
            non_ddi_pairs[drug_pair] = {}

fa.close()



# Obtain all the results subfolders of the results main folder
results_dir_list = [join(main_directory, f) for f in listdir(main_directory) if isdir(join(main_directory, f))]

for directory in results_dir_list:

    # Obtain all the result files inside each result subfolder
    results_file_list = [join(directory, f) for f in listdir(directory) if isfile(join(directory, f))]

    results_regex = re.compile('results_.+/results_.+.txt$')
    extended_regex = re.compile('/general_extended_analysis')
    linker_regex = re.compile('/linker_extended_analysis')

    general = {}
    extended = {}
    extended["node"] = []
    extended["edge"] = []
    extended["func"] = []
    linker = {}
    linker["node"] = []
    linker["edge"] = []
    linker["func"] = []
    node = False
    edge = False
    func = False
    seednode = False
    seedfunc = False

    for results_file in results_file_list:

        # The file matches with the GENERAL RESULTS pattern
        if results_regex.search(results_file) != None:
            #print("GENERAL")
            fg = open(results_file,"r")
            for line in fg:
                nline = '## Node profile comparison ##\n'
                eline = '## Edge profile comparison ##\n'
                fline = '## Functional profile comparison ##\n'
                snline = '## Seeds node profile comparison ##\n'
                sfline = '## Seeds functional profile comparison ##\n'
                if line == nline:
                    node = True
                elif line == eline:
                    edge = True
                elif line == fline:
                    func = True
                elif line == snline:
                    seednode = True
                elif line == sfline:
                    seedfunc = True
                if line[0] != "#":
                    if node == True:
                        words = line.split()
                        general["node"] = words # introducing [Spearman, P-value, Dot Product]
                        node = False
                    elif edge == True:
                        words = line.split()
                        general["edge"] = words
                        edge = False
                    elif func == True:
                        words = line.split()
                        general["func"] = words 
                        func = False
                    elif seednode == True:
                        words = line.split()
                        general["seednode"] = words
                        seednode = False
                    elif seedfunc == True:
                        words = line.split()
                        general["seedfunc"] = words
                        seedfunc = False

            # Introduce the general results inside the main dictionary
            if directory in ddi_pairs:
                ddi_pairs[directory]['general'] = general
            elif directory in non_ddi_pairs:
                non_ddi_pairs[directory]['general'] = general

        # The file matches with the GENERAL EXTENDED RESULTS pattern or the LINKER EXTENDED RESULTS pattern
        elif extended_regex.search(results_file) != None or linker_regex.search(results_file) != None:
            if extended_regex.search(results_file) != None:
                ext = True
                link = False
                #print("EXTENDED")
            elif linker_regex.search(results_file) != None:
                link = True
                ext = False
                #print("LINKER")
            fe = open(results_file,"r")
            for line in fe:
                nline = '# NODE PROFILE\n'
                eline = '# EDGE PROFILE\n'
                fline = '# FUNCTIONAL PROFILE\n'
                end = '# END DATAFRAME\n'
                if line == nline:
                    node = True
                elif line == eline:
                    edge = True
                elif line == fline:
                    func = True
                elif line == end:
                    node = edge = func = False
                if line[0] != "#":
                    if node == True:
                        words = line.split()
                        if ext == True:
                            extended["node"].append(words)
                        elif link == True:
                            linker["node"].append(words)
                    elif edge == True:
                        words = line.split()
                        if ext == True:
                            extended["edge"].append(words)
                        elif link == True:
                            linker["edge"].append(words)
                    elif func == True:
                        words = line.split()
                        if ext == True:
                            extended["func"].append(words)
                        elif link == True:
                            linker["func"].append(words)

            # Introduce the general results inside the main dictionary
            if directory in ddi_pairs:
                ddi_pairs[directory]['extended'] = extended
                ddi_pairs[directory]['linker'] = linker
            elif directory in non_ddi_pairs:
                non_ddi_pairs[directory]['extended'] = extended
                non_ddi_pairs[directory]['linker'] = linker


#print(ddi_pairs)
#print(non_ddi_pairs)
print("The number of DDI results introduced is %d" %(len(ddi_pairs)))
print("The number of NON-DDI results introduced is %d" %(len(non_ddi_pairs)))


###############
#### PLOTS ####
###############

# Print the desired results in a plot

results_general = {}
results_general['ddi'] = {}
# results_general['ddi']['node'] = {}
# results_general['ddi']['node']['spearman'] = []
# results_general['ddi']['node']['dot_product'] = []
# results_general['ddi']['edge'] = {}
# results_general['ddi']['edge']['spearman'] = []
# results_general['ddi']['edge']['dot_product'] = []
# results_general['ddi']['func'] = {}
# results_general['ddi']['func']['spearman'] = []
# results_general['ddi']['func']['dot_product'] = []
results_general['non_ddi'] = {}
# results_general['non_ddi']['node'] = {}
# results_general['non_ddi']['node']['spearman'] = []
# results_general['non_ddi']['node']['dot_product'] = []
# results_general['non_ddi']['edge'] = {}
# results_general['non_ddi']['edge']['spearman'] = []
# results_general['non_ddi']['edge']['dot_product'] = []
# results_general['non_ddi']['func'] = {}
# results_general['non_ddi']['func']['spearman'] = []
# results_general['non_ddi']['func']['dot_product'] = []
results_general['ddi']['seednode'] = {}
results_general['ddi']['seednode']['spearman'] = []
results_general['ddi']['seednode']['dot_product'] = []
results_general['ddi']['seedfunc'] = {}
results_general['ddi']['seedfunc']['spearman'] = []
results_general['ddi']['seedfunc']['dot_product'] = []
results_general['non_ddi']['seednode'] = {}
results_general['non_ddi']['seednode']['spearman'] = []
results_general['non_ddi']['seednode']['dot_product'] = []
results_general['non_ddi']['seedfunc'] = {}
results_general['non_ddi']['seedfunc']['spearman'] = []
results_general['non_ddi']['seedfunc']['dot_product'] = []


for result in ddi_pairs:
#     results_general['ddi']['node']['spearman'].append(float(ddi_pairs[result]['general']['node'][0]))
#     results_general['ddi']['node']['dot_product'].append(float(ddi_pairs[result]['general']['node'][2]))
#     results_general['ddi']['edge']['spearman'].append(float(ddi_pairs[result]['general']['edge'][0]))
#     results_general['ddi']['edge']['dot_product'].append(float(ddi_pairs[result]['general']['edge'][2]))
#     results_general['ddi']['func']['spearman'].append(float(ddi_pairs[result]['general']['func'][0]))
#     results_general['ddi']['func']['dot_product'].append(float(ddi_pairs[result]['general']['func'][2]))
    results_general['ddi']['seednode']['spearman'].append(float(ddi_pairs[result]['general']['seednode'][0]))
    results_general['ddi']['seedfunc']['spearman'].append(float(ddi_pairs[result]['general']['seedfunc'][0]))
    results_general['ddi']['seednode']['dot_product'].append(float(ddi_pairs[result]['general']['seednode'][2]))
    results_general['ddi']['seedfunc']['dot_product'].append(float(ddi_pairs[result]['general']['seedfunc'][2]))

for result in non_ddi_pairs:
#     results_general['non_ddi']['node']['spearman'].append(float(non_ddi_pairs[result]['general']['node'][0]))
#     results_general['non_ddi']['node']['dot_product'].append(float(non_ddi_pairs[result]['general']['node'][2]))
#     results_general['non_ddi']['edge']['spearman'].append(float(non_ddi_pairs[result]['general']['edge'][0]))
#     results_general['non_ddi']['edge']['dot_product'].append(float(non_ddi_pairs[result]['general']['edge'][2]))
#     results_general['non_ddi']['func']['spearman'].append(float(non_ddi_pairs[result]['general']['func'][0]))
#     results_general['non_ddi']['func']['dot_product'].append(float(non_ddi_pairs[result]['general']['func'][2]))
    results_general['non_ddi']['seednode']['spearman'].append(float(non_ddi_pairs[result]['general']['seednode'][0]))
    results_general['non_ddi']['seedfunc']['spearman'].append(float(non_ddi_pairs[result]['general']['seedfunc'][0]))
    results_general['non_ddi']['seednode']['dot_product'].append(float(non_ddi_pairs[result]['general']['seednode'][2]))
    results_general['non_ddi']['seedfunc']['dot_product'].append(float(non_ddi_pairs[result]['general']['seedfunc'][2]))


# # Node spearman histogram 
# n, bins, patches = plt.hist([ results_general['ddi']['node']['spearman'], results_general['non_ddi']['node']['spearman'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Node\ profile\ comparison\ (Spearman)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Spearman\'s correlation coefficient')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/node_spear_hist.png')
# #plt.show()

# plt.clf()

# # Node dot product histogram
# n, bins, patches = plt.hist([ results_general['ddi']['node']['dot_product'], results_general['non_ddi']['node']['dot_product'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Node\ profile\ comparison\ (dot\ product)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Dot product')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/node_dotpr_hist.png')
# #plt.show()

# plt.clf()

# # Edge spearman histogram 
# n, bins, patches = plt.hist([ results_general['ddi']['edge']['spearman'], results_general['non_ddi']['edge']['spearman'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Edge\ profile\ comparison\ (Spearman)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Spearman\'s correlation coefficient')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/edge_spear_hist.png')
# #plt.show()

# plt.clf()

# # Node dot product histogram
# n, bins, patches = plt.hist([ results_general['ddi']['edge']['dot_product'], results_general['non_ddi']['edge']['dot_product'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Edge\ profile\ comparison\ (dot\ product)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Dot product')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/edge_dotpr_hist.png')
# #plt.show()

# plt.clf()

# # Functional spearman histogram 
# n, bins, patches = plt.hist([ results_general['ddi']['func']['spearman'], results_general['non_ddi']['func']['spearman'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Functional\ profile\ comparison\ (Spearman)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Spearman\'s correlation coefficient')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/func_spear_hist.png')
# #plt.show()

# plt.clf()

# # Functional dot product histogram
# n, bins, patches = plt.hist([ results_general['ddi']['func']['dot_product'], results_general['non_ddi']['func']['dot_product'] ], 
#                             histtype='bar',
#                             label=['DDI', 'non-DDI'],
#                             color=['green','red'])
# plt.title(r'$\mathrm{Functional\ profile\ comparison\ (dot\ product)}$')
# plt.legend(loc='upper left')
# plt.xlabel('Dot product')
# plt.ylabel('Frequency')
# plt.savefig('100_analyses/func_dotpr_hist.png')
# #plt.show()

# plt.clf()

# Seeds node spearman histogram
n, bins, patches = plt.hist([ results_general['ddi']['seednode']['spearman'], results_general['non_ddi']['seednode']['spearman'] ], 
                            histtype='bar',
                            label=['DDI', 'non-DDI'],
                            color=['green','red'])
plt.title(r'$\mathrm{Seeds\ node\ profile\ comparison\ (Spearman)}$')
plt.legend(loc='upper right')
plt.xlabel('Spearman\'s correlation coefficient')
plt.ylabel('Frequency')
plt.savefig('100_analyses/seednode_spear_hist.png')
#plt.show()

plt.clf()

# Seeds functional spearman histogram
n, bins, patches = plt.hist([ results_general['ddi']['seedfunc']['spearman'], results_general['non_ddi']['seedfunc']['spearman'] ], 
                            histtype='bar',
                            label=['DDI', 'non-DDI'],
                            color=['green','red'])
plt.title(r'$\mathrm{Seeds\ functional\ profile\ comparison\ (Spearman)}$')
plt.legend(loc='upper right')
plt.xlabel('Spearman\'s correlation coefficient')
plt.ylabel('Frequency')
plt.savefig('100_analyses/seedfunc_spear_hist.png')
#plt.show()

plt.clf()

# Seeds node dot product histogram
n, bins, patches = plt.hist([ results_general['ddi']['seednode']['dot_product'], results_general['non_ddi']['seednode']['dot_product'] ], 
                            histtype='bar',
                            label=['DDI', 'non-DDI'],
                            color=['green','red'])
plt.title(r'$\mathrm{Seeds\ node\ profile\ comparison\ (dot\ product)}$')
plt.legend(loc='upper right')
plt.xlabel('Dot product')
plt.ylabel('Frequency')
plt.savefig('100_analyses/seednode_dotpr_hist.png')
#plt.show()

plt.clf()

# Seeds functional dot product histogram
n, bins, patches = plt.hist([ results_general['ddi']['seedfunc']['dot_product'], results_general['non_ddi']['seedfunc']['dot_product'] ], 
                            histtype='bar',
                            label=['DDI', 'non-DDI'],
                            color=['green','red'])
plt.title(r'$\mathrm{Seeds\ functional\ profile\ comparison\ (dot\ product)}$')
plt.legend(loc='upper right')
plt.xlabel('Dot product')
plt.ylabel('Frequency')
plt.savefig('100_analyses/seedfunc_dotpr_hist.png')
#plt.show()


##############################
#### Super Vector Machine ####
##############################

##### Preparing the vector with information (X) and the vector with the classification (y) #####

# mStart = time.time()

# Xvect= []
# Yvect = []

# for pair in ddi_pairs:
#     Yvect.append(1)
#     subvector = []
#     for profile_type in ddi_pairs[pair]['general']:
#         if profile_type != 'seednode' and profile_type != 'seedfunc':
#             subvector.append(float(ddi_pairs[pair]['general'][profile_type][0]))
#             subvector.append(float(ddi_pairs[pair]['general'][profile_type][2]))
#     Xvect.append(subvector)

# for pair in non_ddi_pairs:
#     Yvect.append(0)
#     subvector = []
#     for profile_type in non_ddi_pairs[pair]['general']:
#         if profile_type != 'seednode' and profile_type != 'seedfunc':
#             subvector.append(float(non_ddi_pairs[pair]['general'][profile_type][0]))
#             subvector.append(float(non_ddi_pairs[pair]['general'][profile_type][2]))
#     Xvect.append(subvector)

# X = np.array(Xvect, np.float64)
# y = np.array(Yvect, np.float64)

##### Choosing the best parameters for building the model #####

# # Setting of the parameters example: http://scikit-learn.org/stable/auto_examples/model_selection/grid_search_digits.html#sphx-glr-auto-examples-model-selection-grid-search-digits-py

# # Split the dataset in two equal parts
# X_train, X_test, y_train, y_test = train_test_split(
#     X, Yvect, test_size=0.4, random_state=0)

# # Set the parameters by cross-validation
# tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
#                      'C': [1, 10, 100, 1000]},
#                     {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

# scores = ['precision', 'recall']

# for score in scores:
#     print("# Tuning hyper-parameters for %s" % score)
#     print()

#     clf = GridSearchCV(svm.SVC(C=1), tuned_parameters, cv=5,
#                        scoring='%s_macro' % score)
#     clf.fit(X_train, y_train)

#     print("Best parameters set found on development set:")
#     print()
#     print(clf.best_params_)
#     print()
#     print("Grid scores on development set:")
#     print()
#     means = clf.cv_results_['mean_test_score']
#     stds = clf.cv_results_['std_test_score']
#     for mean, std, params in zip(means, stds, clf.cv_results_['params']):
#         print("%0.3f (+/-%0.03f) for %r"
#               % (mean, std * 2, params))
#     print()

#     print("Detailed classification report:")
#     print()
#     print("The model is trained on the full development set.")
#     print("The scores are computed on the full evaluation set.")
#     print()
#     y_true, y_pred = y_test, clf.predict(X_test)
#     print(classification_report(y_true, y_pred))
#     print()


##### 10-fold cross validation and generation of the ROC curve #####

# Example: http://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py

# cv = StratifiedKFold(n_splits=10)
# classifier = svm.SVC(kernel='linear', probability=True,
#                      C=1000)

# mean_tpr = 0.0
# mean_fpr = np.linspace(0, 1, 100)

# colors = cycle(['cyan', 'indigo', 'seagreen', 'yellow', 'blue', 'darkorange', 'brown', 'grey', 'pink', 'red'])
# lw = 2

# i = 0
# for (train, test), color in zip(cv.split(X, y), colors):
#     print(train)
#     print(test)
#     probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
#     # Compute ROC curve and area the curve
#     fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
#     mean_tpr += interp(mean_fpr, fpr, tpr)
#     mean_tpr[0] = 0.0
#     roc_auc = auc(fpr, tpr)
#     plt.plot(fpr, tpr, lw=lw, color=color,
#              label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

#     i += 1
# plt.plot([0, 1], [0, 1], linestyle='--', lw=lw, color='k',
#          label='Luck')

# mean_tpr /= cv.get_n_splits(X, y)
# mean_tpr[-1] = 1.0
# mean_auc = auc(mean_fpr, mean_tpr)
# plt.plot(mean_fpr, mean_tpr, color='g', linestyle='--',
#          label='Mean ROC (area = %0.2f)' % mean_auc, lw=lw)

# plt.xlim([-0.05, 1.05])
# plt.ylim([-0.05, 1.05])
# plt.xlabel('False Positive Rate')
# plt.ylabel('True Positive Rate')
# plt.title('Receiver operating characteristic example')
# plt.legend(loc="lower right")
# plt.savefig('100_analyses/roc_curve.png')
# #plt.show()








# svm_model = svm.SVC(C=1.5, gamma = 0.1, kernel = 'rbf', probability = False, cache_size=2000 , 
#                     verbose=True)

# svm_model.fit(X, Yvect)

# mEnd = time.time()
# print "Model created in " + repr(mEnd - mStart) + " seconds"

# # Now I should use the model to do predictions
# for pair in testing_set:
#             Xvect= []
#             ann= ""
#                 Xt = np.array(Xvect, np.float64)
#                 #Xt_scaled = scaler.transform(Xt)
#                 predictions = svm_model.predict(Xt) # This fuction uses model to generate a prediction, using the testing set data
#                 for j in range(len(Xt)):
#                     ann += predictions[j]
#                 pred.setAnnotation(ann)
#                 resList.append(pred)



