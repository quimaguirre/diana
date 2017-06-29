import cPickle
import mysql.connector

def main():

    best_predictions()

    return


def best_predictions():

    # Folder containing everything about the analysis
    current_dir = '/home/quim/project/diana_results'
    results_dir = current_dir + '/3_targets_analysis'

    dump_file = current_dir + "/dcdb2targets.pcl"
    dcdb2targets = cPickle.load(open(dump_file))

    greater_or_smaller = 'greater'
    me_too = True
    me_too_str = ''
    if me_too:
        me_too_str = '_metoo'

    input_file = results_dir + '/tables/drug_comb_scores_list_targets_and_method_{}{}.txt'.format(greater_or_smaller, me_too_str)
    #input_file = results_dir + '/tables/drug_comb_scores_list_classification_and_method{}.txt'.format(me_too_str)
    cnx = mysql.connector.connect(user='quim', password='',
                                  host='localhost',
                                  database='BIANA_JAN_2017')
    dcdb2name = obtain_dcdb2name(cnx)
    dcdb2atc = obtain_dcdb2atc(cnx)

    output_file = results_dir + '/tables/drug_comb_scores_list_targets_and_method_{}{}_final.txt'.format(greater_or_smaller, me_too_str)
    #output_file = results_dir + '/tables/drug_comb_scores_list_classification_and_method{}_final.txt'.format(me_too_str)

    input_fd = open(input_file, 'r')
    output_fd = open(output_file, 'w')

    for line in input_fd:
        if line[0] == '#':
            output_fd.write(line)
            continue
        pair, mean = line.strip().split('\t')
        dcdb1, dcdb2 = pair.split('---')
        name1 = dcdb2name[dcdb1]
        try:
            atcs1 = ','.join(list(obtain_set_of_first_letter_ATC(dcdb2atc[dcdb1])))
        except:
            atcs1 = '-'
        name2 = dcdb2name[dcdb2]
        try:
            atcs2 = ','.join(list(obtain_set_of_first_letter_ATC(dcdb2atc[dcdb2])))
        except:
            atcs1 = '-'
        num_tar1 = len(dcdb2targets[dcdb1])
        num_tar2 = len(dcdb2targets[dcdb2])
        output_fd.write('{}\t{}\t{}\t{}\t{}\t{}\t{:.3f}\n'.format(name1, atcs1, num_tar1, name2, atcs2, num_tar2, float(mean)))

    input_fd.close()
    output_fd.close()

    return

def obtain_dcdb2name(cnx):
    """
    Obtain dictionary DCDB_drugID : drug_name
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, N.value FROM externalEntityDCDB_drugID D, externalEntityName N WHERE D.externalEntityID = N.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2name = {}
    for dcdb, name in cursor:
        dcdb2name[dcdb] = name

    cursor.close()

    return dcdb2name

def obtain_dcdb2atc(cnx):
    """
    Obtain dictionary DCDB_drugID : set(ATCs)
    """

    cursor = cnx.cursor()

    query = (''' SELECT D.value, A.value FROM externalEntityDCDB_drugID D, externalEntityATC A WHERE D.externalEntityID = A.externalEntityID
             ''')

    cursor.execute(query)

    dcdb2atc = {}
    for dcdb, atc in cursor:
        dcdb2atc.setdefault(dcdb, set())
        dcdb2atc[dcdb].add(atc)

    cursor.close()

    return dcdb2atc

def obtain_set_of_first_letter_ATC(ATC_set):
    """ set(['C10AB05','J01XE01']) --> set(['C','J']) """
    one_letter = set()
    for atc in ATC_set:
        one_letter.add(atc[0])
    return one_letter


if  __name__ == "__main__":
    main()
