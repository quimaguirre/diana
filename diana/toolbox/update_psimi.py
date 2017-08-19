import mysql.connector
import re

database = {}
user = {}
password = {}
host = {}
output_file = {}
affinity_file = {}
complementation_file = {}

database['APR_2017'] = 'BIANA_APR_2017'
user['APR_2017'] = 'quim'
password['APR_2017'] = None
host['APR_2017'] = 'localhost'
output_file['APR_2017'] = 'psimi_file_APR_2017.txt'
affinity_file['APR_2017'] = 'affinity_file_APR_2017.txt'
complementation_file['APR_2017'] = 'complementation_file_APR_2017.txt'

database['JAN_2017'] = 'BIANA_JAN_2017'
user['JAN_2017'] = 'quim'
password['JAN_2017'] = None
host['JAN_2017'] = 'localhost'
output_file['JAN_2017'] = 'psimi_file_JAN_2017.txt'
affinity_file['JAN_2017'] = 'affinity_file_JAN_2017.txt'
complementation_file['JAN_2017'] = 'complementation_file_JAN_2017.txt'

database['2016'] = 'BIANA_2016'
user['2016'] = 'quim'
password['2016'] = None
host['2016'] = 'localhost'
output_file['2016'] = 'psimi_file_2016.txt'
affinity_file['2016'] = 'affinity_file_2016.txt'
complementation_file['2016'] = 'complementation_file_2016.txt'

database['2013'] = 'BIANA_MARCH_2013'
user['2013'] = 'biana_user'
password['2013'] = 'biana_password'
host['2013'] = 'sbi.upf.edu'
output_file['2013'] = 'psimi_file_2013.txt'
affinity_file['2013'] = 'affinity_file_2013.txt'
complementation_file['2013'] = 'complementation_file_2013.txt'

key = 'APR_2017'


affinity_dict={
'492':'in vivo',
'493':'in vitro',
'0':'molecular interaction',
'4':'affinity chromatography technology',
'6':'anti bait coimmunoprecipitation',
'7':'anti tag coimmunoprecipitation',
'8':'array technology',
'9':'bacterial display',
'19':'coimmunoprecipitation',
'28':'cosedimentation in solution',
'29':'cosedimentation through density gradient',
'30':'cross-linking',
'34':'display technology',
'47':'far western blotting',
'48':'filamentous phage display',
'49':'filter binding',
'66':'lambda phage display',
'71':'molecular sieving',
'73':'mrna display',
'81':'peptide array',
'84':'phage display',
'89':'protein array',
'92':'protein in situ array',
'95':'proteinchip(r) on a surface-enhanced laser desorption/ionization',
'96':'pull down',
'98':'ribosome display',
'108':'t7 phage display',
'115':'yeast display',
'225':'chromatin immunoprecipitation array',
'400':'affinity technology',
'402':'chromatin immunoprecipitation assay',
'405':'competition binding',
'411':'enzyme linked immunosorbent assay',
'412':'electrophoretic mobility supershift assay',
'413':'electrophoretic mobility shift assay',
'440':'saturation binding',
'657':'systematic evolution of ligands by exponential enrichment',
'676':'tandem affinity purification',
'678':'antibody array',
'695':'sandwich immunoassay',
'729':'luminescence based mammalian interactome mapping',
'813':'proximity enzyme linked immunosorbent assay',
'858':'immunodepleted coimmunoprecipitation',
'892':'solid phase assay',
'899':'p3 filamentous phage display',
'900':'p8 filamentous phage display',
'921':'surface plasmon resonance array',
'946':'ping',
'947':'bead aggregation assay',
'963':'interactome parallel affinity capture',
'1017':'rna immunoprecipitation',
'1028':'modified chromatin immunoprecipitation',
'1029':'proteomics of isolated chromatin segments',
'1031':'protein folding/unfolding',
'1087':'monoclonal antibody blockade'
}

affinity = list(affinity_dict.values())

complementation_dict={
'492':'in vivo',
'493':'in vitro',
'0':'molecular interaction',
'10':	'beta galactosidase complementation',
'11':	'beta lactamase complementation',
'14':	'adenylate cyclase complementation',
'18':	'two hybrid',
'90':	'protein complementation assay',
'97':	'reverse ras recruitment system',
'111':	'dihydrofolate reductase reconstruction',
'112':	'ubiquitin reconstruction',
'228':	'cytoplasmic complementation assay',
'229':	'green fluorescence protein complementation assay',
'230':	'membrane bound complementation assay',
'231':	'mammalian protein protein interaction trap',
'232':	'transcriptional complementation assay',
'369':	'lex-a dimerization assay',
'370':	'tox-r dimerization assay',
'397':	'two hybrid array',
'398':	'two hybrid pooling approach',
'399':	'two hybrid fragment pooling approach',
'432':	'one hybrid',
'437':	'protein tri hybrid',
'438':	'rna tri hybrid',
'588':	'3 hybrid method',
'655':	'lambda repressor two hybrid',
'726':	'reverse two hybrid',
'727':	'lexa b52 complementation',
'728':	'gal4 vp16 complementation',
'809':	'bimolecular fluorescence complementation',
'895':	'protein kinase A complementation',
'916':	'lexa vp16 complementation',
'1037':	'Split renilla luciferase complementation',
}

complementation = list(complementation_dict.values())

cnx = mysql.connector.connect(user=user[key], password=password[key],
                              host=host[key],
                              database=database[key])

cursor = cnx.cursor()

query = ("SELECT k.value, e.value FROM key_attribute_2 AS k JOIN externalEntitypsimi_name AS e ON e.externalEntityID = k.externalEntityID ")

cursor.execute(query)

f = open(output_file[key],"w")
fa = open(affinity_file[key], "w")
fc = open(complementation_file[key], "w")

fa.write('affinity_dict={\n')
fc.write('complementation_dict={\n')

# Pattern of affinity
pa = re.compile('affinity|precipitation')
# Pattern of complementation
pc = re.compile('complementation|hybrid')

for (externalEntityID, value) in cursor:
    f.write("{}\t{}\n".format(externalEntityID, value))
    if value in affinity or pa.search(value) != None:
        fa.write("'{}':\t'{}',\n".format(externalEntityID, value))
    if value in complementation or pc.search(value) != None:
        fc.write("'{}':\t'{}',\n".format(externalEntityID, value))

fa.write('}')
fc.write('}')

f.close()
fa.close()
fc.close()

cursor.close()

cnx.close()