import sys, os, re


class DCDB(object):

    def __init__(self, input_dir):

        self.components_file = os.path.join(input_dir, "COMPONENTS.txt")
        self.dc2components_file = os.path.join(input_dir, "DC_TO_COMPONENTS.txt")
        self.drug_combination_file = os.path.join(input_dir, "DRUG_COMBINATION.txt")
        self.dcc2targets_file = os.path.join(input_dir, "DCC_TO_TARGETS.txt")
        self.targets_file = os.path.join(input_dir, "TARGETS.txt")
        self.dcc2atc_file = os.path.join(input_dir, "DCC_TO_ATC.txt")
        self.atc_file = os.path.join(input_dir, "ATC_CODES.txt")
        self.dcc2outlink_file = os.path.join(input_dir, "DCC_TO_OUTLINK.txt")
        self.outlink_file = os.path.join(input_dir, "DRUG_OUTLINK.txt")
        self.dc2interaction_file = os.path.join(input_dir, "DC_TO_INTERACTION.txt")
        self.interaction_file = os.path.join(input_dir, "DRUG_INTERACTION.txt")
        self.dc2usage_file = os.path.join(input_dir, "DC_TO_DCU.txt")
        self.usage_file = os.path.join(input_dir, "DC_USAGE.txt")

        self.components = set()
        self.component2name = {}
        self.component2target = {}
        self.component2atcid = {}
        self.component2outlink = {}

        self.combinations = set()
        self.combination2component = {}
        self.combination2name = {}
        self.combination2mechanism = {}
        self.combination2usage = {}

        self.targets = set()
        self.target2genename = {}
        self.target2genesymbol = {}
        self.target2geneid = {}
        self.target2uniprotacc = {}
        self.target2taxid = {}
        self.target2hgnc = {}

        self.atcid2atccode = {}

        self.outlink2sourcelink = {}

        self.interactions = set()
        self.interaction2combination = {}
        self.interaction2type = {}
        self.interaction2classification = {}
        self.interaction2component = {}

        self.usages = set()
        self.usage2icd10codes = {}
        self.usage2description = {}
        self.usage2efficacy = {}
        self.usage2effecttype = {}
        self.usage2source = {}

        self.diseases = set()
        self.icd10code2name = {}

        return

    def parse(self):

        ###############################
        #### PARSE COMPONENTS FILE ####
        ###############################

        print("\n.....PARSING DCDB COMPONENTS FILE.....\n")

        with open(self.components_file,'r') as components_file_fd:

            first_line = components_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the component ID and the component name
            # Introduce them in the dictionary component2name: "component_id" => "component_name"
            for line in components_file_fd:
                fields = line.strip().split("\t")
                component = fields[ fields_dict['dcc_id'] ]
                name = fields[ fields_dict['generic_name'] ].lower()
                self.component2name[component] = name
                self.components.add(component)

        #print(self.components)
        #print(self.component2name)


        ####################################
        #### PARSE DC 2 COMPONENTS FILE ####
        ####################################

        print("\n.....PARSING DCDB DC_TO_COMPONENTS FILE.....\n")

        with open(self.dc2components_file,'r') as dc2components_file_fd:

            first_line = dc2components_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the component ID and the combination id
            # Introduce them in the dictionary combination2component: "combination_id" => "[component_id1, ... , component_idn]"
            for line in dc2components_file_fd:
                fields = line.strip().split("\t")
                combination = fields[ fields_dict['dc_id'] ]
                component = fields[ fields_dict['dcc_id'] ]
                self.combination2component.setdefault(combination, set()).add(component)

        #print(self.combination2component)


        ################################
        #### PARSE COMBINATION FILE ####
        ################################

        print("\n.....PARSING DCDB COMBINATION FILE.....\n")

        with open(self.drug_combination_file,'r') as drug_combination_file_fd:

            first_line = drug_combination_file_fd.readline()

            fields_dict = self.obtain_header_fields(first_line)

            for line in drug_combination_file_fd:
                fields = line.strip().split("\t")
                combination = fields[ fields_dict['dc_id'] ]
                comb_names = fields[ fields_dict['brand_name'] ].lower()
                mechanism = fields[ fields_dict['mechanism'] ]
                if mechanism != "null" and mechanism != "":
                    self.combination2mechanism[combination] = mechanism
                if comb_names != "null" and comb_names != "":
                    for comb_name in comb_names.split(';'):
                        self.combination2name.setdefault(combination, set()).add(comb_name)
                self.combinations.add(combination)

        #print(self.combinations)
        #print(self.combination2mechanism)
        #print(self.combination2name)


        ##################################
        #### PARSE DCC 2 TARGETS FILE ####
        ##################################

        print("\n.....PARSING DCDB DCC_TO_TARGETS FILE.....\n")

        with open(self.dcc2targets_file,'r') as dcc2targets_file_fd:

            first_line = dcc2targets_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the component ID and the combination id
            # Introduce them in the dictionary combination2component: "combination_id" => "[component_id1, ... , component_idn]"
            for line in dcc2targets_file_fd:
                fields = line.strip().split("\t")
                component = fields[ fields_dict['dcc_id'] ]
                target = fields[ fields_dict['tar_id'] ]
                self.component2target.setdefault(component, set()).add(target)

        #print(self.component2target)


        ############################
        #### PARSE TARGETS FILE ####
        ############################

        print("\n.....PARSING DCDB TARGETS FILE.....\n")

        with open(self.targets_file,'r') as targets_file_fd:

            first_line = targets_file_fd.readline()

            fields_dict = self.obtain_header_fields(first_line)

            for line in targets_file_fd:
                fields = line.strip().split("\t")
                # There is one target per line and no multiple target attributes
                target = fields[ fields_dict['tar_id'] ]
                genename = fields[ fields_dict['genename'] ]
                gene_symbol = fields[ fields_dict['gene_symbol'] ]
                geneid = fields[ fields_dict['gene_id'] ]
                uniprotacc = fields[ fields_dict['uniprot_accessionnumber'] ]
                taxid = fields[ fields_dict['taxon_number'] ]
                hgnc = fields[ fields_dict['hgnc_id'] ][5:]
                self.targets.add(target)
                if genename != "null" and genename != "":
                    self.target2genename[target] = genename
                if gene_symbol != "null" and gene_symbol != "":
                    self.target2genesymbol[target] = gene_symbol
                if geneid != "null" and geneid != "":
                    self.target2geneid[target] = geneid
                if uniprotacc != "null" and uniprotacc != "":
                    self.target2uniprotacc[target] = uniprotacc
                if taxid != "null" and taxid != "":
                    self.target2taxid[target] = taxid
                if hgnc != "null" and hgnc != "":
                    self.target2hgnc[target] = hgnc

        #print(self.target2genename)
        #print(self.target2genesymbol)
        #print(self.target2geneid)
        #print(self.target2uniprotacc)
        #print(self.target2taxid)
        #print(self.target2hgnc)


        ##############################
        #### PARSE DCC 2 ATC FILE ####
        ##############################

        print("\n.....PARSING DCDB DCC_TO_ATC FILE.....\n")

        with open(self.dcc2atc_file,'r') as dcc2atc_file_fd:

            first_line = dcc2atc_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the component ID and the combination id
            # Introduce them in the dictionary combination2component: "combination_id" => "[component_id1, ... , component_idn]"
            for line in dcc2atc_file_fd:
                fields = line.strip().split("\t")
                component = fields[ fields_dict['dcc_id'] ]
                atcid = fields[ fields_dict['dat_id'] ]
                self.component2atcid.setdefault(component, set()).add(atcid)

        #print(self.component2atcid)


        ##############################
        #### PARSE ATC_CODES FILE ####
        ##############################

        print("\n.....PARSING DCDB ATC CODES FILE.....\n")

        with open(self.atc_file,'r') as atc_file_fd:

            first_line = atc_file_fd.readline()

            fields_dict = self.obtain_header_fields(first_line)

            for line in atc_file_fd:
                fields = line.strip().split("\t")
                atcid = fields[ fields_dict['dat_id'] ]
                atccode = fields[ fields_dict['code'] ]
                self.atcid2atccode[atcid] = atccode

        #print(self.atcid2atccode)


        ##################################
        #### PARSE DCC 2 OUTLINK FILE ####
        ##################################

        print("\n.....PARSING DCDB DCC_TO_OUTLINK FILE.....\n")

        with open(self.dcc2outlink_file,'r') as dcc2outlink_file_fd:

            first_line = dcc2outlink_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the component ID and the combination id
            # Introduce them in the dictionary combination2component: "combination_id" => "[component_id1, ... , component_idn]"
            for line in dcc2outlink_file_fd:
                fields = line.strip().split("\t")
                component = fields[ fields_dict['dcc_id'] ]
                outlink = fields[ fields_dict['dol_id'] ]
                self.component2outlink.setdefault(component, set()).add(outlink)

        #print(self.component2outlink)


        #################################
        #### PARSE DRUG_OUTLINK FILE ####
        #################################

        print("\n.....PARSING DCDB DRUG_OUTLINK FILE.....\n")

        with open(self.outlink_file,'r') as outlink_file_fd:

            first_line = outlink_file_fd.readline()

            fields_dict = self.obtain_header_fields(first_line)

            # Distinct sources --> 'bindingdb', 'pubchem substance', 'pharmgkb', 'wikipedia', 'drugbank', 'kegg compound', 'pubchem compound', 'rxlist', 'chebi', 'pdb', 'kegg drug'

            for line in outlink_file_fd:
                fields = line.strip().split("\t")
                outlink = fields[ fields_dict['dol_id'] ]
                source = fields[ fields_dict['source'] ]
                link = fields[ fields_dict['link'] ]
                self.outlink2sourcelink.setdefault(outlink, {})
                self.outlink2sourcelink[outlink]['source'] = source.lower()
                self.outlink2sourcelink[outlink]['link'] = link

        #print(self.outlink2sourcelink)


        #####################################
        #### PARSE DC 2 INTERACTION FILE ####
        #####################################

        print("\n.....PARSING DCDB DC_TO_INTERACTION FILE.....\n")

        with open(self.dc2interaction_file,'r') as dc2interaction_file_fd:

            first_line = dc2interaction_file_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            # Obtain the interaction ID and the combination id
            # Introduce them in the dictionary combination2interaction: "combination_id" => "[interaction_id1, ... , interaction_idn]"
            for line in dc2interaction_file_fd:
                fields = line.strip().split("\t")
                combination = fields[ fields_dict['dc_id'] ]
                interaction = fields[ fields_dict['di_id'] ]
                self.interaction2combination.setdefault(interaction, [])
                self.interaction2combination[interaction].append(combination)

        #print(self.interaction2combination)


        #####################################
        #### PARSE DRUG INTERACTION FILE ####
        #####################################

        print("\n.....PARSING DCDB DRUG INTERACTION FILE.....\n")

        with open(self.interaction_file,'r') as drug_interaction_file_fd:

            first_line = drug_interaction_file_fd.readline()

            fields_dict = self.obtain_header_fields(first_line)

            for line in drug_interaction_file_fd:
                fields = line.strip().split("\t")
                interaction = fields[ fields_dict['di_id'] ]
                type_int = fields[ fields_dict['type'] ]
                #set(['Pharmacodynical', 'Pharmacokinetical', 'Pharmacodynamical'])
                if type_int == 'Pharmacodynical':
                    type_int = 'Pharmacodynamical'
                classification = fields[ fields_dict['classification'] ]
                #set(['Enhancement of metabolism', 'Same target', 'Inhibition of metabolism', 'Different targets in same biological process', 'anti-acid in blood', 'Different targets of different biological processes', 'Different targets in related biological processes', 'Different targets in different biological processes', 'Inhibtion of metabolism'])
                if classification == 'Inhibtion of metabolism':
                    classification = 'Inhibition of metabolism'
                if classification == 'Different targets of different biological processes':
                    classification = 'Different targets in different biological processes'

                if type_int != "null" and type_int != "":
                    self.interaction2type[interaction] = type_int
                if classification != "null" and classification != "":
                    self.interaction2classification[interaction] = classification
                self.interactions.add(interaction)

                components = fields[ fields_dict['component'] ]
                (component1, component2) = components.split(' - ')

                compid1 = None
                compid2 = None

                component1 = self.check_exception(component1)
                component2 = self.check_exception(component2)

                for component in self.component2name:
                    if self.component2name[component].lower() == component1.lower():
                        compid1 = component
                    elif self.component2name[component].lower() == component2.lower():
                        compid2 = component
                    elif component1.lower() in self.component2name[component].lower():
                        compid1 = component
                    elif component2.lower() in self.component2name[component].lower():
                        compid2 = component

                if compid1 != None and compid2 != None:
                    self.interaction2component.setdefault(interaction, set())
                    self.interaction2component[interaction].add(compid1)
                    self.interaction2component[interaction].add(compid2)
                else:
                    if compid1 == None and compid2 == None:
                        print('The components {} and {} of the drug interaction {} do not have ID\n'.format(component1, component2, interaction))
                    elif compid1 == None and compid2 != None:
                        print('The component {} of the drug interaction {} does not have ID\n'.format(component1, interaction))
                    elif compid2 == None and compid1 != None:
                        print('The component {} of the drug interaction {} does not have ID\n'.format(component2, interaction))
                    sys.exit(10)

        #print(self.interactions)
        #print(self.interaction2component)
        #print(self.interaction2type)
        #print(self.interaction2classification)


        ###############################
        #### PARSE DC 2 USAGE FILE ####
        ###############################

        print("\n.....PARSING DCDB DC_TO_USAGE FILE.....\n")

        with open(self.dc2usage_file,'r') as dc2usage_fd:

            first_line = dc2usage_fd.readline()

            # Obtain a dictionary: "field_name" => "position"
            fields_dict = self.obtain_header_fields(first_line)

            for line in dc2usage_fd:
                fields = line.strip().split("\t")
                combination = fields[ fields_dict['dc_id'] ]
                usage = fields[ fields_dict['dcu_id'] ]
                self.combination2usage.setdefault(combination, set()).add(usage)

        #print(self.combination2usage)


        #############################
        #### PARSE DC_USAGE FILE ####
        #############################

        print("\n.....PARSING DCDB DC_USAGE FILE.....\n")

        sources = set()
        with open(self.usage_file,'r') as usage_fd:

            first_line = usage_fd.readline()

            # DCU_ID	DOSE	ADMINISTRATION	DISEASE	EFFICACY	EFFECT_TYPE	STAGE	SOURCE	ICD10_S	ICD10_L	EFFECTS_T	TOXICITY	OVERALL
            # DCU01786	Didanosine: 250 mg; Lopinavir: 400 mg; Nevirapine: 200 mg; Ritonavi: 100 mr; Zidovudine: 300 mg	null	HIV Infections	Efficacious	Unclear	Phase 2	NCT00109590	B20  Human immunodeficiency virus [HIV] disease	A00-B99  Certain infectious and parasitic diseases; B20-B24  Human immunodeficiency virus (HIV) disease; B20  Human immunodeficiency virus (HIV) disease resulting in infectious and parasitic diseases	At entry, the 169 participants had a median CD4 cell count of 456 cells/mcL and an HIV load of 3.49 log(10) copies/mL. The incidence of mutations in each of the 3 P1032 arms was 0% by sequencing and 1.8%, 7.1%, and 5.3% by OLA in arms A, B, and C, respectively, compared with 13.4% by sequencing and 29.4% by OLA in the comparison group (P < .001 for each study arm vs comparison group).	null	1 month of dual therapy after SD-NVP prevents most NVP resistance to minimal toxicity.
            fields_dict = self.obtain_header_fields(first_line)

            for line in usage_fd:
                fields = line.strip().split("\t")
                usage = fields[ fields_dict['dcu_id'] ]
                # They can be null!!!
                usage_description = fields[ fields_dict['disease'] ].lower() # e.g. ear infections
                efficacy = fields[ fields_dict['efficacy'] ].lower() # set(['non-efficacious', 'efficacious', 'need further study'])
                effect_type = fields[ fields_dict['effect_type'] ].lower() # set(['additive to synergistic', 'unclear', 'additive', 'antagonistic (injection)', 'antagonistic', 'reductive', 'potentiative', 'synergistic'])
                stage = fields[ fields_dict['stage'] ] # e.g. Approved
                source = fields[ fields_dict['source'] ] # e.g. fda orange book, literature curated, clinicaltrials.gov IDs (e.g. NCT00423657, NCT00512707, NCT00825149...)
                icd10_s = fields[ fields_dict['icd10_s'] ] # e.g. I10 Essential (primary) hypertension  # ==> "I10" is the ICD10 ID, "Essential (primary) hypertension", is the IDC10 name
                icd10_l = fields[ fields_dict['icd10_l'] ] # e.g. I00-I99  Diseases of the circulatory system; I10-I15  Hypertensive diseases; I10  Essential (primary) hypertension
                effects_t = fields[ fields_dict['effects_t'] ] # e.g. Mean reductions in SBP and DBP were significantly greater with candesartan/HCTZ 32/25 mg (21/14 mmHg) than with candesartan 32 mg (13/9 mmHg), HCTZ 25 mg (12/8 mmHg) or placebo (4/3 mmHg) [p < 0.001 for all comparisons].
                toxicity = fields[ fields_dict['toxicity'] ] #
                overall = fields[ fields_dict['overall'] ] # e.g. more effective than monotherapy
                #print(usage, disease, efficacy, effect_type)
                #print(stage, source, icd10_s, icd10_l)
                #print(effects_t, toxicity, overall)

                if usage_description == 'null' or usage_description == '' or icd10_s == 'null' or icd10_s == '':
                    continue

                self.usages.add(usage)
                self.usage2description[usage] = usage_description

                icd10_results = icd10_s.lstrip('"').rstrip('"').split('; ')
                for icd10_result in icd10_results:
                    icd10_result = icd10_result.lstrip()
                    if '  ' in icd10_result:
                        icd10_fields = icd10_result.split('  ')
                        icd10_id = icd10_fields[0].upper()
                        icd10_name = icd10_fields[1].lower().rstrip(';') # Remove names with ; in the end
                        self.usage2icd10codes.setdefault(usage, set()).add(icd10_id)
                        self.icd10code2name[icd10_id] = icd10_name
                        self.diseases.add(icd10_id)
                    elif ' ' in icd10_result:
                        icd10_fields = icd10_result.split(' ')
                        icd10_id = icd10_fields[0].upper()
                        icd10_name = ' '.join(icd10_fields[1:]).lower().rstrip(';') # Remove names with ; in the end
                        self.usage2icd10codes.setdefault(usage, set()).add(icd10_id)
                        self.icd10code2name[icd10_id] = icd10_name
                        self.diseases.add(icd10_id)
                    elif '-' in icd10_result:
                        icd10_id = icd10_result.upper()
                        self.usage2icd10codes.setdefault(usage, set()).add(icd10_id)
                        self.diseases.add(icd10_id)
                    else:
                        print('ICD10 in unrecognized format for usage {}: {}'.format(usage, icd10_result))
                        sys.exit(10)

                if efficacy != 'null' and efficacy != '':
                    self.usage2efficacy[usage] = efficacy

                if effect_type != 'null' and effect_type != '':
                    if effect_type == 'antagonistic (injection)': # Edit this type of effect
                        effect_type == 'antagonistic'
                    self.usage2effecttype[usage] = effect_type

                if source != 'null' and source != '':
                    if source.upper().startswith('NCT'):
                        self.usage2source.setdefault(usage, {})
                        self.usage2source[usage]['clinicaltrials.gov'] = source.upper()
                        pass
                    else:
                        self.usage2source[usage] = source.lower()
                        sources.add(source.lower())

        #print(self.usage2icd10codes)

        return


    def check_exception(self, drug_name):
        """ 
        Obtain a dictionary: "field_name" => "position" 
        """
        exception = {
            'Cholecacliferol':'Cholecalciferol',
            'Alemtuzumab':'Campath 1H',
            'Etanercept':'Enbrel',
            'LY294002 ':'LY294002',
            'Progesterone':'Progestrone',
            'Epoetin-alfa':'Epoetin-alpha',
            'Medroxyprogesterone acetate':'Medroxyprogeterone acetate',
            'Risperidone':'Risperidoene',
            'Mecasermin':'Iplex',
            'Insulin lispro':'Humalog',
            'Ec107 fabI':'Ec107fabI',
            '"Fluorouracil': 'Fluorouracil',
            'Idronoxil':'Phenoxodiol'
        }

        if drug_name in exception:
            drug_name = exception[drug_name]
        else:
            if len(drug_name.split(' ')) > 1:
                if drug_name.split(' ')[1] in ('hydrochloride','cypionate','trifenatate','sodium','acetate','kamedoxomil','hydrobromide'):
                    drug_name = drug_name.split(' ')[0]

        return drug_name


    def obtain_header_fields(self, first_line):
        """ 
        Obtain a dictionary: "field_name" => "position" 
        """
        fields_dict = {}

        header_fields = first_line.strip().split("\t")
        for x in range(0, len(header_fields)):
            fields_dict[header_fields[x].lower()] = x

        return fields_dict


if __name__ == "__main__":
    main()