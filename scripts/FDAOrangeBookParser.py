import sys, os, re
import csv


class FDAOrangeBook(object):

    def __init__(self, input_path):

        self.input_path = input_path
        self.products_file = os.path.join(self.input_path, 'products.txt')
        self.patent_file = os.path.join(self.input_path, 'patent.txt')
        self.exclusivity_file = os.path.join(self.input_path, 'exclusivity.txt')

        self.appl_numbers = set()
        self.ingredients = set()
        self.appl_no_to_ingredients = {}
        self.appl_no_to_trade_name = {}
        self.appl_no_to_product_no = {}
        self.appl_no_to_product_type = {}
        self.appl_no_to_approval_date = {}

        self.drugname_to_id = {}
        self.drugname_to_id_to_types = {}

        return


    def parse_products_file(self):
        """
        Parsing of the products file of FDA Orange Book
        """

        print("\n.....PARSING FDA ORANGE BOOK PRODUCTS FILE.....\n")

        num_line = 0

        with open(self.products_file, 'r') as products_fd:

            csvreader = csv.reader(products_fd, delimiter='~')

            for fields in csvreader:

                num_line += 1

                # Obtain a dictionary: "field_name" : "position"
                if num_line == 1:
                    # Ingredient~DF;Route~Trade_Name~Applicant~Strength~Appl_Type~Appl_No~Product_No~TE_Code~Approval_Date~RLD~RS~Type~Applicant_Full_Name
                    # BUDESONIDE~AEROSOL, FOAM;RECTAL~UCERIS~SALIX~2MG/ACTUATION~N~205613~001~~Oct 7, 2014~Yes~Yes~RX~SALIX PHARMACEUTICALS INC
                    fields_dict = self.obtain_header_fields( '~'.join(fields), separator='~')
                    continue

                # Get useful fields
                current_ingredients = fields[ fields_dict['Ingredient'] ].lower() # There can be one ingredient (e.g. BUDESONIDE) or multiple (e.g. SACUBITRIL; VALSARTAN)
                trade_name = fields[ fields_dict['Trade_Name'] ].lower() # UCERIS
                appl_type = fields[ fields_dict['Appl_Type'] ].upper() # The type of new drug application approval. New Drug Applications (NDA or innovator) are "N". Abbreviated New Drug Applications (ANDA or generic) are "A".
                appl_no = fields[ fields_dict['Appl_No'] ] # The FDA assigned number to the application. Format is nnnnnn (e.g. 205613)
                product_no = fields[ fields_dict['Product_No'] ] # The FDA assigned number to identify the application products. Each strength is a separate product.  May repeat for multiple part products. If there are multiple dosages, there will be multiple product numbers because each dosage is part of a different product, but they will still have the same application number (e.g. SOFOSBUVIR; VELPATASVIR). Format is nnn (e.g. 001).
                approval_date = fields[ fields_dict['Approval_Date'] ] # The date the product was approved as stated in the FDA approval letter to the applicant.  The format is Mmm dd, yyyy.  Products approved prior to the January 1, 1982 contain the phrase: "Approved prior to Jan 1, 1982". (e.g. Oct 7, 2014). There can be multiple lines with different dates depending on the dosage (e.g. SOFOSBUVIR; VELPATASVIR)
                product_type = fields[ fields_dict['Type'] ].upper() # The group or category of approved drugs. Format is RX, OTC, DISCN.

                #print(current_ingredients, trade_name, appl_no, product_no)

                # Check if application numbers, ingredients and trade names are available
                if appl_no == '':
                    print('Trade name product {} without application number!'.format(trade_name))
                    sys.exit(10)
                if current_ingredients == '' or trade_name == '': 
                    print('Missing trade name / ingredients for application number {}'.format(appl_no))
                    sys.exit(10)

                # Add application number
                self.appl_numbers.add(appl_no)

                # Check if multiple ingredients
                if '; ' in current_ingredients:
                    current_ingredients = current_ingredients.split('; ')
                elif ';' in current_ingredients:
                    # Examples:
                    # triple sulfa (sulfabenzamide;sulfacetamide;sulfathiazole)
                    # trisulfapyrimidines (sulfadiazine;sulfamerazine;sulfamethazine)
                    # liotrix (t4;t3)
                    #print('Multiple ingredients separated differently: {}'.format(current_ingredients))
                    current_ingredients = current_ingredients.split('; ')
                else:
                    current_ingredients = current_ingredients.split('; ')

                # Add ingredients
                for current_ingredient in current_ingredients:
                    self.appl_no_to_ingredients.setdefault(appl_no, set()).add(current_ingredient)
                    self.ingredients.add(current_ingredient)

                # Add trade name
                self.appl_no_to_trade_name[appl_no] = trade_name

        #print(self.appl_numbers)
        #print(self.ingredients)

        return


    def obtain_header_fields(self, first_line, separator='\t'):
        """ 
        Obtain a dictionary: "field_name" => "position" 
        """
        fields_dict = {}

        header_fields = first_line.strip().split(separator)
        for x in range(0, len(header_fields)):
            fields_dict[header_fields[x]] = x

        return fields_dict


if __name__ == "__main__":
    main()

