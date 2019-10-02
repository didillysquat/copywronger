import pandas as pd
import re
import sys
from collections import defaultdict
import argparse
import os

class rrnaNorm:
    def __init__(self):
        self.args = self._define_args()
        self.db_file_path = self.args.db_path
        self.input_count_table_path = self.args.input_path
        self.taxa_to_rrna_copy_dict = {}
        self.taxa_reg_exs = self._generate_taxa_reg_exs()
        self._curate_db()
        self.count_df = None
        self.input_taxa_map = {}
        self._curate_input()
        if self.args.output_path:
            self.output_path = self.args.output_path
        else:
            self.output_path = os.path.join(os.path.dirname(self.input_count_table_path), 'counts_out.tsv')

    def _define_args(self):
        parser = argparse.ArgumentParser(
            description='Simple script to normalise 16S count tables according to an external db',
            epilog='For support email: benjamin.hume@kaust.edu.sa')
        parser.add_argument(
            '--input_path',
            help='Full path to the input count table. Should be given in tab delim format with samples '
                 'in columns and OTUs in rows. The last column of the count table will be excluded from '
                 'normalisation. Order of rows will be maintained', required=True)
        # TODO to improve automatically look for the rrna_db.tsv file in the same directory
        parser.add_argument(
            '--db_path',
            help='Full path to the file that is tab delim file mapping taxa annotations to rrna copy numbers.'
                 ' Additional columns will be ignored.', required=True)

        parser.add_argument(
            '--output_path',
            help='Full path for the output table to be written to.', required=False)
        return parser.parse_args()

    def normalise(self):
        """Do the normalisation of the input count table according to the reference db"""

        print('\nNORMALISING')
        # calculate the grand average of otu copy numbers to divide by if an index map is not possible
        grand_total = sum(self.taxa_to_rrna_copy_dict.values())/len(self.taxa_to_rrna_copy_dict)
        # drop the otu column and add it back on later
        last_col_name = list(self.count_df)[-1]
        OTUs = self.count_df[last_col_name]
        self.count_df.drop(columns=last_col_name, inplace=True)

        # correct on a taxa by taxa basis (index by index)
        # take into account the indexes are not unique
        list_of_index_names = self.count_df.index.values.tolist()
        for i in range(len(list_of_index_names)):
            sys.stdout.write(f'\rCorrecting {list_of_index_names[i]} ({OTUs.iat[i]})')
            if list_of_index_names[i] not in self.input_taxa_map:
                self.count_df.iloc[i] /= grand_total
            else:
                self.count_df.iloc[i] /= self.taxa_to_rrna_copy_dict[self.input_taxa_map[list_of_index_names[i]]]

        # now make relative
        print('\nNormalising count data')
        self.count_df = self.count_df.div(self.count_df.sum(axis=0), axis=1)

        #now reattach the last col
        self.count_df[last_col_name] = OTUs

        # now output
        print(f'Writing out to {self.output_path}')
        self.count_df.to_csv(self.output_path, sep='\t', header=True, index=True)

        print('\n\n\nFINISHED')

    def _curate_db(self):
        """ Read in the db file and keep only rows that contain annotation to rrna copy mappings.
        Discard all other rows.
        Curate the text db so that keys are in a standard format.
        This will be each of the taxa level strings separated by a semicolon. If a taxa level does not have
        a string it will not be part of the key"""

        print('STARTING REFERENCE DB CURATION')

        # read in the file
        with open(self.db_file_path, 'r') as f:
            input_file = f.readlines()

        # keep only the lines that are annotations (will start with a k for kingdom)
        keep_input_file = []
        for line in input_file:
            if line[0] == 'k':
                keep_input_file.append(line.rstrip())

        two_d_list = [line.split('\t')[:2] for line in keep_input_file]
        taxa_str_to_correction_factor = {}
        for two_i_list in two_d_list:
            taxa_str_to_correction_factor[two_i_list[0]] = float(two_i_list[1])

        print('Curating the reference database')
        # now for each key clean it up and add it to the lookup dict
        # kingdom_re = re.compile('k__\[?([\w\s]*)\]?;')
        count_added = 0
        already_present_dict = defaultdict(lambda : 1)

        k_list = list(taxa_str_to_correction_factor.keys())
        for i in range(len(k_list)):
            key_items_list = []
            for reg_ex in self.taxa_reg_exs:
                match = reg_ex.search(k_list[i])
                if match:
                    group = match.groups()[0]

                    # then we found a match and we can add this to the key_items_list
                    key_items_list.append(group)
                else:
                    # we didn't find a match and we need to break out of this and add what we have as the key
                    break
            if key_items_list:
                str_key_to_add = ';'.join(key_items_list)
                if str_key_to_add not in self.taxa_to_rrna_copy_dict:
                    sys.stdout.write('\rAdding: ' + str_key_to_add)
                    self.taxa_to_rrna_copy_dict[str_key_to_add] = taxa_str_to_correction_factor[k_list[i]]
                    count_added += 1
                else:
                    already_present_dict[str_key_to_add] += 1

        print(f'\n\n{count_added} keys were successfuly curated and added to the database')
        print(f'{sum(already_present_dict.values())} of the key instance were duplicated in the database.')
        print(f'These instances represented {len(already_present_dict)} uique keys')
        print('For these keys the first associated value was used.')
        print('The non-unique keys and their abundances are given below:')
        for k, v in already_present_dict.items():
            print(f'\t{k}: {v}')

        print('REFERENCE DB CURATION COMPLETE')

    def _generate_taxa_reg_exs(self):
        reg_exs = [re.compile(r'k__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r'p__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r'c__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r'o__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r'f__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r'g__\[?([\w\s]+)(?=\]?;)'),
                   re.compile(r's__\[?([\w\s]+)(?=\]?)')]
        return reg_exs

    def _curate_input(self):
        """Here we will make a mapping dict of the current tax annotations to the annotation
        that should be used in the rrna copy look up value."""

        print('STARTING COUNT TABLE CURATION')
        print('\n\nNow mapping the taxanomic annotations of your input to those of the newly curated reference db')
        # Read in the input as df and do some curation
        self.count_df = pd.read_csv(self.input_count_table_path, sep='\t')

        self.count_df.set_index('Taxonomy', drop=True, inplace=True)

        # go through the indexes (the taxa annotations) and format them in the same way as the dict

        for taxa_a in self.count_df.index:
            # if taxa_a == 'k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Pseudoalteromonadaceae; g__Pseudoalteromonas; s__porphyrae  ':
            #     poo = 'asdf'
            if taxa_a in self.input_taxa_map:
                # we have already made a map for this taxa annotation
                continue
            key_items_list = []
            for reg_ex in self.taxa_reg_exs:
                match = reg_ex.search(taxa_a)
                if match:
                    group = match.groups()[0].rstrip()
                    # then we found a match and we can add this to the key_items_list
                    key_items_list.append(group)
                else:
                    # we didn't find a match and we need to break out of this and add what we have as the key
                    break
            if key_items_list:
                str_key_to_add =';'.join(key_items_list)
                if str_key_to_add not in self.taxa_to_rrna_copy_dict:
                    # then we are not going to be able to perform a look up
                    short_match_found = False
                    for i in range(-1, -1*len(key_items_list), -1):
                        new_str_to_test = ';'.join(key_items_list[:i])
                        if new_str_to_test in self.taxa_to_rrna_copy_dict:
                            print(f'{str_key_to_add} had to be shortened to {new_str_to_test} (shortened by {-1*i} taxanomic levels)')
                            self.input_taxa_map[taxa_a] = new_str_to_test
                            short_match_found = True
                            break
                    if not short_match_found:
                        print(f'Despite shortening the taxa annotation, no match was found for {taxa_a}')
                        print(f'This taxa will be corrected according to the grand average number of rrna copies')

                else:
                    # put the mapping in
                    self.input_taxa_map[taxa_a] = str_key_to_add
            else:
                raise RuntimeError('Key conversion was not possible.')
        print('COUNT TABLE CURATION COMPLETE')

if __name__ == "__main__":
    rrnanorm = rrnaNorm()
    rrnanorm.normalise()