#!/home/fin/anaconda/bin

#####################################################
#       ATB  -  Automated Phylogenetic Binner       #
#####################################################
#                                                   #
#     Copyright (C) - 2014 -  Finlay Maguire        #
#                                                   #
#  Released under an MIT Licence - see LICENSE.txt  #
#                                                   #
#                Use at own risk                    #
#                                                   #
#####################################################

import os                                       #to handle system interactions
import shutil                                   #mainly for rmtree utility
import subprocess                               #to handle calling external programs and multithreading
import ete2                                     #phylogenetic tree parsing
from Bio import Entrez, SeqIO, SearchIO         #to handle sequence data, parse and query taxonomic db
import sklearn.neighbors as skn                 #to categorise sequences using knn classifier
import numpy as np                              #to handle vector operations provided in training dataset
import accessories                              #assistant functions e.g. call_multi, glob etc


def initialise_classification_(user_supplied_folders, categories):

    os.chdir(user_supplied_folders['starting_folder'])

    required_folders = ['sequences',
            'alignment',
            'mask',
            'trees']

    for folder in required_folders:
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)


    classifications_folder = 'classification_categories'

    os.mkdir(classifications_folder)

    os.chdir(classifications_folder)

    for category in categories:
        os.mkdir(category)


class TreeParser():
    """
    Take a tree file and return
    """

    def __init__(self, phylogeny_file, categories, taxdb, logger, n=10):

        self.file_name = phylogeny_file
        self.categories = categories
        self.taxdb = taxdb
        self.logger = logger
        self.n = n

    def parse_phylogeny(self, input_tree):
        '''Function to parse phylogeny using ete2'''

        #ensure tree file is real
        try:
            open(input_tree)
        except IOError:
            print('Tree File Does Not Exist: '+input_tree)

        #read tree file in as ete tree object
        tree = ete2.Tree(input_tree)

        #find node of seed sequence
        seed_node = tree.search_nodes(name="SEED SEQUENCE ")

        #ensure correct node is found and exists
        if len(seed_node) != 0:
            self.logger.error('Seed Node Not Found for Phlyogeny: {}'.format(self.file_name))
            raise ValueError('Seed Node Not Found for Phlyogeny: {}'.format(self.file_name))

        #ensure only one node has been found
        if len(seed_node) > 1:
            self.logger.error('Too Many Seed Nodes Found for Phlyogeny: {}'.format(self.file_name))
            raise ValueError('Too Many Seed Nodes Found for Phlyogeny: {}'.format(self.file_name))

        # return the nearest n relatives (returning a list of node,distance pairs)
        # n being at least number specified or greater or max possible with tree size
        # i.e. if a tree has a 5 leaves, and we specify n of 10, we will return a list of 4 (-1 for seed itself)
        self.seed_node = seed_node
        self.phylogeny = input_tree


    def get_nearest_n_nodes(self):
        '''traverse the tree to get the nearest n nodes'''

        current_node = self.seed_node

        # increase relatives to return by 1 as this will be inclusive of seed
        #in case a tree has fewer leaves than specified number to return
        relatives_to_return = min(len(tree), self.n + 1)

        #a slightly non-pythonic do while
        while True:
            parent_node = current_node.up
            relatives = parent_node.get_leaves()
            current_node = parent_node

            # to be modified as chosen
            # can be more than number of relatives to return
            if len(relatives) >= relatives_to_return:
                relatives_and_seed = relatives

        #get node, distances from seed as k,v pairs and remove seed
        relatives_and_distances = [(node, seed_node.get_distance(node)) \
                for node in relatives_and_seed \
                if seed_node.get_distance(node)!= 0.0]

        #in order to return only the nearest n relatives - sort list of tuples by distance
        #from seed (second element)
        sorted_relative_and_distances = sorted(relatives_and_distances, key = lambda x: x[1])

        #and trim list from closest to n-closest by distance
        nearest_relatives_and_distances = sorted_relative_and_distances[:minimum_relatives_to_return]

        #make sure we haven't fucked up and the list actually contains elements
        if len(nearest_relatives_and_distances) == 0,
            self.logger.error('Error Retrieving Neighbour Nodes: {}'.format(self.file_name))
            raise ValueError('Error Retrieving Neighbour Nodes: {}'.format(self.file_name))

        if len(nearest_relatives_and_distances) > relatives_to_return:
            self.logger.error('Too Many Neighbour Nodes: {}'.format(self.file_name))
            raise ValueError('Too Many Neighbour Nodes: {}'.format(self.file_name))

        return nearest_relatives_and_distances

    def taxonomy_look_up_of_relatives(self, taxdb):
        '''use ete2 NCBITaxa to look up closest relatives to seed'''

        for relative_distance_tuple in relatives:
            nearest_name = taxdb.get_fuzzy_name_translation(relatives_and_distance_tuple[0])

            taxid = taxdb.get_name_translator(nearest_name)

            id_lineage = taxdb.get_lineage(taxid)

            lineage_names = taxdb.get_taxid_translator(id_lineage)

            # get_taxid translator returns a dict so isn't ordered
            # therefore need to iterate through
            species_lineage = []
            for lineage_taxid in id_lineage:
                species_lineage.append(lineage_names[lineage_taxid])

            # in case we have nested id in our look up table want to check taxonomic levels from most specific to least
            species_lineage.reverse()

            # if the taxon level in the lineage is in the look up replace it (in a new list to be safe) with the value in the look up
            for taxonomic_level in species_lineage:
                if taxonomic_level in taxa_lookup.keys():
                    taxonomic_distance_list.append((taxa_lookup[taxonomic_level], relative_distance_tuple[1]))
                    break
                else:
                    continue

        #ensure we haven't lost any elements in taxa lookup
        assert len(taxonomic_distance_list)==len(relatives), 'Relatives Have Been Lost in Taxonomic Lookup: '+relatives

        return taxonomic_distance_list

def create_classifier_vector(nearest_neighbors, categories):
    '''convert the list of tuples of nearest neighbors and distances to n-vectors depending on category number'''

    #init vector of correct length
    state_vec = np.zeros((1, len(categories)))

    for category in categories:
        for neighbor in nearest_neighbors:
            if neighbor[0] == category:
                state_vec[0][index] = state_vec[0][index]+1.0/neighbor[1]
            else:
                assert False, 'Classification not belonging to categories'
        index+=1

    return state_vec






























def parse_training_data(training_data_folder, classification_categories, k, taxa_lookup, user_email):
    '''Read the phylogenies and their classifications into a list of tuples of (state_vec, classification)'''

    training_data = []

    for category in classification_categories:
        folder = training_data_folder + "/" + category
        for tree in os.listdir(folder):
            nearest_relatives = parse_phylogeny(tree, k)
            taxonomic_neighbors = taxonomy_lookup_of_relatives(nearest_relatives, taxa_lookup, user_email)
            tree_state_vector = create_classifier_vector(taxonomic_neighbors, classification_categories)

            training_data.append((tree_state_vector, category))


    return training_data

def train_classifier(training_data_folder, classification_categories):
    '''Function to train classifier using manually checked trees'''
#    skn.KNeighborsClassifier





'''
weights= 'distance'
'''

def classify(classifier_instance, nearest_neighbors, categories, seed_node):
    '''Use classifier on a tree (relatives), classify and move'''

    prediction = classifier_instance.predict([nearest_neighbors])

    #assert prediction is in categories.values(), 'KNN prediction of category failed: ' + seed_node




if __name__=='__main__':
       categories = ['endosymbiont', 'host', 'food']

        #lookup table for taxa querying to label nearest neighbors for classification
        taxa_lookup = {'Archaeplastida':'endosymbiont',
                       'Ciliate':'host',
                       'Bacteria':'food'}


        intialise_analysis()

        train_classifier(training_data_folder, categories)

        for tree in tree_folder:

            #parse_phylogeny()
            #et_nearest_relatives(seed_node, tree, n_neighbors_to_recover)


            # read the tree using ete2 and return a list of the n_nearest_relatives and distances from seed
            nearest_seed_relatives = parse_phylogeny(tree, 5)

            # perform a taxonomic lookup of nearest nodes
            taxonomy_of_seed_relatives = taxonomy_look_up_of_relatives(nearest_seed_relatives,
                    taxa_lookup,
                    user_email)


            # use k-nearest neighbors classifier to
            classification = knn_classification(taxonomy_of_seed_relatives)

            # move phylogeny to appropriate bin based on classification
            binning(tree, classification)


