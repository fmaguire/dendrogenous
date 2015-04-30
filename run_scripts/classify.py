#!/usr/bin/env python

from sklearn.externals import joblib
import dendrogenous as dg
import dendrogenous.settings
import dendrogenous.utils
import dendrogenous.core

def main(modelpath, unlabelled_tree_dir):

    classifier = joblib.load(modelpath)

    parser = dg.core.TreeParser(unlabelled_tree_dir)

    tree_matrix = parser.build_matrix()

    file_name_index = parse.file_name_index

    predicted_labels = classifier.predict(tree_matrix)

    dg.utils.bin_trees(zip(index, predicted_labels))

if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])




