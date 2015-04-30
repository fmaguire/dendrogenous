#!/usr/bin/env python

from sklearn.externals import joblib
import dendrogenous as dg
import dendrogenous.settings
import dendrogenous.utils
import dendrogenous.core
import dendrogenous.svm

def main(labelled_tree_directory, modelpath):

    parser = dg.core.TreeParser(labelled_tree_directory)

    tree_matrix, labels = parser.build_training()

    classifier = dg.svm.classify(tree_matrix, labels)

    classifier = joblib.dump(classifier, modelpath)

if __name__=='__main__':

    main(sys.arv[1])

