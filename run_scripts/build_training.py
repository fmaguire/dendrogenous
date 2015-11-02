#!/usr/bin/env python
# Create training data

import dendrogenous as dg
import dendrogenous.settings
import dendrogenous.utils
import dendrogenous.core
import json
import pickle
import sys


if __name__=='__main__':


    with open("run_scripts/training.json") as fh:
        settings = json.load(fh)

    label_def = settings['class_defs']
    label_loc = settings['class_locs']

    tree_parser = dg.core.BuildTraining(label_def, label_loc)

    X, y, encoded_labels = tree_parser.build_training()

    print(X.shape)
    print(y.shape)
    print(encoded_labels)

    with open("X_train.pkl", 'wb') as fh:
        pickle.dump(X, fh, protocol=0)

    with open("y_train.pkl", 'wb') as fh:
        pickle.dump(y, fh, protocol=0)

    with open("label_encoding.pkl", 'wb') as fh:
        pickle.dump(encoded_labels, fh, protocol=0)

