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


    with open("run_scripts/test.json") as fh:
        settings = json.load(fh)

    label_def = settings['class_defs']
    test_dir = settings['test_dir']

    tree_parser = dg.core.BuildTest(label_def, test_dir)

    X = tree_parser.build_test()

    print(X)

    with open("X_test.pkl", 'wb') as fh:
        pickle.dump(X, fh, protocol=0)


