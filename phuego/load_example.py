# -*- coding: utf-8 -*-

import os
import pandas as pd

def load_test_example():
    test_path = os.path.join(os.path.dirname(__file__), 'data/EGF_vs_Untreated_@min_15_63_240.txt')
    df = pd.read_csv(test_path, sep='\t', header=None)
    return test_path, df