from collections import defaultdict
import pandas as pd
from scipy.stats import gaussian_kde
#import matplotlib.pyplot as plt
import numpy as np
from genbank import gen_bank

'''
    We're using the Operon model here built from the e.coli 
    baceterium's integernic distances to predict operon membership
    between genes inside of the agro bacterium genome
'''
# build our model
# the threshld for classifying operon membership is 0.60
threshold = 0.60
# recreate the operon model
pos_ctrl, neg_ctrl = [], []
# load the negative and positive controls first
with open('pos_ctrl.txt') as f:
    f = f.read()
    pos_ctrl = [int(num) for num in f.strip('\n').split('\t')]
with open('neg_ctrl.txt') as f:
    f = f.read()
    neg_ctrl = [int(num) for num in f.strip('\n').split('\t')]
    
# generate the log liklihood out of the positive and negative control
LL_h1 = gaussian_kde(pos_ctrl)
LL_h0 = gaussian_kde(neg_ctrl)
def model(x):
    num = LL_h1(x)*0.60
    den = LL_h0(x)*0.40 + num
    return (num/den)
def opern_pred(path:str):
    # get all the genes from each of the replicons and order them by their left position of their first exons. 
    # Then get the distances between each of the gene.
    predictions = []
    result = gen_bank(path)[1]
    # make the results more nicer to access
    result = sorted(result, key=lambda x: x['left'])
    
    # get the distance between every pair of two adjacent genes
    for i in range(len(result)-1): 
        geneA = result[i]
        geneB = result[i+1]
        
        # make sure they're within the same strand
        if geneA['strand'] != geneB['strand']:
            continue

        dist = geneB['left'] - geneA['right'] + 1
        
        pred = True if model(dist)[0] >= threshold else False
        predictions.append({
                'Gene1': geneA['gene_id'],
                'Gene2': geneB['gene_id'],
                'dist': dist,
                'strand': geneA['strand'],
                'pred': pred
            })

    return predictions
