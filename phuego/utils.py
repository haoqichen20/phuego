#!/usr/bin/env python

import math
import numpy as np
from sklearn.neighbors import KernelDensity

def load_gene_names(uniprot_to_gene_path):
	uniprot_to_gene={}
	gene_to_uniprot={}
	localization={}
	f1=open(uniprot_to_gene_path,"r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		uniprot_to_gene[seq[0]]=seq[1].split(";")[0].strip()
		seq=f1.readline()
	return uniprot_to_gene

def denoise_square(G):
	weight=G.strength(G.vs, mode='all', loops=False, weights='weight')
	for i in G.es():
		node_A=i.tuple[0]
		node_B=i.tuple[1]
		den=math.sqrt(weight[node_A]*weight[node_B])
		num=i["weight"]
		G.es[i.index]["weight"]=num/den
	return (G)

def calc_kde(vector):
	obs = len(vector)
	sigma = np.std(vector, ddof=1)
	IQR = (np.percentile(vector, q=75) - np.percentile(vector, q=25)) / 1.3489795003921634
	sigma = min(sigma, IQR)
	if sigma > 0:
		bw= sigma * (obs * 3 / 4.0) ** (-1 / 5)
	else:
		IQR = (np.percentile(vector, q=99) - np.percentile(vector, q=1)) / 4.6526957480816815
		if IQR > 0:
			bw = IQR * (obs * 3 / 4.0) ** (-1 / 5)
	if bw:
		kde = KernelDensity(kernel = 'gaussian', bandwidth=bw).fit(vector.reshape(-1,1))
		grid = np.linspace(min(vector)-1,max(vector)+1,len(vector)*100).reshape(-1,1)
		log_dens = kde.score_samples(grid)
		pdf=np.exp(log_dens)
		grid=grid.ravel()
		normalization=sum(pdf.ravel())
		cdf=np.cumsum(pdf)/normalization
	else:
		onevec=np.ones(len(vector))
		return onevec,onevec

	return cdf,grid