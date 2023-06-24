# -*- coding: utf-8 -*-

import numpy as np
import igraph as ig

def rwr_values(network, graph_nodes, ini_pos, ini_neg, seeds, seeds_pos, 
               seeds_neg, network_path, network_random_path, damping, res_folder):
    
	number_of_nodes=network.vcount()
	# rwr_random_values={}
	empirical_rwr=np.zeros((6,number_of_nodes),dtype=float)

	empirical_values={}
	pvalues={}
	for i in network.vs:
		pvalues[i["name"]]=np.zeros(6,dtype=int)
		empirical_values[i["name"]]=np.zeros(6,dtype=float)

	pos=False
	if ini_pos[0]!="False":
		pos=True
		to_delete=[]
		for j in network.vs.select(name_in=ini_pos):
			to_delete.append(j.index)

	neg=False
	if ini_neg[0]!="False":
		neg=True
		to_delete=[]
		for j in network.vs.select(name_in=ini_neg):
			to_delete.append(j.index)
	flag_pos=0
	flag_neg=0
	for i in enumerate(seeds):
		if len(i[1])>0:
			if pos==True and flag_pos==0 and i[0]<=2:
				network.delete_vertices(to_delete)
				flag_pos=1
			if neg==True and flag_neg==0 and i[0]>2:
				network.delete_vertices(to_delete)
				flag_neg=1

			number_of_nodes=network.vcount()
			reset_vertex=np.zeros(number_of_nodes)

			for j in network.vs.select(name_in=i[1]):
				if i[0]<=2:
					reset_vertex[j.index]=seeds_pos[j["name"]]
				else:
					reset_vertex[j.index]=seeds_neg[j["name"]]
			pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight',implementation='prpack'))
			if flag_pos==1 and i[0]<=2:
				for j in to_delete:
					pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

			if flag_neg==1 and i[0]>2:
				for j in to_delete:
					pagerank=np.concatenate((pagerank[:j], [0], pagerank[j:]))

			empirical_rwr[i[0]]=pagerank

		if i[0]==2 and pos==True:
			network = ig.Graph.Read_Ncol(network_path, weights=True, directed=False)
			pos=False
			flag_pos=0
			to_delete=[]


	for i in enumerate(empirical_rwr.T):
		empirical_values[graph_nodes[i[0]]]=i[1]

	for ii in range(1000):
		network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
		nodes=network_random.vs["name"]
		number_of_nodes=network_random.vcount()
		random_rwr=np.zeros((6,number_of_nodes),dtype=float)
		pos=False
		if ini_pos[0]!="False":
			to_delete=[]
			for j in network_random.vs.select(name_in=ini_pos):

				to_delete.append(j.index)
			pos=True
		neg=False
		if ini_neg[0]!="False":
			to_delete=[]
			for j in network_random.vs.select(name_in=ini_neg):
				to_delete.append(j.index)
			neg=True
		flag_pos=0
		flag_neg=0
		#st = time.time()

		for jj in enumerate(seeds):
			if len(jj[1])>0:
				if pos==True and flag_pos==0 and jj[0]<=2:
					network_random.delete_vertices(to_delete)
					flag_pos=1
				if neg==True and flag_neg==0 and jj[0]>2:
					network_random.delete_vertices(to_delete)
					flag_neg=1

				number_of_nodes=network_random.vcount()
				reset_vertex=np.zeros(number_of_nodes)

				for j in network_random.vs.select(name_in=jj[1]):
					if jj[0]<=2:
						reset_vertex[j.index]=seeds_pos[j["name"]]
					else:
						reset_vertex[j.index]=seeds_neg[j["name"]]
				prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight',implementation='prpack'))

				if flag_pos==1 and jj[0]<=2:
					for j in to_delete:
						prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

				if flag_neg==1 and jj[0]>2:
					for j in to_delete:
						prandom=np.concatenate((prandom[:j], [0], prandom[j:]))

				if jj[0]==2 and pos==True:
					network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
					flag_pos=0
					pos=False
					to_delete=[]

				random_rwr[jj[0]]=prandom
			else:
				if jj[0]==2:
					network_random = ig.Graph.Read_Ncol(network_random_path+str(ii)+".txt", weights=True, directed=False)
					flag_pos=0
					pos=False
					to_delete=[]
		random_rwr=random_rwr.T
		for i in enumerate(nodes):
			pvalues[i[1]]=np.add(pvalues[i[1]],np.greater(empirical_values[i[1]],random_rwr[i[0]]))

	'''
	Output.
 	'''
	#writing the results for global propagation
	f1=open(res_folder+str(damping)+"_rwr_scores.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+"\t".join(map(str,empirical_values[i]))+"\n")
	f1.close()

	f1=open(res_folder+str(damping)+"_pvalues.txt","w")
	for i in pvalues:
		f1.write(i+"\t"+"\t".join(map(str,pvalues[i]))+"\n")
	f1.close()

	f1=open(res_folder+str(damping)+"_start_seeds.txt","w")
	for i in seeds[0]+seeds[1]+seeds[2]:
		f1.write(i+"\t"+str(seeds_pos[i])+"\n")

	for i in seeds[3]+seeds[4]+seeds[5]:
		f1.write(i+"\t"+str(-seeds_neg[i])+"\n")
	f1.close()
 
	# # Returned variables.
	# return empirical_values,pvalues