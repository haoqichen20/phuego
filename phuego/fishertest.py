# -*- coding: utf-8 -*-

from scipy.stats import fisher_exact
import sys
def fisher_test(protein_list,threshold,component,path_def,starting_proteins,
                uniprot_to_gene,geneset_path, fname):
    protein=protein_list
    temp={}
    for ii in component:
        descr={}
        temp[ii]=[]
        f1=open(geneset_path+ii+"_descr.txt","r")
        seq=f1.readline()
        while (seq!=""):
            seq=seq.strip().split("\t")
            descr[seq[0]]=seq[1]
            seq=f1.readline()
        f1.close()
        f1=open(geneset_path+ii+".txt","r")
        seq=f1.readline()
        fisher={}
        fisher_count={}
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher[seq[0]]=seq[1:]
            seq=f1.readline()
        f1.close()
        f1=open(geneset_path+ii+"_count.txt","r")
        seq=f1.readline()
        while(seq!=""):
            seq=seq.strip().split("\t")
            fisher_count[seq[0]]=int(seq[1])
            seq=f1.readline()
        f1.close()
        fisherset={}
        fisherset_count=[]
        protein_annotation={}
        for i in protein:
            if i in fisher:
                for j in fisher[i]:
                    if j in protein_annotation:
                        protein_annotation[j].append(i)
                    else:
                        protein_annotation[j]=[]
                        protein_annotation[j].append(i)

                    if j in fisherset:
                        fisherset[j]=fisherset[j]+1
                    else:
                        fisherset[j]=1

            else:
                fisherset_count.append(i)


        totalfisher=len(fisher)
        numberofproteins=len(protein)-len(list(set(fisherset_count)))
        fisher={}
        fisher_value_ic={}
        fisher_value={}
        fisher_value_no={}
        lenfisherset=len(fisherset)
        for i in fisherset:

            a=fisherset[i]
            b=numberofproteins-a
            c=fisher_count[i]-a
            d=totalfisher-a-b-c
            table=[[a,b],[c,d]]
            fisher[i]=fisher_exact(table,alternative ="greater")[1]
            if fisher[i]<(threshold/lenfisherset):
                if fisher[i] in fisher_value:
                    #fisher_value_ic[ic[i]].append(i)
                    fisher_value[fisher[i]].append(i)
                else:
                    #fisher_value_ic[ic[i]]=[]
                    fisher_value[fisher[i]]=[]
                    #fisher_value_ic[ic[i]].append(i)
                    fisher_value[fisher[i]].append(i)
            """
            else:
                if fisher[i] in fisher_value_no:
                    fisher_value_no[fisher[i]].append(i)
                else:
                    fisher_value_no[fisher[i]]=[]
                    fisher_value_no[fisher[i]].append(i)
            """
        
        # Allow input name.
        f2=open(path_def+fname+"_"+ii+"fisher.txt","w")
        
        # Add a header to the fisher output.
        f2.write("Geneset"+"\t"+"adjusted_p"+"\t"+"Geneset_size"+"\t"+"N_nodes"+"\t"
        +"Description"+"\t"+"Nodes"+"\n")
        
        for i in sorted(fisher_value):
            for j in fisher_value[i]:
                temp_gene=[]
                summation=0.0
                for k in list(set(protein_list).intersection(set(protein_annotation[j]))):
                    #summation+=pagerank.get(k,0.0)
                    temp_gene.append(uniprot_to_gene.get(k,k))
                #summation=summation/float(len(temp_gene))
                f2.write(j+"\t"+str(i)+"\t"+str(fisherset[j])+"\t"+str(len(list(set(starting_proteins).intersection(set(protein_annotation[j])))))+"\t"+descr[j]+"\t"+"\t".join(temp_gene)+"\n")
"""
        for i in sorted(fisher_value_no):
            for j in fisher_value_no[i]:
                temp_gene=[]
                summation=0.0
                for k in list(set(protein_list).intersection(set(protein_annotation[j]))):
                    temp_gene.append(uniprot_to_gene.get(k,k))
                    summation+=pagerank[k]
                summation=summation/float(len(temp_gene))
                f2.write("*"+"\t"+str(summation)+"\t"+j+"\t"+str(i)+"\t"+str(fisherset[j])+"\t"+str(len(list(set(starting_proteins).intersection(set(protein_annotation[j])))))+"\t"+descr[j]+"\t"+"\t".join(temp_gene)+"\n")

"""
