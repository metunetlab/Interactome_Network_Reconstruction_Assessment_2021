# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:52:55 2020

@author: mkaanarici
"""

import os, numpy as np,gc, pandas as pd, networkx as nx, ROC
from PCSF import PCSF as pcsf

gc.set_threshold(700,10,10)


Tunning_paths=["Wnt","TGFbetaReceptor","TNFalpha","TCR"]

Interactomelist=[i for i in os.listdir("../Interactomes/") if i.endswith("oi2") ]
print (Interactomelist)

path_list=[i for i in os.listdir("../Sample") if i not in Tunning_paths ]
print (path_list,"\n"*3)




# =============================================================================
# Pamaeter tunning over Wnt, TGFBR, TNF alpha and TCR
# # =============================================================================


#Ws = list(np.arange(0,5,0.5))
#Bs = list(np.arange(0,5,0.5))
#Gs = list(np.arange(0.5,10,1))
#
#
#
#
#validation_path_prizes=[]
#for path in Tunning_paths:
#    for i in range (1,6):
#        validation_path_prizes.append(f'../Sample/{path}/{path}_chunk{i}_prize.tsv')
#int_name='OmniPath'
#interactome_file=f'../Interactomes/{int_name}.oi2'
#
#for interactome_file in Interactomelist:
#     int_name=interactome_file.split(".")[0]
#     print ("**********",int_name,"************")
#     os.makedirs(f'../Tunning/{int_name}', exist_ok=True )
#     for path in Tunning_paths:
#         # if int_name=='Hippie' and path!='TCR':
#             # continue
#
#         print(int_name," - ",path)
#         os.makedirs(f'../Tunning/{int_name}/{path}', exist_ok=True )
#         for i in range (1,6):
#             print(int_name, path,"chunk",i)
#             if  (int_name=='STRING' and path in ["Wnt",	"TGFbetaReceptor","TNFalpha"]) or (path=="TCR" and  i in [1,2]):
#                 continue
#          
#             OI2=pcsf(f'../Interactomes/{interactome_file}')
#             OI2.get_optimim_parameters(prize_file=f'../Sample/{path}/{path}_chunk{i}_prize.tsv',Ws=Ws,Bs=Bs,Gs=Gs)
#             OI2.write_optimum_parameters(f'../Tunning/{int_name}/{path}/{int_name}_{path}_chunk{i}')
#             OI2=None



# =============================================================================
# getting validation parameters for each interactome"
# =============================================================================


#parameter_pool=[]
#
#for interactome in Interactomelist:
#     int_name=interactome.split(".")[0] 
#
#     for path in Tunning_paths:
#         for i in range(1,6):
#             with open (f'../Tunning/{int_name}/{path}/{int_name}_{path}_chunk{i}.parameters',"r") as f:
#                 parameters={i.strip() for i in f.readlines()}
#                 parameter_pool.append(parameters)
#     parameter_intersection=parameter_pool[0]
#     for parameters in parameter_pool:
#         parameter_intersection=set.intersection(parameter_intersection,parameters)
#     print (int_name,"  -  ",len(parameter_intersection))
#     with open (f'../Tunning/00params/{int_name}.parameters',"w") as f:
#         for parameter in parameter_intersection:
#             f.writelines(parameter+"\n")

 

# =============================================================================
# 
# =============================================================================


#for interactome in Interactomelist:    
#    int_name=interactome.split(".")[0] 
#    print ("**"*5," - ", int_name," - ","**"*5)
#    parameter_file=f'../Tunning/00params/{int_name}.parameters'
#    oi2=pcsf(f'../Interactomes/{interactome}')
#    os.makedirs(f'../Results/{int_name}/', exist_ok=True)
#    for path in path_list:
#        os.makedirs(f'../Results/{int_name}/{path}', exist_ok=True)
#        for i in range(1,2):
#            
#            output_name=f'../Results/{int_name}/{path}/{int_name}_{path}_chunk{i}'            
#            print (int_name," - ",path, "chunk ",i)
#            
#            prize_file=f'../Sample/{path}/{path}_chunk{i}_prize.tsv'
#            oi2.write_augmented_networks_intersection(generated_network_name=output_name,
#                                                      prize_file=prize_file,
#                                                      parameter_file=parameter_file)
#            
# =============================================================================
# Performance
# =============================================================================


for interactome in Interactomelist:
    int_name=interactome.split(".")[0]
    print (int_name)
#    
    os.makedirs('../Results/00_performance/',exist_ok=True)
    
    interactome_df=pd.read_csv(f'../Interactomes/{interactome}',sep="\t")
    col=interactome_df.columns
    if len(col)==3:
        network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[2])
    elif len(col)==2:
        network=nx.from_pandas_edgelist(interactome_df,col[0],col[1])
    
    elif(len(col))==7:
        network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[6])
    
    roc_edge=ROC.ROC_edge_based(network)
    roc_node=ROC.ROC_node_based(network)
    
    with open (f'../Results/00_performance/{int_name}.performance',"w") as f:
        f.writelines("Pathway\tedge_len\tE_F1\tE_mcc\tE_Recall\tE_Precision\tE_FPR\tlen_Node\tN_F1\tN_mcc\tN_Recall\tN_Precision\tN_FPR\n")
    for path_name in path_list:
        print (path_name)
#        
        pathway_pd=pd.read_csv(f'../../Netpath/{path_name}-edges.txt',sep="\t")        
        pathway_pd=pd.read_csv(f'../../Netpath/{path_name}-edges.txt',sep="\t")[["#tail","head"]]
        pathway_edges=set(zip(list(pathway_pd['#tail']),list(pathway_pd['head'])))
        
        pathway=nx.from_edgelist(pathway_edges)
        path_nodes=list(pathway.nodes())
        
        E_temp_mcc,E_temp_recall,E_temp_fpr,E_temp_preision,E_temp_F1=[],[],[],[],[]
        N_temp_mcc,N_temp_recall,N_temp_fpr,N_temp_precision,N_temp_F1=[],[],[],[],[]
        for i in range (1,6):
            if os.path.isfile(f'../Results/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.res') is False:
                continue
            with open (f'../Sample/{path_name}/{path_name}_chunk{i}_prize.tsv','r') as f:
                initial_nodes=[i.split("\t")[0] for i in f.readlines()[1:]]
        
            Outcome_pd=pd.read_csv(f'../Results/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.res',sep="\t",names=["node1","node2","scores"])
#            print (Outcome_pd.head())
#        
#
            outcome_nx=nx.from_pandas_edgelist(Outcome_pd,source="node1",target="node2")
            returned_nodes=outcome_nx.nodes()
            returned_edges=outcome_nx.edges()
            
#        
            roc_node.define_initial_nodes(initial_nodes)            
            roc_node.define_pathways_nodes(path_nodes)
#            
            roc_node.define_outcome_nodes(returned_nodes)
            roc_node.calculate_all_values()
#            print (roc_node.calculate_recall())
##            break
#            
            N_temp_recall.append(roc_node.calculate_recall())
            N_temp_fpr.append(roc_node.calculate_false_positive_rate())
            N_temp_precision.append(roc_node.calculate_precision())
            N_temp_F1.append(roc_node.calculate_F1_score())
            N_temp_mcc.append(roc_node.calculate_MCC())
   
        
            roc_edge.define_pathway_edges(pathway_edges)
        
            roc_edge.define_outcome_edges(returned_edges)
            roc_edge.calculate_all_values()
        
            E_temp_recall.append(roc_edge.calculate_recall())
            E_temp_fpr.append(roc_edge.calculate_false_positive_rate())
            E_temp_precision.append(roc_edge.calculate_precision())
            E_temp_F1.append(roc_edge.calculate_F1_score())
            E_temp_mcc.append(roc_edge.calculate_MCC())
        print (f'{path_name}\t{len(pathway_edges)}\t{round(np.mean(E_temp_F1),2)}\t{round(np.mean(E_temp_mcc),2)}\t{round(np.mean(E_temp_recall),2)}\t{round(np.mean(E_temp_precision),2)}\t{round(np.mean(E_temp_fpr),2)}\t{len(path_nodes)}\t{round(np.mean(N_temp_F1),2)}\t{round(np.mean(N_temp_mcc),2)}\t{round(np.mean(N_temp_recall),2)}\t{round(np.mean(N_temp_precision),2)}\t{round(np.mean(N_temp_fpr),2)}\n')
        with open (f'../Results/00_performance/{int_name}.performance',"a") as f:
            f.writelines(f'{path_name}\t{len(pathway_edges)}\t{round(np.mean(E_temp_F1),2)}\t{round(np.mean(E_temp_mcc),2)}\t{round(np.mean(E_temp_recall),2)}\t{round(np.mean(E_temp_precision),2)}\t{round(np.mean(E_temp_fpr),2)}\t{len(path_nodes)}\t{round(np.mean(N_temp_F1),2)}\t{round(np.mean(N_temp_mcc),2)}\t{round(np.mean(N_temp_recall),2)}\t{round(np.mean(N_temp_precision),2)}\t{round(np.mean(N_temp_fpr),2)}\n')

             
