# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 23:18:12 2020

@author: mkaanarici
"""

import os , pandas as pd, networkx as nx, numpy as np,ROC


from ShortestPath import Shortest_Paths as sp

from multiprocessing import Process

Interactomelist=[i for i in os.listdir("../../Interactomes/") if i.endswith("with_UniProt.txt")]
#Interactomelist=['OmniPath_with_UniProt.txt']
print (Interactomelist)


path_list=[i for i in os.listdir("../../NetPath")  if i.endswith("-edges.txt")]
print (path_list)
print ("************"*5)

# =============================================================================
# shortest path
# =============================================================================


#for interactome in Interactomelist:
#     name=interactome.split("_")[0]
#     print (name)
#    
#     interactome_df=pd.read_csv(f'../../Interactomes/{interactome}')
#     print (interactome_df)
#    
#     print (f'{name} is loading into networkx')
#     network=nx.from_pandas_edgelist(interactome_df, "UniprotID_1", "UniprotID_2",create_using=nx.Graph)
##    
#     if os.path.isdir(f'../Results/{name}') is False:
#         os.makedirs(f'../Results/{name}')
##    
#     shrt_path=sp(network)
#     procs = [] 
#     for path_name in path_list:
#         path_name=path_name.split("-")[0]
#         if os.path.isdir(f'../Results/{name}/{path_name}') is False:
#             os.makedirs(f'../Results/{name}/{path_name}')
#         print (path_name,"  -  ", name)
#         for i in range(1,6):
#             with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes',"r") as f :
#                 nodes=f.readline().split("\t")
#                
#                
#                 proc = Process(target=shrt_path.write_shortest_path_edges_among,args=( nodes, nodes,f'../Results/{name}/{path_name}/{name}_{path_name}_chuk{i}.res'))
#                 procs.append(proc)
#                 proc.start()
#         for proc in procs:
#             proc.join()
            
            
            


            
        
        
    



# =============================================================================
#   Performance     
# =============================================================================

for interactome in Interactomelist:
    int_name=interactome.split("_")[0]
    print (int_name)
    
    
    interactome_df=pd.read_csv(f'../../Interactomes/{interactome}')
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
        path_name=path_name.split("-")[0]
        print (path_name)
        pathway_pd=pd.read_csv(f'../../NetPath/{path_name}-edges.txt',sep="\t")        
        pathway_pd=pd.read_csv(f'../../NetPath/{path_name}-edges.txt',sep="\t")[["#tail","head"]]
        pathway_edges=set(zip(list(pathway_pd['#tail']),list(pathway_pd['head'])))
        
        pathway=nx.from_edgelist(pathway_edges)
        path_nodes=list(pathway.nodes())
        
        E_temp_mcc,E_temp_recall,E_temp_fpr,E_temp_precision,E_temp_F1=[],[],[],[],[]
        N_temp_mcc,N_temp_recall,N_temp_fpr,N_temp_precision,N_temp_F1=[],[],[],[],[]
        for i in range (1,6):            
            with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes','r') as f:
                initial_nodes=f.readline().split("\t")
        
            Outcome_pd=pd.read_csv(f'../Results/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.res',sep="\t",names=["node1","node2","scores"])
        

            outcome_nx=nx.from_pandas_edgelist(Outcome_pd,"node1","node2")
            returned_nodes=outcome_nx.nodes()
            returned_edges=outcome_nx.edges()
        
            roc_node.define_initial_nodes(initial_nodes)            
            roc_node.define_pathways_nodes(path_nodes)
            
            roc_node.define_outcome_nodes(returned_nodes)
            roc_node.calculate_all_values()
            
            N_temp_recall.append(roc_node.calculate_recall())
            N_temp_fpr.append(roc_node.calculate_false_positive_rate())
            N_temp_precision.append(roc_node.calculate_precision())
            N_temp_F1.append(roc_node.calculate_F1_score())
            N_temp_mcc.append(roc_node.calculate_MCC())
                
                
            roc_edge.define_pathway_edges(pathway_edges)
        
            roc_edge.define_pathway_edges(pathway_edges)
            roc_edge.define_outcome_edges(returned_edges)
            roc_edge.calculate_all_values()
        
            E_temp_recall.append(roc_edge.calculate_recall())
            E_temp_fpr.append(roc_edge.calculate_false_positive_rate())
            E_temp_precision.append(roc_edge.calculate_precision())
            E_temp_F1.append(roc_edge.calculate_F1_score())
            E_temp_mcc.append(roc_edge.calculate_MCC())
        with open (f'../Results/00_performance/{int_name}.performance',"a") as f:
            f.writelines(f'{path_name}\t{len(pathway_edges)}\t{round(np.mean(E_temp_F1),2)}\t{round(np.mean(E_temp_mcc),2)}\t{round(np.mean(E_temp_recall),2)}\t{round(np.mean(E_temp_precision),2)}\t{round(np.mean(E_temp_fpr),2)}\t{len(path_nodes)}\t{round(np.mean(N_temp_F1),2)}\t{round(np.mean(N_temp_mcc),2)}\t{round(np.mean(N_temp_recall),2)}\t{round(np.mean(N_temp_precision),2)}\t{round(np.mean(N_temp_fpr),2)}\n')
                    









