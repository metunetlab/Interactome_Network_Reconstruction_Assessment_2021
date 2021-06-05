# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 13:30:37 2020

@author: mkaanarici
"""

import os, pandas as pd, networkx as nx, numpy as np
import ROC, HeatDiffusion_F

from multiprocessing import Process,Pool

Tunning_paths=["Wnt","TGFbetaReceptor","TNFalpha","TCR"]

Interactomelist=[i for i in os.listdir("../../Interactomes/") if i.endswith("with_UniProt.txt") ]

print (Interactomelist)

path_list=[i for i in os.listdir("../Sample") if i not in Tunning_paths ]
print (path_list,"\n"*3)

alphas=list(np.arange(0.0,1,0.05))

thresholds=list(np.arange(0.05,1,0.05))



# =============================================================================
# Parameter Tunning
# =============================================================================



#
#def opt_helper(network,hd,alpha,path_name,int_name):  
#
#        
#     roc_node=ROC.ROC_node_based(network)
#     roc_edge=ROC.ROC_edge_based(network)  
#    
#     os.makedirs(f'../Tunning/{int_name}/{path_name}',exist_ok=True)
#     with open (f'../../NetPath/{path_name}-edges.txt') as f:
#         path_edges=[i.split("\t")[:2] for i in f.readlines()[1:]]
#          
#     pathway=nx.from_edgelist(path_edges)
#     path_nodes=pathway.nodes()       
#    
#     roc_edge.define_pathway_edges(path_edges)
#    
#     for i in range(1,6):
#         print (int_name, " - ",alpha," - ",path_name," - chunk",i )
#         with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes') as f:
#             initial_nodes=f.readlines()[0].split()
#         hd.load_initial_nodes(initial_nodes)
#         returned_dict=hd.node_heat_diffusion(alpha=alpha)
#        
#        
#         roc_node.define_initial_nodes(initial_nodes)
#         roc_node.define_pathways_nodes(path_nodes)
#        
#         for threshould in thresholds:
#             returned_edges=hd.edge_diffusion(node_diffusion_weights=returned_dict,threshould=threshould)   
#             returned_nodes=[]
#             for k in returned_edges:
#                 if k[0] not in returned_nodes:
#                     returned_nodes.append(k[0])
#                 if k[1] not in returned_nodes:
#                     returned_nodes.append(k[1])
#            
#            
#             roc_node.define_outcome_nodes(returned_dict.keys())
#             roc_node.calculate_all_values()
#             N_Recall=roc_node.calculate_recall()
#             N_FPR=roc_node.calculate_false_positive_rate()
#             N_Precision=roc_node.calculate_precision()
#             N_F1=roc_node.calculate_F1_score()
#             N_mcc=roc_node.calculate_MCC()
#            
#             roc_edge.define_outcome_edges(returned_edges.keys())
#            
#             roc_edge.calculate_all_values()
#             Recall=roc_edge.calculate_recall()
#             FPR=roc_edge.calculate_false_positive_rate()
#             Precision=roc_edge.calculate_precision()
#             F1=roc_edge.calculate_F1_score()
#             mcc=roc_edge.calculate_MCC()  
#            
#             with open (f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt',"a") as f:
#                 f.writelines(f'{alpha}\t{threshould}\t{F1}\t{Recall}\t{Precision}\t{FPR}\t{mcc}\t{len(returned_edges)}\t{len(returned_nodes)}\t{N_F1},{N_Recall}\t{N_Precision}\t{N_FPR}\t{N_mcc}\n')
#                        
#print ('multiprocessing')
#
#for interactome in Interactomelist:
#     int_name=interactome.split("_")[0]
#     os.makedirs(f'../Tunning/{int_name}',exist_ok=True)
#     for path_name in Tunning_paths:
#         os.makedirs(f'../Tunning/{int_name}/{path_name}',exist_ok=True)
#         for i in range(1,6):
#             with open (f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt',"w") as f:
#                 f.writelines("Alpha\tThreshold\tEdge_F1\tEdge_Recall\tEdge_Precision\tEdge_FPR\tEdge_MCC\tEdges_number\tNodes_Number\tNode_F1\tNode_Recall\tNode_Precision\tNode_FPR\tMCC\n")
#     
#
#
#for interactome in Interactomelist:
#     int_name=interactome.split("_")[0]
#     print ("--------"*2,int_name,"--------"*2) 
#     os.makedirs(f'../Tunning/{int_name}',exist_ok=True)
#    
#     interactome_df=pd.read_csv(f'../../Interactomes/{interactome}')
#     col=interactome_df.columns
#     if len(col)==3:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[2])
#     elif len(col)==2:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1])
#    
#     elif(len(col))==7:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[6])
#        
#    
#    
#     hd=HeatDiffusion_F.HeatDiffusion(network)    
#     # roc_node=ROC.ROC_node_based(network)
#     # roc_edge=ROC.ROC_edge_based(network)
#    
#    
#    
#     for alpha in alphas:    
#        
#                  
#         hd.set_parameters(alpha=alpha)  
#              
#         procs=[]
#         for path_name in Tunning_paths:
#             proc=Process(target=opt_helper,args=(network,hd,alpha,path_name,int_name))
#             proc.start()
#             procs.append(proc)
#        
#         for proc in procs:
#             proc.join()
            
            
                
# =============================================================================
# PathwayReconstruction with HDF and tunned parameters
# =============================================================================
            
#Tunning_res={}
#for interactome in Interactomelist:
#    int_name=interactome.split("_")[0]
#    
#    Best_Res_all=pd.DataFrame(columns=['Pathway','Alpha', 'Threshold', 'Edge_F1', 'Edge_Recall', 'Edge_Precision',
#            'Edge_FPR', 'Edge_MCC', 'Edges_number', 'Nodes_Number', 'Node_F1',
#            'Node_Recall', 'Node_Precision', 'Node_FPR', 'MCC'])
#    
#    for path_name in Tunning_paths:
#        
#        
#        for i in range(1,6):
#            print(f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt')
#            Opt_Res=pd.read_csv(f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt',sep="\t")
#            Opt_Res['Pathway']=path_name
#            print (Opt_Res)
#            Best_Res_all=Best_Res_all.append(Opt_Res[Opt_Res.Edge_F1.isin(Opt_Res.Edge_F1.nlargest(10))],ignore_index=True)
#        Alpha=sum(Best_Res_all.Alpha)/len(Best_Res_all.index)
#        Thres=sum(Best_Res_all.Threshold)/len(Best_Res_all.index)
#    Tunning_res[int_name]={"alpha":round(Alpha,2),"threshould":round(Thres,2)}
#
#with open ("../Tunning/HDF_optimized.parameters","w") as f:
#    f.writelines("interactome\talpha\tthreshold\n")
#    for int_name in Tunning_res.keys():
#        print (int_name, Tunning_res[int_name]["alpha"])
#        f.writelines(f'{int_name}\t{Tunning_res[int_name]["alpha"]}\t{Tunning_res[int_name]["threshould"]}\n')    
#
#print (Tunning_res)


# =============================================================================
# Pathway Reconstruction with HDF
# =============================================================================


#for interactome in Interactomelist:  
#     int_name=interactome.split("_")[0]
#     print ("--------"*2,int_name,"--------"*2)    
#     os.makedirs(f'../Results/{int_name}',exist_ok=True)    
#     interactome_df=pd.read_csv(f'../Interactomes/{interactome}')
#     col=interactome_df.columns
#     if len(col)==3:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[2])
#     elif len(col)==2:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1])
#    
#     elif(len(col))==7:
#         network=nx.from_pandas_edgelist(interactome_df,col[0],col[1],col[6])  
#        
#     hd=HeatDiffusion_F.HeatDiffusion(network)
#     hd.set_parameters(alpha=Tunning_res[int_name]["alpha"])
#    
#     for path_name in path_list:
#         for i in range(1,6):
#             os.makedirs(f'../Results/{int_name}/{path_name}',exist_ok=True)  
#             with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes') as f:
#                 initial_nodes=f.readlines()[0].split()
#             hd.load_initial_nodes(initial_nodes)
#             #returned_dict=hd.node_heat_diffusion()
#             returned_edges=hd.edge_diffusion(threshould=Tunning_res[int_name]["threshould"])
#             with open (f'../Results/{int_name}/{path_name}/{int_name}_{int_name}_chunk{i}.res',"w") as f:
#                 for key in returned_edges.keys():
#                     f.writelines(f'{key[0]}\t{key[1]}\t{returned_edges[key]}\n')
    

# =============================================================================
# performance
# =============================================================================

for interactome in Interactomelist:
     int_name=interactome.split("_")[0]
     print (int_name)
    
     os.makedirs('../Results/00_performance/',exist_ok=True)
    
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
                    


