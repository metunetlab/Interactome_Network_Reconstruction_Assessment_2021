
"""
Created on Sat Nov  7 00:02:24 2020

@author: mkaanarici
"""

import os, pandas as pd, networkx as nx, numpy as np
from multiprocessing import Process
import ROC, RWR_F



Tunning_paths=["Wnt","TGFbetaReceptor","TNFalpha","TCR"]

Interactomelist=[i for i in os.listdir("../../Interactomes/") if i.endswith("with_UniProt.txt") ]
print (Interactomelist)

path_list=[i for i in os.listdir("../Sample") if i not in Tunning_paths ]
print (path_list,"\n"*3)
#
#
#
#
alphas=list(np.arange(0.0,1,0.05))

thresholds=list(np.arange(0.05,1.05,0.05))


# =============================================================================
# Optimization
# =============================================================================

#def opt_helper(rwalk,initial_nodes,file_name):    
#      rwalk.load_initial_nodes(initial_nodes)
#      output_df=rwalk.optimize_threshold(alphas=alphas,thresholds=thresholds,pathway_edges=path_edges)
#      output_df.to_csv(file_name,sep="\t",index_label=False,index=False)
#
#
#for interactome in Interactomelist:
#     int_name=interactome.split("_")[0]
#     print (int_name)
#    
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
#     rwalk=RWR_F.RandomWalk(network)
#     procs = []
#    
#    
#     for path_name in Tunning_paths:
#         os.makedirs(f'../Tunning/{int_name}/{path_name}',exist_ok=True)
#         with open (f'../../NetPath/{path_name}-edges.txt') as f:
#             path_edges=[i.split("\t")[:2] for i in f.readlines()[1:]]
#        
#        
#        
#         for i in range(1,6):
#             with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes') as f:
#                 initial_nodes=f.readlines()[0].split()
#                
#                
#             file_name=f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt'
#            
#             proc = Process(target=opt_helper,args=(rwalk,initial_nodes,file_name))
#             procs.append(proc)
#             proc.start()
#            
#             # rwalk.load_initial_nodes(initial_nodes)
#             # output_df=rwalk.optimize_threshold(alphas=alphas,thresholds=thresholds,pathway_edges=path_edges)
#             # output_df.to_csv(file_name)
#            
#     for proc in procs:
#             # print ('ending')
#         proc.join()
            
        
        
# =============================================================================
# Get optimum alpha and threshold and run y
# =============================================================================

#def help_rw_run(network,int_name,path_name):
#     rw=RWR_F.RandomWalk(network)
#            
#     os.makedirs(f'../Results/{int_name}/{path_name}',exist_ok=True)
#    
#     print (int_name, path_name)    
#     for i in range(1,6):
#         with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes') as f:
#             nodes=f.readlines()[0].split()
#         # print (path_name,nodes)
#         init=rw.load_initial_nodes(nodes)
#         return_edges=rw.random_walk(alpha = Validation_res[int_name]['alpha'], threshould=Validation_res[int_name]['threshould'])
#         # print (return_edges)
#         with open (f'../Results/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.res',"w") as f:
#             for prot1,prot2 in return_edges.keys():
#                 f.writelines(f'{prot1}\t{prot2}\t{return_edges[(prot1,prot2)]}\n')
#                
#
#
#
#Validation_res={}
#for interactome in Interactomelist:
#     int_name=interactome.split("_")[0]
#    
#     Best_Res_all=pd.DataFrame(columns=['Pathway','Alpha', 'Threshold', 'Edge_F1', 'Edge_Recall', 'Edge_Precision',
#             'Edge_FPR', 'Edge_MCC', 'Edges_number', 'Nodes_Number', 'Node_F1',
#             'Node_Recall', 'Node_Precision', 'Node_FPR', 'MCC'])
#    
#     for path_name in Tunning_paths:
#        
#        
#         for i in range(1,6):
#             Opt_Res=pd.read_csv(f'../Tunning/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.opt',sep="\t")
#             Opt_Res['Pathway']=path_name
#             Best_Res_all=Best_Res_all.append(Opt_Res[Opt_Res.Edge_F1.isin(Opt_Res.Edge_F1.nlargest(10))],ignore_index=True)
#         Alpha=sum(Best_Res_all.Alpha)/len(Best_Res_all.index)
#         Thres=sum(Best_Res_all.Threshold)/len(Best_Res_all.index)
#     Validation_res[int_name]={"alpha":round(Alpha,2),"threshould":round(Thres,2)}
#    
#
#
#with open ("../Tunning/RWF_optimized.parameters","w") as f:
#     f.writelines("interactome\talpha\tthreshold\n")
#     for int_name in Validation_res.keys():
#         print (int_name, Validation_res[int_name]["alpha"])
#         f.writelines(f'{int_name}\t{Validation_res[int_name]["alpha"]}\t{Validation_res[int_name]["threshould"]}\n')    
#
#print (Validation_res)
#
#
#for interactome in Interactomelist:
#     int_name=interactome.split("_")[0]
#     print ("--------"*2,int_name,"--------"*2)    
#     os.makedirs(f'../Results/{int_name}',exist_ok=True)    
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
#    
#    
#     procs=[]
#     for path_name in path_list:
#         proc=Process(target=help_rw_run,args=(network,int_name,path_name))
#         proc.start()
#         procs.append(proc)
#    
#     for proc in procs:
#         proc.join()
#        
#        
#        
#         rw=RWR_F.RandomWalk(network)
#        
#         os.makedirs(f'../Results/{int_name}/{path_name}',exist_ok=True)
#        
#         print (int_name, path_name)    
#         for i in range(1,6):
#             with open (f'../Sample/{path_name}/chuk{i}_{path_name}_initial.nodes') as f:
#                 nodes=f.readlines()[0].split()
#             # print (path_name,nodes)
#             init=rw.load_initial_nodes(nodes)
#             return_edges=rw.random_walk(alpha = Validation_res[int_name]['alpha'], threshould=Validation_res[int_name]['threshould'])
#             # print (return_edges)
#             with open (f'../Results/{int_name}/{path_name}/{int_name}_{path_name}_chunk{i}.res',"w") as f:
#                 for prot1,prot2 in return_edges.keys():
#                     f.writelines(f'{prot1}\t{prot2}\t{return_edges[(prot1,prot2)]}\n')
                
                                        


# =============================================================================
# performance
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
            
            print (roc_node.calculate_recall())
                
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
                    










            
            
            
            
            
            
            
            