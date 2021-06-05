# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 12:27:46 2020

@author: mkaanarici
"""

import os, networkx as nx, pandas as pd, numpy as np, sys,  ROC

#from multiprocessing import Process


class HeatDiffusion():
    def __init__(self,network,alpha=None,N=None,edge_weight=None):
        self.network=network 
        self.nodes=list(self.network.nodes())  
        self.dim=len(self.nodes)         
        self.node_thrshold=1/self.dim
        
        if N is None:
            N=3
        if alpha is None:
            alpha=0.8
        self.alpha=alpha
        self.N=N
        self.part1=None
        
        self.network_node_lngth=self.network.number_of_nodes()
        self.edge_weight=edge_weight
        #self.parameters_dict={}
        
        
        
        
    
    def set_parameters(self,alpha=None,N=None):
        if N is None:
            N=self.N
        if alpha is None:
            alpha=self.alpha
        self.__prepare_part1(alpha,N)
        
        
    def __prepare_part1(self,alpha=None,N=None):    
        if N is None:
            N=self.N
        if alpha is None:
            alpha=self.alpha
        if self.alpha==alpha and self.N==N:
            if self.part1:
                return self.part1
            
        print ("parameter matrix is being calculated")   
        degree_dict=dict(self.network.degree)
        degree=np.zeros((self.dim,self.dim))
        c=0
        for node in self.nodes:
            degree[c][c]=degree_dict[node]
            c+=1
        #print(degree)
        inver_degree=np.linalg.inv(degree)
        adjacency=nx.adjacency_matrix(self.network).toarray()
        L=np.matmul(inver_degree,adjacency)
        const=alpha/N
        identity=np.identity(self.dim)
        self.part1=np.linalg.matrix_power(identity-(const*L),N)
        return self.part1
    
    
    def load_initial_nodes(self, initial_nodes):
        weight=1/len(initial_nodes)
        self.initial_vector=np.zeros(self.dim)
        for i in initial_nodes:
            if i in self.nodes:
                index=self.nodes.index(i)
                self.initial_vector[index]=weight
        
        self.initial_nodes=[i for i in initial_nodes if i in self.nodes]
        weight=[1/(len(self.initial_nodes))]*len(self.initial_nodes)
        self.initial_weight=dict(zip(self.initial_nodes,weight))        
        return self.initial_vector
    
    def node_heat_diffusion(self, initial_nodes=None, alpha=None,N=None):  
        print("Node_heat_diffusion")
        if N != None or alpha != None:
            self.__prepare_part1(alpha=alpha,N=N)
        self.__prepare_part1(alpha=alpha,N=N)
        if initial_nodes:
            self.load_initial_nodes(initial_nodes)
        self.diffusued_weights=dict(zip(self.nodes,abs(self.initial_vector.dot(self.part1))))
        return self.diffusued_weights
    
    def edge_diffusion(self, initial_nodes=None, alpha=None,N=None,threshould=0.5): 
        
        print("edge_heat_diffusion")
        if N != None or alpha != None:
            self.__prepare_part1(alpha=alpha,N=N)
        node_diffusion_weights=self.node_heat_diffusion(initial_nodes=initial_nodes, alpha=alpha,N=N)
        heatdiffused=node_diffusion_weights
        returned_nodes=[i for i,j in heatdiffused.items() if j >=self.node_thrshold]
        Resolved_u=[]
        for u in heatdiffused:
            denom = self.network.degree(u)
            for v in self.network.neighbors(u): #  neighbors
                
                if (u not in returned_nodes) or (v not in returned_nodes):
                    self.network[u][v]['flux'] = 0
                    continue
                    
                if v in Resolved_u:
                    continue    
                    
                if self.edge_weight:                    
                    flux1 = heatdiffused[u]*self.network[u][v][self.edge_weight]/denom
                    flux2 = heatdiffused[v]*self.network[u][v][self.edge_weight]/denom
                    self.network[u][v]['flux']=min([flux1,flux2])

                else:
                    flux1=heatdiffused[u]*1/denom
                    flux2=heatdiffused[v]*1/denom                    
                    self.network[u][v]['flux'] = min([flux1,flux2])                    
            Resolved_u.append(u)
            
        tot_flux = 0
        for u,v in self.network.edges(): 
            if self.network[u][v]['flux'] == 0:
                self.network[u][v]['neglog_flux'] = sys.float_info.max ## maximum float value
            else:
                self.network[u][v]['neglog_flux'] = - (np.log10(self.network[u][v]['flux'])) 
                tot_flux+=self.network[u][v]['flux']
        
        Returned_edges = {} #dictionary of (u,v): -log(flux*flux)
        flux_sum = 0
        
        for u,v,d in sorted(self.network.edges(data=True), key=lambda t: t[2]['neglog_flux']):   
            if u==v:
                continue
            Returned_edges[(u,v)] = self.network[u][v]['neglog_flux']
            flux_sum+=self.network[u][v]['flux']
            if flux_sum/tot_flux > threshould:
                print('theshold of %f limits predictions to %d edges'  %(threshould,len(Returned_edges)))
                break        
        return Returned_edges 
        
    
    def heat_diffusion(self, alpha=None, N=None,threshould=0.20):     
        if alpha is None:alpha=self.alpha 
        if N is None :  N=self.N 
        heatdiffused=self.node_heat_diffusion(alpha=alpha,N=N)    
        return self.edge_diffusion(self, node_diffusion_weights=heatdiffused,threshould=threshould)

    
        
    
    def optimize_threshold(self,alphas, thresholds, pathway_edges):  
        
        pathway=nx.from_edgelist(pathway_edges)
        path_nodes=pathway.nodes()       
        
        output_df=pd.DataFrame(columns=["Alpha","Threshold","Edge_F1", "Edge_Recall", 
                                            "Edge_Precision", "Edge_FPR", "Edge_MCC", 
                                            "Edges_number", "Nodes_Number",
                                            "Node_F1","Node_Recall","Node_Precision",
                                            "Node_FPR", "MCC"])
        
        roc_node=ROC.ROC_node_based(self.network)
        
        
        for alpha in alphas: 
            
            heatdiffused=self.node_heat_diffusion(alpha=alpha,N=self.N)  
            returned_nodes=[i for i,j in heatdiffused.items() if j >= self.node_thrshold]
        
            roc_node.define_initial_nodes(self.initial_weight.keys())
            roc_node.define_pathways_nodes(path_nodes)
            roc_node.define_outcome_nodes(returned_nodes)
            roc_node.calculate_all_values()
            N_Recall=roc_node.calculate_recall()
            N_FPR=roc_node.calculate_false_positive_rate()
            N_Precision=roc_node.calculate_precision()
            N_F1=roc_node.calculate_F1_score()
            N_mcc=roc_node.calculate_MCC()
            
            
            Resolved_u=[]
            for u in heatdiffused:
                denom = self.network.degree(u)
                for v in self.network.neighbors(u): #  neighbors

                    if v in Resolved_u:
                        continue      
                        
                    if (u not in returned_nodes) or (v not in returned_nodes):
                        self.network[u][v]['flux'] = 0
                        continue
                        
                    if self.edge_weight:
                        flux1 = heatdiffused[u]*self.network[u][v][self.edge_weight]/denom
                        flux2 = heatdiffused[v]*self.network[u][v][self.edge_weight]/denom
                        self.network[u][v]['flux']=min([flux1,flux2])
                        
                    else:
                        flux1=heatdiffused[u]*1/denom
                        flux2=heatdiffused[v]*1/denom
                        self.network[u][v]['flux'] = min([flux1,flux2])

                Resolved_u.append(u)

            tot_flux = 0
            for u,v in self.network.edges(): 

                if self.network[u][v]['flux'] == 0 :
                    self.network[u][v]['neglog_flux'] = sys.float_info.max ## maximum float value
                else:
                    self.network[u][v]['neglog_flux'] = - (np.log10(self.network[u][v]['flux'])) 
                    tot_flux+=self.network[u][v]['flux']

            ROC_edge=ROC.ROC_edge_based(self.network)
            ROC_edge.define_pathway_edges(pathway_edges)

            

            for thres in thresholds: 
                
                Returned_edges = {} #dictionary of (u,v): -log(flux*flux)
                flux_sum = 0
                
                for u,v,d in sorted(self.network.edges(data=True), key=lambda t: t[2]['neglog_flux']):  

                    Returned_edges[(u,v)] = self.network[u][v]['neglog_flux']
                    flux_sum+=self.network[u][v]['flux']
                    if flux_sum/tot_flux > thres:
                        print('theshold of %f limits predictions to %d edges'  %(thres,len(Returned_edges)))
                        break
                    
                ROC_edge.define_outcome_edges(Returned_edges.keys())
                ROC_edge.calculate_all_values()
                Recall=ROC_edge.calculate_recall()
                FPR=ROC_edge.calculate_false_positive_rate()
                Precision=ROC_edge.calculate_precision()
                F1=ROC_edge.calculate_F1_score()
                mcc=ROC_edge.calculate_MCC()                
                
                output_df=output_df.append({"Alpha":alpha,"Threshold":thres, "Edge_F1":F1, "Edge_Recall":Recall, 
                                            "Edge_Precision":Precision, "Edge_FPR":FPR, "Edge_MCC":mcc, 
                                            "Edges_number":len(Returned_edges), "Nodes_Number":len(returned_nodes),
                                            "Node_F1":N_F1,"Node_Recall":N_Recall,"Node_Precision":N_Precision,
                                            "Node_FPR":N_FPR, "MCC":N_mcc},ignore_index=True)
                

                print ("Alpha:",alpha,"thres:",thres, "F1:",F1,"mcc:",mcc,"Recall:", Recall, "Precision:",Precision, 
                       "FPR:",FPR,"edge_len", len(Returned_edges))
        return output_df
        
       
    
# =============================================================================
# 
# =============================================================================

# int_name='HuRI'
# int_df=pd.read_csv(f'../Interactomes/{int_name}_with_UniProt.txt',sep=",")

# int_nx=nx.from_pandas_edgelist(int_df,"UniprotID_1","UniprotID_2","confidence")
# path_name='Wnt'
# with open (f'../Sample/{path_name}/chuk1_{path_name}_initial.nodes') as f:
#     nodes=f.readlines()[0].split()
#     #print (nodes)    
    
# with open (f'../Pathways/{path_name}-edges.txt') as f:
#     path_edges=[i.split("\t")[:2] for i in f.readlines()[1:]]
#     #print (len(path_edges))    
    
# hd=HeatDiffusion(int_nx)

# hd.load_initial_nodes(nodes)

# Opt_Res_Df=hd.optimize_threshold(alphas=[0.4,0.8],thresholds=[0.2,0.8],pathway_edges=path_edges)

# Opt_Res_Df.to_csv("optiziationtrial.tsv")

   
    