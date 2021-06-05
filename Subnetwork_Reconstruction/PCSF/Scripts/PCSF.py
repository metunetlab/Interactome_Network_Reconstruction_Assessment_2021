# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 13:47:21 2020

@author: mkaanarici
"""

import numpy as np
import pandas as pd
import networkx as nx,gc
import OmicsIntegrator as oi
import multiprocessing

gc.set_threshold(1,2,2)
class PCSF:
    def __init__(self, network_file, parameters=None):
        self.network_file=network_file
        if parameters==None:
            self.OI2_graph=oi.Graph(self.network_file)
        if parameters:
            self.OI2_graph=oi.Graph(self.network_file,parameters)
        self.augmented_pool,self.forest_pool=[],[]
    def load_initial_nodes(self, prize_file):
        self.prize_file=prize_file
        self.OI2_graph.prepare_prizes(prize_file)
        
    def get_optimim_parameters(self, prize_file, Ws, Bs, Gs,parameters=None):
        if parameters is None:
            parameters = {"noise": 0.1, 
                          "dummy_mode": "terminals", 
                          "exclude_terminals": False, 
                          "seed": 1}
        #int_oi2=oi.Graph(interactome_file,parameters)
        results = self.OI2_graph.grid_search(prize_file, Ws, Bs, Gs)
        membership_df = oi.summarize_grid_search(results, "membership")

        prize=pd.read_csv(prize_file,sep="\t")
        initial_nodes=list(prize.name)
                
        results_with_terminals=membership_df[membership_df.index.isin(initial_nodes)]
        Initial_node_covers=results_with_terminals.sum().sort_values(ascending=False).to_frame(name="Covering_nodes")
        self.parameters=set(Initial_node_covers[Initial_node_covers["Covering_nodes"]==max(Initial_node_covers["Covering_nodes"])].index)
        del results,membership_df,prize,Initial_node_covers
        return self.parameters
    def get_optimim_parameters_among(self, prize_files, Ws, Bs, Gs,parameters=None):
        self.parameters=set()
        for prize_file in prize_files:
            temp_parameters=self.get_optimim_parameters(prize_file, Ws, Bs, Gs,parameters)
            self.parameters=self.parameters.intersection(temp_parameters)
        #self.parameters=set(self.parameters)
        return self.parameters
    
    def write_optimum_parameters(self,f_name=None):
        if f_name is None:
            f_name="parameters"
        with open (f'{f_name}.parameters',"w") as f:
            for param in self.parameters:
                f.writelines(param+"\n")
    def load_parameter_files(self, parameter_file=None):
        if parameter_file is None:
            parameter_file="parameters.parameters"
        with open (f'{parameter_file}',"r") as f:
            self.parameters={i.strip() for i in f.readlines()}
        return self.parameters
    
    def _get_network(self,parameter=None):
        if parameter!=None:
            graph = oi.Graph(self.network_file, parameter)
            graph.prepare_prizes(parameter["prize_file"])
            vertex_indices, edge_indices = graph.pcsf()
            forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
            #self.augmented_pool.append(augmented_forest)
            #self.forest_pool.append(forest)
            return{"forest":forest,"augmented_forest":augmented_forest}
        else:
            parameter={i:float(j) for i,j in parameter.items() if i=='w' and j!="0.00"}
            graph = oi.Graph(self.network_file, parameter)
            graph.prepare_prizes(parameter["prize_file"])
            vertex_indices, edge_indices = graph.pcsf()
            forest, augmented_forest = graph.output_forest_as_networkx(vertex_indices, edge_indices)
            #self.augmented_pool.append(augmented_forest)
            #self.forest_pool.append(forest)
            return{"forest":forest,"augmented_forest":augmented_forest}
    
    def get_augmented_networks_intersections(self,prize_file,parameter_file=None):
        if parameter_file:
            self.load_parameter_files(parameter_file)
            
        Ws,Gs,Bs=set(),set(),set()
        parameters=list()
        for param in self.parameters:
            item=param.split("_")
            parameters.append({item[0].lower():item[1],item[2].lower():item[3],item[4].lower():item[5]})
        for param in parameters:
            Ws.add(float(param['w']))
            Gs.add(float(param['g']))
            Bs.add(float(param['b']))
            
        results = self.OI2_graph.grid_search(prize_file, Ws, Bs, Gs)
        self.augmented_pool_edges=[]
        
        for param in self.parameters:
            #print (param,len(set(results[param]["augmented_forest"].edges())),len(results[param]["augmented_forest"].edges()))
            temp_edges=results[param]["augmented_forest"].edges()
            temp_edges_set=set()
            for each in temp_edges:
                temp_edges_set.add((each[0],each[1]))
                temp_edges_set.add((each[1],each[0]))
            self.augmented_pool_edges.append(temp_edges_set)
        intersect=self.augmented_pool_edges.pop()
        for each_edge_set in self.augmented_pool_edges:
            intersect=intersect.intersection(each_edge_set)
        del self.augmented_pool_edges,temp_edges,results,temp_edges_set,parameters,parameter_file, Ws,Gs,Bs
        # n = gc.collect()
        network=nx.from_edgelist(intersect)
        network.remove_edges_from(nx.selfloop_edges(network))
        return network.edges()
    
    def write_augmented_networks_intersection(self,generated_network_name,prize_file,parameter_file=None):
        edges=self.get_augmented_networks_intersections(prize_file,parameter_file)
        with open (f'{generated_network_name}.res',"w") as f:
            for edge in edges:
                f.writelines(f'{edge[0]}\t{edge[1]}\n')
                

        