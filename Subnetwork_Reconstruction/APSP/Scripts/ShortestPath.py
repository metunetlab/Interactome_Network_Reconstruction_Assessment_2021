# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 06:17:51 2020

@author: mkaanarici
"""

import networkx as nx, numpy as np, seaborn as sns, matplotlib.pyplot as plt, os,pandas as pd,time,scipy.stats as st, pickle


# =============================================================================
# shortest pathway
# =============================================================================


class Shortest_Paths():
    def __init__(self, network):
        self.network=network
        self.shortest_path_dict={}
        self.nodes=list(self.network.nodes())
    
    def find_shortest_path_between(self,source,target):
        if source not in self.nodes or target not in self.nodes:
            return []
        if f'{source}_{target}' in self.shortest_path_dict.keys():
            #print ("we have")
            return self.shortest_path_dict[f'{source}_{target}']
        if f'{target}_{source}' in self.shortest_path_dict.keys():
            #print ("we have")
            return self.shortest_path_dict[f'{target}_{source}']
    

        try:
            shrt_path=list(nx.algorithms.shortest_paths.generic.all_shortest_paths(G=self.network,source=source,target=target))
        except nx.NetworkXNoPath:
            print (f'{source} cannot reach to {target}')
            self.shortest_path_dict[f'{source}_{target}']=[]
            return []
        self.shortest_path_dict[f'{source}_{target}']=shrt_path
        return shrt_path
    
    def find_shortest_path_among(self, sources, targets):
        paths=[]
        done=[]
        for source in sources:
            for target in targets:
                if source==target or f'{source}-{target}' in done or f'{target}-{source}' in done:
                    continue
                shrt_path=self.find_shortest_path_between(source=source,target=target)
                done.append(f'{source}-{target}')
                for i in shrt_path:
                    paths.append(i)
        return paths
    
    def write_shortest_path_edges_among(self, sources, targets, file_name):
        with open (file_name,"w") as f:
            edges=[]
            for source in sources:
                for target in targets:
                    if source==target:
                        continue
                    shrt_path=self.find_shortest_path_between(source=source,target=target)
                    
                    for pth in shrt_path:  
                        i=0
                        #print ("*",pth,len(pth))
                        while i+1<len(pth):
                            node1,node2=pth[i],pth[i+1]
                            # print (i,node1,node2)
                            i+=1                        
                            if [node1,node2] not in edges and [node2,node1] not in edges:
                                f.writelines(f'{node1}\t{node2}\n')
                                edges.append([node1,node2])
        return edges
        
    
    


# =============================================================================
# 
# =============================================================================





















