

import os, networkx as nx, pandas as pd, numpy as np, sys, multiprocessing as mp

import ROC

class RandomWalk():
    def __init__(self, network,edge_weight=None):
        self.network=network
        
        self.network_node_lngth=self.network.number_of_nodes()
        self.nodes=list(self.network.nodes())
        self.edge_weight=edge_weight
        
    def load_initial_nodes(self,initial_nodes):
        self.initial_nodes=[i for i in initial_nodes if i in self.nodes ]
        weight=[1/(len(self.initial_nodes))]*len(self.initial_nodes)
        self.initial_weight=dict(zip(self.initial_nodes,weight))
        return self.initial_weight
    
    def random_walk(self, alpha=0.50, threshould=0.20):        
        pagerank = nx.pagerank(self.network,alpha=alpha,personalization=self.initial_weight,weight=self.edge_weight)        
        node_thrshold=1/self.network_node_lngth
        returned_nodes=[i for i,j in pagerank.items() if j >= node_thrshold]
        Resolved_u=[]
        for u in pagerank:
            denom = self.network.degree(u)
            for v in self.network.neighbors(u): #  neighbors
                
                if (u not in returned_nodes) or (v not in returned_nodes):
                    self.network[u][v]['flux'] = 0
                    continue
                    
                if v in Resolved_u:
                    continue    
                    
                if self.edge_weight:                    
                    flux1 = pagerank[u]*self.network[u][v][self.edge_weight]/denom
                    flux2 = pagerank[v]*self.network[u][v][self.edge_weight]/denom
                    self.network[u][v]['flux']=min([flux1,flux2])

                else:
                    flux1=pagerank[u]*1/denom
                    flux2=pagerank[v]*1/denom                    
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
            Returned_edges[(u,v)] = self.network[u][v]['neglog_flux']
            flux_sum+=self.network[u][v]['flux']
            if flux_sum/tot_flux > threshould:
                print('theshold of %f limits predictions to %d edges'  %(threshould,len(Returned_edges)))
                break        
        return Returned_edges  
        
        
    
    # def __random_walk_rs_ss(self, pathway_edges, alpha=0.80, thres=0.80):
        
    #     pathway=nx.from_edgelist(pathway_edges)
    #     path_nodes=pathway.nodes()
        
    #     pagerank = nx.pagerank(self.network,alpha=alpha,personalization=self.initial_weight,weight=self.edge_weight)
        
    #     node_thrshold=1/self.network_node_lngth
    #     returned_nodes=[i for i,j in pagerank.items() if j >= node_thrshold]
        
    #     roc_node=ROC.ROC_node_based(self.network)
        
    #     roc_node.define_initial_nodes(self.initial_weight.keys())
    #     roc_node.define_pathways_nodes(path_nodes)
    #     roc_node.define_outcome_nodes(returned_nodes)
    #     roc_node.calculate_all_values()
    #     N_Recall=roc_node.calculate_recall()
    #     # N_FPR=roc_node.calculate_false_positive_rate()
    #     N_Precision=roc_node.calculate_precision()
    #     N_F1=roc_node.calculate_F1_score()
    #     # N_mcc=roc_node.calculate_MCC()
        
    #     print ("Node :","Recall", N_Recall,"N_Precision", N_Precision,"N_F1:",N_F1)


    #     Resolved_u=[]
    #     for u in pagerank:
    #         denom = self.network.degree(u)
    #         for v in self.network.neighbors(u): #  neighbors

    #             if v in Resolved_u:
    #                 continue      


    #             if (u not in returned_nodes) or (v not in returned_nodes):
    #                 self.network[u][v]['flux'] = 0
    #                 continue


    #             if self.edge_weight:
    #                 flux1 = pagerank[u]*self.network[u][v][self.edge_weight]/denom
    #                 flux2 = pagerank[v]*self.network[u][v][self.edge_weight]/denom
    #                 self.network[u][v]['flux']=min([flux1,flux2])


    #             else:
    #                 flux1=pagerank[u]*1/denom
    #                 flux2=pagerank[v]*1/denom
    #                 self.network[u][v]['flux'] = min([flux1,flux2])

    #         Resolved_u.append(u)

    #     tot_flux = 0
    #     for u,v in self.network.edges(): 

    #         if self.network[u][v]['flux'] == 0 :
    #             self.network[u][v]['neglog_flux'] = sys.float_info.max ## maximum float value
    #         else:
    #             self.network[u][v]['neglog_flux'] = - (np.log10(self.network[u][v]['flux'])) 
    #             tot_flux+=self.network[u][v]['flux']

    #     ROC_edge=ROC.ROC_edge_based(self.network)
    #     ROC_edge.define_pathway_edges(pathway_edges)
        
    #     Returned_edges = {} #dictionary of (u,v): -log(flux*flux)
    #     flux_sum = 0


    #     for u,v,d in sorted(self.network.edges(data=True), key=lambda t: t[2]['neglog_flux']):  

    #         Returned_edges[(u,v)] = self.network[u][v]['neglog_flux']
    #         flux_sum+=self.network[u][v]['flux']
    #         if flux_sum/tot_flux > thres:
    #             print('theshold of %f limits predictions to %d edges'  %(thres,len(Returned_edges)))
    #             break
    #     ROC_edge.define_outcome_edges(Returned_edges)

    #     ROC_edge.calculate_all_values()

    #     Recall=ROC_edge.calculate_recall()
    #     # FPR=ROC_edge.calculate_false_positive_rate()
    #     Precision=ROC_edge.calculate_precision()
    #     F1=ROC_edge.calculate_F1_score()
    #     # mcc=ROC_edge.calculate_MCC()
    #     print ("Recall:", Recall,"F1:",F1, "Precision:", Precision)
        
        
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
            
            pagerank = nx.pagerank(self.network,alpha=alpha,personalization=self.initial_weight,weight=self.edge_weight)

            node_thrshold=1/self.network_node_lngth
            returned_nodes=[i for i,j in pagerank.items() if j >= node_thrshold]
        
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
            for u in pagerank:
                denom = self.network.degree(u)
                for v in self.network.neighbors(u): #  neighbors

                    if v in Resolved_u:
                        continue      
                        
                    if (u not in returned_nodes) or (v not in returned_nodes):
                        self.network[u][v]['flux'] = 0
                        continue
                        
                    if self.edge_weight:
                        flux1 = pagerank[u]*self.network[u][v][self.edge_weight]/denom
                        flux2 = pagerank[v]*self.network[u][v][self.edge_weight]/denom
                        self.network[u][v]['flux']=min([flux1,flux2])
                        
                    else:
                        flux1=pagerank[u]*1/denom
                        flux2=pagerank[v]*1/denom
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
                

                #print ("thres:",thres, "F1:",F1,"mcc:",mcc,"Recall:", Recall, "Precision:",Precision, 
                 #      "FPR:",FPR,"edge_len", len(Returned_edges))
        return output_df
