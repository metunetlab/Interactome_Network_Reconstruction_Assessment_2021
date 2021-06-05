

class ROC_node_based():
    def __init__(self, network):
        self.network=network
        self.nodes=set(self.network.nodes())
        
    def define_initial_nodes(self,initial_nodes):
        self.initial_nodes=set(initial_nodes).intersection(self.nodes)
        #print ("initial nodes   :", len(initial_nodes))
        #print ("intial nodes found in interactome:",len(self.initial_nodes) )
        self.space=self.nodes.difference(self.initial_nodes)
        #print (" considered nodes in interactome :", len(self.space))
        #print (" considered nodes in interactome  are composed of nodes apart from initial nodes")
        
        
    def define_pathways_nodes(self,nodes):
        #print ("pathway nodes :", len(nodes))
        self.pathway_nodes=self.space.intersection(set(nodes).difference(self.initial_nodes))
        #print ("considered pathway nodes :", len(self.pathway_nodes))
        self.pathway_negative_nodes=self.space.difference(self.pathway_nodes)
        #print ("nodes out of pathway :",len(self.pathway_negative_nodes))
        return self.pathway_nodes
    
    def define_outcome_nodes(self, nodes):       
        self.positive_nodes=set(nodes).difference(self.initial_nodes)
        #print ("nodes found by propagation :", len(self.positive_nodes))
        self.negative_nodes=self.space.difference(self.positive_nodes)
        #print ("nodes out of propagation :", len(self.negative_nodes))
        
    def find_true_positives(self):
        self.true_positives=self.positive_nodes.intersection(self.pathway_nodes)
        self.TP=len(self.true_positives)
        return self.true_positives
    
    def find_false_positives(self):
        self.false_positives=self.positive_nodes.difference(self.pathway_nodes)
        self.FP=len(self.false_positives)
        return self.false_positives
    
    def find_true_negatives(self):
        self.true_negatives=self.negative_nodes.intersection(self.pathway_negative_nodes)
        self.TN=len(self.true_negatives)
        return self.true_negatives
    
    def find_false_negatives(self):
        self.false_negatives=self.negative_nodes.intersection(self.pathway_nodes)
        self.FN=len(self.false_negatives)
        return self.false_negatives
    
    def calculate_all_values(self):
        self.find_false_negatives()
        self.find_false_positives()
        self.find_true_negatives()
        self.find_true_positives()
        return True
    
    def calculate_MCC(self):
        divedend= self.TP* self.TN- self.FP* self.FN
        enumarater=((self.TP+self.FP)*(self.TP+self.FN)*(self.TN+self.FP)*(self.TN+self.FN))**0.5
        if enumarater==0:
            return 0
        return divedend/enumarater
    


    
    def calculate_recall(self): # Recall: TP/(TP+FN)
        if (self.TP+self.FN)==0 or self.TP==0:
            return 0
        return self.TP/(self.TP+self.FN)
    
    def calculate_precision(self): # Precision TP/(TP+FP)
        if self.TP+self.FP==0:
            return 0
        else:
            return self.TP/(self.TP+self.FP)
    
    def calculate_F1_score(self): # F1-score： 2/(1/P+1/R)
        Precision=self.calculate_precision()
        Recall=self.calculate_recall()
        if Recall + Precision==0:
            return 0
        return 2*(Recall * Precision) / (Recall + Precision)
    
    def calculate_false_positive_rate(self): # FPR = FP / ( FP + TN )
        if self.FP==0 or (self.FP+self.TN)==0:
            return 0
        return self.FP/(self.FP+self.TN)
    
    
    

# =============================================================================
# edge based ROC
# =============================================================================

class ROC_edge_based():
    def __init__(self,network):
        self.network=network
        self.edges=set(self.network.edges())
        
    def define_pathway_edges(self,edges):
        pathway_PE1={(i[0],i[1]) for i in edges}.intersection(self.edges)        
        pathway_PE2={(i[1],i[0]) for i in edges}.intersection(self.edges)
        self.pathway_positive_edges=pathway_PE1.union(pathway_PE2)
        self.pathway_negative_edges=self.edges.difference(self.pathway_positive_edges)
        
        
    def define_outcome_edges(self,edges):
        PE1={(i[0],i[1]) for i in edges}.union({(i[1],i[0]) for i in edges})
        self.positive_edges=PE1.intersection(self.edges)
        self.negative_edges=self.edges.difference(self.positive_edges)
        
        
    def find_true_positives(self):
        self.true_positives=self.positive_edges.intersection(self.pathway_positive_edges)
        self.TP=len(self.true_positives)
        return self.true_positives
    
    def find_false_positives(self):
        self.false_positives=self.positive_edges.difference(self.pathway_positive_edges)
        self.FP=len(self.false_positives)
        return self.false_positives
    
    def find_true_negatives(self):
        self.true_negatives=self.negative_edges.intersection(self.pathway_negative_edges)
        self.TN=len(self.true_negatives)
        return self.true_negatives
    
    def find_false_negatives(self):
        self.false_negatives=self.negative_edges.intersection(self.pathway_positive_edges)
        self.FN=len(self.false_negatives)
        return self.false_negatives
    
    def calculate_all_values(self):
        self.find_false_negatives()
        self.find_false_positives()
        self.find_true_negatives()
        self.find_true_positives()
        return True
    
    
    
    
    def calculate_recall(self): # Recall: TP/(TP+FN)
        if self.TP==0 or (self.TP+self.FN)==0:
            return 0
        return self.TP/(self.TP+self.FN)
    
    def calculate_precision(self): # Precision TP/(TP+FP)
        if self.TP+self.FP==0:
            return 0
        else:
            return self.TP/(self.TP+self.FP)
    
    def calculate_F1_score(self): # F1-score： 2/(1/P+1/R)
        Precision=self.calculate_precision()
        Recall=self.calculate_recall()
        if Recall + Precision==0:
            return 0
        return 2*(Recall * Precision) / (Recall + Precision)
    
    def calculate_false_positive_rate(self): # FPR = FP / ( FP + TN )
        if self.FP==0 or (self.FP+self.TN)==0:
            return 0
        return self.FP/(self.FP+self.TN)        
        
    def calculate_MCC(self):
        divedend= self.TP* self.TN- self.FP* self.FN
        enumarater=((self.TP+self.FP)*(self.TP+self.FN)*(self.TN+self.FP)*(self.TN+self.FN))**0.5
        if enumarater==0:
            return 0
        return divedend/enumarater


