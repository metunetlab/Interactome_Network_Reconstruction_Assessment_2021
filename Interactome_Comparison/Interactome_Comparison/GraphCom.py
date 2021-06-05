# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:13:58 2020

@author: mkaanarici
"""

import pandas as pd, numpy as np, networkx as nx,seaborn as sns,matplotlib.pyplot as plt, collections,os
from scipy import stats



class GraphCom():
    """
    GraphCom object is composed of graphs, which are defined with networkx.Graph. 
    
    """
    def __init__(self, Interactomes=None):
        """
        Graphs can be defined initially in the form zipped tuple in list.  
        Interactomes: the list in which the name of graph and graph in tuples or None. 
        If the Interactome list is not be provided by setting  as None, graphes are defined with 'GraphCom.add_graph' function.
        """
        self.Graphs=[] # each graph is listed here
        self.index=0 # unique index number is given for each graph
        self.names=[] # graph names are listed here
        self.indexer_reverse,self.indexer={},{} #index and name is converted into each other by using dictionaries
        self.attributes={} # edge atributes are kept in dictionary with graph name 
        
        self.node_sim_matrix_check=False
        self.edge_sim_matrix_check=False
        self.correlation_matrix_check=False
        self.Intersection_matrix_check=False
        
        if Interactomes: # interactome can be given in zipped name and grah in list or can be empty 
            for name,graph in Interactomes:
                self.add_graph(name,graph)

            
    def get_indexer(self):
        """ @the dictionary is the output in which index is key, graph name is value """
        return self.indexer 
    
    def get_index(self,interactome): # the name of graph is given, index is output
        """ interactome (str):the name of interactome return the index """
        return self.indexer_reverse[interactome]
    
    def get_network(self,interactome): # the name of graph is given, networkx is output
        """ interactome (str):the name of interactome return networkx graph """
        return self.Graphs[self.get_index(interactome)]
    
    def get_name_list(self):
        """ return the name list of interactome"""
        return self.names
    
    def get_name(self,index): # the index of graph is given and the name of graph is returned
        """index(int) : the index of list, return the name of interactome in the given index"""
        return self.indexer[index]
    
    def add_graph_from_CSV(self,file,name,sep=",",source="UniprotID_1", target="UniprotID_2",edge_attr=None):
        """"""
        DataFrame=pd.read_csv(file, sep=sep).dropna()
        NX=nx.from_pandas_edgelist(DataFrame, source, target,edge_attr=edge_attr,create_using=nx.Graph)
        self.add_graph(name,NX)
        

    def add_graph(self,name,graph):         
        """
        name: the name of graph is given in str
        graph: the graph is given in netwrokx.Graph 
        """
        if isinstance(graph,nx.classes.graph.Graph):
                self.Graphs.append(graph)                
                self.names.append(name)
                self.indexer[self.index]=name
                self.indexer_reverse[name]=self.index
                self.index+=1
        else:
            raise ValueError("items in Interactomes must be networkx.classes.graph.Graph")
            
    def remove_graph(self,interactome):
        """
        to delete graph
        """
        index=self.get_index(interactome)
        del(self.Graphs[index])
        self.names.remove(interactome)        
        self.index-=1
        self.indexer=dict(zip(list(range(0,len(self.names))),self.names))
        self.indexer_reverse=dict(zip(self.names,list(range(0,len(self.names)))))
        self.node_sim_matrix_check=False
        self.edge_sim_matrix_check=False
        self.correlation_matrix_check=False
        self.Intersection_matrix_check=False
    
    def get_nodes_edges_count(self,graph=None):
        """
        if the certain network name is not given,  the number of nodes, the number of edges are returned in the dictionary
        where graph names construct key. 
        Providing the certain netwrok name; the number of nodes and the number of edges are respectively returned in tuple
        """
        if graph is None:
            Node_edge={}
            for name in self.names:
                Node_edge[name]=(self.get_network(name).number_of_nodes(),self.get_network(name).number_of_edges())
            return np.array(Node_edge)
        else:
            return (self.get_network(graph).number_of_nodes(),self.get_network(graph).number_of_edges())
    
    def get_edges_between(self,interactome=None,min_score=None,max_score=None):
        """
        find out edges in interactome having score between min_score and max_score
        interactome (str) : the name of interactome
        min_score (int): minimum score for edges
        max_score (int): maximum score fot edges
        """
        if interactome:
            NX=self.get_network(interactome)
            PD_DF=nx.convert_matrix.to_pandas_edgelist(NX)
            if min_score:
                PD_DF=pandas.loc[pandas['mi-score']>min_score]
            if max_score:
                PD_DF=pandas.loc[pandas['mi-score']<max_score]
            return PD_DF
        
    def visualize_edge_dispersion_with_violin_plot(self,name=None):
        DataFrame=pd.DataFrame(columns=['Databases','Scores'])
        
        for nm in self.names:    
            temp_df=pd.DataFrame(columns=['Databases','Scores'])
            graph=self.get_network(nm)
            score1=None
            graph_df=nx.convert_matrix.to_pandas_edgelist(graph)
            #print(nm,graph.number_of_edges())
            Columns1=graph_df.columns
            if len(Columns1)>2:
                score1=Columns1[2]
            else:
                #print (f'{nm} doesn not have any attribute score')
                continue
            #print(name,score1)
            temp_df['Scores']=graph_df[score1].astype(float)
            lngth=len(temp_df.index)
            names=[nm]*lngth
            #print (lngth,len(nm))
            temp_df['Databases']=names
            DataFrame=DataFrame.append(temp_df)
        # print(DataFrame)
        fig=plt.figure(figsize=(30,20),dpi=300)
        ax=sns.violinplot(x='Databases',y='Scores', data=DataFrame,split=True,scale='width',  inner="quartiles")
        #ax.set_xticklabels(ax.get_ticklabels(),fontsize=40,fontweight='bold',rotation=90)
        #ax.set_yticklabels(ax.get_yticklabels(),fontsize=20,fontweight='bold')
        plt.ylabel("Score",fontsize=30)
        plt.xticks(fontsize=40,fontweight='bold',rotation=90)
        plt.yticks(fontsize=30)
        plt.xlabel('')
        if name is None:
            plt.savefig("edge_dispersion_with_violin_plot.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
            

            
    
    def visualize_categorized_edges_count(self,x,y,z,name=None):
        """
        visulalize the categorized edge counts in mathplotlib bar plot
        x,y,z are the boundaries for bins
        """
        x,y,z=sorted([x,y,z])
        Bin1,Bin2,Bin3=[],[],[]
        names = []
        barWidth = 0.5
        for grph in self.names:            
            graph=self.get_network(grph)
            score1=None
            graph_df=nx.convert_matrix.to_pandas_edgelist(graph)
            #print(graph,graph.number_of_edges())
            Columns1=graph_df.columns
            if len(Columns1)>2:
                score1=Columns1[2]
            else:
                print (f'{grph} doesn not have any attribute score')
                continue
            #print(name,score1)
            graph_df[score1]=graph_df[score1].astype(float)
            graph_df1=graph_df.loc[(graph_df[score1]<y) & (graph_df[score1]>=x)]
            Bin1.append(len(graph_df1.index))
            graph_df2=graph_df.loc[(graph_df[score1]<z) & (graph_df[score1]>=y)]
            Bin2.append(len(graph_df2.index))
            graph_df3=graph_df.loc[(graph_df[score1]>=z)]
            Bin3.append(len(graph_df3.index))
            names.append(grph)
        r=list(range(0,len(names)))
        bars = np.add(Bin1,Bin2).tolist()
        plt.figure(figsize=(30,15), dpi= 300)
        plt.bar(r, Bin1, color='mistyrose', edgecolor='black', width=barWidth,label=f'{x} - {y}')
        plt.bar(r, Bin2, bottom=Bin1, color='salmon', edgecolor='black', width=barWidth,label=f'{y} - {z}')
        plt.bar(r, Bin3, bottom=bars, color='red', edgecolor='black', width=barWidth,label=f' > {z}')
        plt.xticks(r, names,fontsize=40, fontweight='bold',rotation=90)
        plt.yticks(fontsize=40, fontweight='bold')
        plt.ylabel("The number of edges",fontsize=40,fontweight='bold')
        plt.legend(fontsize=20)
        if name is None:
            plt.savefig("categorized_edges_count.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()

       
    
    def visualize_nodes_edges_counts(self,name=None):
        """
        visulalize both edges and nodes counts in bar plot
        x,y,z are the boundaries for bins
        """
        # Dictionary=self.get_nodes_edges_count()
        # DataFrame=pd.DataFrame()
        Node_count=[]
        Edge_count=[]
        for graph in self.names:
            temp=self.get_nodes_edges_count(graph)
            Node_count.append(temp[0])
            Edge_count.append(temp[1])
#        print (len(self.names),len(Node_count),len(Edge_count))
        DF=pd.DataFrame()
        DF["Databases"]=self.names
        DF["Nodes Count"]=Node_count
        DF["Edge Count"]=Edge_count
        
        plt.figure(figsize=(30,15), dpi= 300)
        DF.plot(kind='bar')
        if name is None:
            plt.savefig("nodes_edges_counts.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
    
    def visualize_nodes_count(self,name=None):
        
        """
        visulalize the node counts in mathplotlib bar plot
        x,y,z are the boundaries for bins
        """
        Count=[]
        for graph in self.names:
            Count.append(self.get_nodes_edges_count(graph)[0])
        Data={"Databases":self.names,"The number of Nodes":Count}               
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.barplot(x="Databases",y="The number of Nodes",data=Data)
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.ylim(0,21000)
        plt.xlabel('')
        plt.ylabel('The number of nodes',fontsize=40)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.0f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')


        if name is None:
            plt.savefig("nodes_count.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
    
    def visualize_edges_count(self,name=None):
        
        """
        visulalize the edge counts in bar plot
        x,y,z are the boundaries for bins
        """
        Counts=[]
        for graph in self.names:
            Counts.append(self.get_nodes_edges_count(graph)[1])
        Data={"Databases":self.names,"The number of Edges":Counts}                      
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.barplot(x="Databases",y="The number of Edges",data=Data)
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel('The number of Edges',fontsize=40)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.0f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')


        if name is None:
            plt.savefig("edges_count.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        #barplot=sns.barplot(x="Databases",y="The number of Edges",data=Data)
        #barplot.set_xticklabels(barplot.get_xticklabels(),rotation=90)
        #barplot.set( ylabel="The number of Edges")
        #return barplot
    
    
                
                    ########################################
                    #                                      #
                    #           Intersections              #
                    #                                      #
                    ########################################
    
    def get_intersection_between(self,interactome1,interactome2,format="networkx"):
        """
        The intersection of two given networks is generated. edge attributes coming from 2 networks are kept        
        output can be in networkx.Graph format or pandas.DataFrame
        
        interactome1 (str) and interactome2 (str) are the name of interactome find out their intersecting edges
        
        format (str): output format can be "networkx" or "Pandas.DataFrame" 
        
        """
        
        graph1=self.get_network(interactome1)
        graph2=self.get_network(interactome2)

        if graph1==graph2:
            return graph1

        DataFrame1=nx.convert_matrix.to_pandas_edgelist(graph1)
        DataFrame2=nx.convert_matrix.to_pandas_edgelist(graph2)
        #DataFrame_Intersection=pd.merge(DataFrame1,DataFrame2,how='inner')

        intersection1=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['source','target'])
        intersection2=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['target','source'])
        intersection2=intersection2.drop(columns=['source_y', 'target_y'])
        intersection2=intersection2.rename(columns={'source_x':"source",'target_x':"target"})     
        DataFrame_Intersection=intersection1.append(intersection2)
        DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()        

        intersection=nx.from_pandas_edgelist(DataFrame_Intersection,edge_attr=True,create_using=nx.Graph)
        if format=="networkx":
            return intersection
        elif format=="Pandas.DataFrame":
            return(nx.convert_matrix.to_pandas_edgelist(intersection))


    
    def get_intersections(self):
        
        """ the all interscetions of networkks are generated and place in numpy matrix and in networkx.graph format"""
        
        self.Intersection_matrix=np.identity(len(self.Graphs),dtype=object)
  #      print (self.Intersection_matrix)
        i1=-1
        for graph1 in self.Graphs:  
            i1+=1
            i2=-1
            for graph2 in self.Graphs:
                i2+=1         
                if graph1==graph2:
                    self.Intersection_matrix[i1,i2]=graph1
                    continue
                DataFrame1=nx.convert_matrix.to_pandas_edgelist(graph1)
                DataFrame2=nx.convert_matrix.to_pandas_edgelist(graph2)
                #DataFrame_Intersection=pd.merge(DataFrame1,DataFrame2,how='inner')
                
                
                intersection1=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['source','target'])
                intersection2=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['target','source'])
                intersection2=intersection2.drop(columns=['source_y', 'target_y'])
                intersection2=intersection2.rename(columns={'source_x':"source",'target_x':"target"})     
                DataFrame_Intersection=intersection1.append(intersection2)
                DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()   
                
                intersection=nx.from_pandas_edgelist(DataFrame_Intersection,create_using=nx.Graph)
                self.Intersection_matrix[i1,i2]=intersection                
                self.Intersection_matrix[i2,i1]=intersection
                if i2>i1:
                    break
        self.Intersection_matrix_check=True        
        return self.Intersection_matrix
        
                    
                    ########################################
                    #                                      #
                    #            NODE & EDGE               #
                    #            SIMILARTIES               #
                    #                                      #
                    ########################################
        
    def get_node_similarity_matrix(self):
        
        """
        The nodes similarity scores are kept in numpy matrix. Similarity score for each intersection 
        is calculated by given formula:
        
             Similarity_Score=the_number_of_node_in_intersection/min(Node_count_in_Network1,Node_count_in_Network2)
        """
        if self.Intersection_matrix_check==False:
            self.get_intersections()
        self.node_sim_matrix=np.identity(len(self.Graphs))
            
        i1=-1
        for graph1 in self.Graphs:  
            i1+=1            
            i2=-1
            for graph2 in self.Graphs:
                i2+=1       
                if graph1==graph2:
                    continue
                node1_count=graph1.number_of_nodes()
                node2_count=graph2.number_of_nodes()
                node_int_count=self.Intersection_matrix[i1,i2].number_of_nodes()
                Similarity_scores=node_int_count/min(node1_count,node2_count)
                self.node_sim_matrix[i1,i2]=Similarity_scores
                self.node_sim_matrix[i2,i1]=Similarity_scores
                if i2>i1:
                    break
        self.node_sim_matrix_check=True
        return self.node_sim_matrix
    
    def get_node_similarity_between(self, interactome1,interactome2):
        """The node similarity score between given interactomes' names is calculated by given formula:
            Similarity_Score=the_number_of_node_in_intersection/min(Node_count_in_Network1,Node_count_in_Network2)"""
        
        i1=self.get_index(interactome1)
        i2=self.get_index(interactome2)
        if self.node_sim_matrix_check:
            return self.node_sim_matrix[i1][i2]
        else:
            self.get_node_similarity_matrix()
            return self.node_sim_matrix[i1][i2]
            
    
    def get_node_similarity_matrix_to_pandas_DataFrame(self):
        """
        The nodes similarity scores are kept in pandas.DataFrame. Similarity score for each intersection 
        is calculated by given formula:
            Similarity_Score=the_number_of_node_in_intersection/min(Node_count_in_Network1,Node_count_in_Network2)
        """ 
        if self.node_sim_matrix_check:
            return pd.DataFrame(self.node_sim_matrix,columns=self.names,index=self.names)
        else:
            self.get_node_similarity_matrix()
            return pd.DataFrame(self.node_sim_matrix,columns=self.names,index=self.names)
     
    def visualize_node_similarity_matrix(self,cmap="Blues", name=None):
        """
        The nodes similarity scores are visualized in seaborn.heatmap. Similarity score for each intersection 
        is calculated by given formula:
            Similarity_Score=the_number_of_node_in_intersection/min(Node_count_in_Network1,Node_count_in_Network2)
        """ 
        if self.node_sim_matrix_check:
            if name is None:
                return self._visualize_matrix(self.node_sim_matrix,name="node_similarity_matrix",cmap=cmap)
            else:
                return self._visualize_matrix(self.node_sim_matrix,name=name,cmap=cmap)
                
        else:
            self.get_node_similarity_matrix()
            if name is None:
                return self._visualize_matrix(self.node_sim_matrix,"node_similarity_matrix",cmap=cmap)
            else:
                return self._visualize_matrix(self.node_sim_matrix,name=name,cmap=cmap)
    
    
    def get_edge_similarity_matrix(self):
        """
        The edges similarity scores are kept in numpy matrix. Similarity score for each intersection 
        is calculated by given formula:
        
             Similarity_Score=the_number_of_edge_in_intersection/min(Edge_count_in_Network1,Edge_count_in_Network2)
        """
        self.edge_sim_matrix=np.identity(len(self.Graphs))
        i1=-1
        if self.Intersection_matrix_check==False:
            self.get_intersections()
        for graph1 in self.Graphs:  
            i1+=1            
            i2=-1
            for graph2 in self.Graphs:
                i2+=1       
                if graph1==graph2:
                    continue
                edge1_count=graph1.number_of_edges()
                edge2_count=graph2.number_of_edges()
                edge_int_count=self.Intersection_matrix[i1,i2].number_of_edges()
                Similarity_scores=edge_int_count/min(edge1_count,edge2_count)
                self.edge_sim_matrix[i1,i2]=Similarity_scores
                self.edge_sim_matrix[i2,i1]=Similarity_scores
                if i2>i1:
                    break
        self.edge_sim_matrix_check=True        
        return self.edge_sim_matrix

    def get_edge_similarity_between(self, interactome1,interactome2):
        """The edge similarity score between given interactomes' names is calculated by given formula:
            Similarity_Score=the_number_of_edge_in_intersection/min(edge_count_in_Network1,edge_count_in_Network2)"""
        i1=self.get_index(interactome1)
        i2=self.get_index(interactome2)
        if self.edge_sim_matrix_check:
            return self.edge_sim_matrix[i1][i2]
        else:
            self.get_edge_similarity_matrix()
            return self.edge_sim_matrix[i1][i2]
    
    
    def get_edge_similarity_matrix_to_pandas_DataFrame(self):
        """
        The edges similarity scores are kept in pandas.DataFrame. Similarity score for each intersection 
        is calculated by given formula:
        
             Similarity_Score=the_number_of_edge_in_intersection/min(Edge_count_in_Network1,Edge_count_in_Network2)
        """
        if self.edge_sim_matrix_check:
            return pd.DataFrame(self.edge_sim_matrix,columns=self.names,index=self.names)
        else:
            self.get_edge_similarity_matrix()
            return pd.DataFrame(self.edge_sim_matrix,columns=self.names,index=self.names)
            
    
    def visualize_edge_similarity_matrix(self,name=None,cmap="Blues"):
        """
        The edges similarity scores are visualized in seaborn.heatmap.. Similarity score for each intersection 
        is calculated by given formula:
        
             Similarity_Score=the_number_of_edge_in_intersection/min(Edge_count_in_Network1,Edge_count_in_Network2)
        """
        if self.edge_sim_matrix_check:
            if name is None:
                return self._visualize_matrix(self.edge_sim_matrix,"edge_similarity_matrix",cmap)
            else:
                return self._visualize_matrix(self.edge_sim_matrix,name,cmap)
        else:
            self.get_edge_similarity_matrix()
            if name is None:
                return self._visualize_matrix(self.edge_sim_matrix,"edge_similarity_matrix",cmap)
            else:
                return self._visualize_matrix(self.edge_sim_matrix,name,cmap)
        
        
    def visualize_node_edge_similarity(self,cmap='bwr',name=None):        
        DF=pd.DataFrame()
        X_name=[]
        Y_name=[]
        Edge_sim=[]
        Node_sim=[]        
        for x in self.get_name_list():
            for y in self.get_name_list():
                if x==y:
                    node_score=0
                    edge_score=0
                else:                   
                    node_score=self.get_node_similarity_between(x,y)
                    edge_score=self.get_edge_similarity_between(x,y)                
                Edge_sim.append(edge_score)
                Node_sim.append(node_score)
                X_name.append(x)
                Y_name.append(y) 
        DF['X']=X_name
        DF['Y']=Y_name
        DF["Edge_sim"]=Edge_sim
        DF["Node_sim"]=Node_sim   
        
        print ("max edge",DF['Edge_sim'].max(),"max_node",DF['Node_sim'].max())
        
        fig, ax = plt.subplots(figsize=(5,5),dpi=300)
        ax=sns.scatterplot(y='Y',x='X', data=DF,size='Edge_sim', hue='Node_sim',sizes=(0,400),palette=(cmap))

        #plt.legend(fontsize=8,title_fontsize=10)
        plt.xticks(rotation=90,fontsize=20,fontweight='bold')
        plt.yticks(fontsize=20,fontweight='bold')
        ax.set_xlabel('', fontsize=40,fontweight='bold')
        ax.set_ylabel("", fontsize=40,fontweight='bold')
        plt.legend(bbox_to_anchor=(1,0, 0.4,0.9),fontsize=16,handlelength=0.8,labelspacing=0.2)
        if name is None:
            plt.savefig("node_edge_similarity.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
        
        
        
        # fig, ax = plt.subplots(figsize=(15,15),dpi=50)
        # scatter=ax.scatter(DF['X'],DF['Y'],c=DF['Edge_sim'],s=DF['Node_sim'], sizes=(0,10,100,1000),alpha=1,cmap=cmap)

        # legend1 = plt.legend(*scatter.legend_elements(),
        #                     bbox_to_anchor=(0.8,-0.5, 0.5,1),fontsize='x-large')
        # legend1.set_title("Edge Similarities",prop = {'size':'xx-large'})
        # ax.add_artist(legend1)

        # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.8)
        # legend2 = plt.legend(handles, labels,bbox_to_anchor=(0.8,0, 0.5,1),fontsize='x-large' )

        # legend2.set_title('Node Similarities',prop = {'size':'xx-large'})
        # plt.xticks(rotation=90,fontsize='xx-large',fontweight='bold')
        # plt.yticks(fontsize='xx-large',fontweight='bold')
        
            
                    ########################################
                    #                                      #
                    #          Score Correlation           #
                    #                                      #
                    ########################################
                    
                    
    def get_correlation_between (self, interactome1,interactome2,output=None):
        """
        The correlation between two interactomes concerning edge attributes are calculated
        """
        
        graph1=self.get_network(interactome1)
        graph1_df=nx.convert_matrix.to_pandas_edgelist(graph1)
        
        Columns1=graph1_df.columns
        if len(Columns1)>2:
            score1=Columns1[-1]
        else:
            # print (f'{interactome1} doesn not have any attribute score')
            return 0      
        
        graph2=self.get_network(interactome2)
        graph2_df=nx.convert_matrix.to_pandas_edgelist(graph2)
        Columns2=graph2_df.columns
        if len(Columns2)>2:
            score2=Columns2[-1]
        else:
            # print (f'{interactome2} doesn not have any attribute score')
            return 0  
            
        if interactome1==interactome2:
            if output=='coefficient':
                return 1
            
            elif output=="p-value":
                return 0
            else:
                return (1,0)
        
        DataFrame1=nx.convert_matrix.to_pandas_edgelist(graph1)
        DataFrame2=nx.convert_matrix.to_pandas_edgelist(graph2)
        
        intersection1=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['source','target'])
        intersection2=DataFrame1.merge(DataFrame2, how='inner',left_on=['source','target'],right_on=['target','source'])
        intersection2=intersection2.drop(columns=['source_y', 'target_y'])
        intersection2=intersection2.rename(columns={'source_x':"source",'target_x':"target"})     
        DataFrame_Intersection=intersection1.append(intersection2)
        #DataFrame_Intersection=DataFrame_Intersection.drop_dublicates()
        score1,score2=DataFrame_Intersection.columns[2:]
        
        if len(DataFrame_Intersection.index)==0:
            return (0,1)
        
        if output=='coefficient':
            # print (DataFrame_Intersection.columns)
            # print(stats.pearsonr(DataFrame_Intersection[score1],DataFrame_Intersection[score2]))
            return stats.pearsonr(DataFrame_Intersection[score1],DataFrame_Intersection[score2])[0]
        elif output=="p-value":
            return stats.pearsonr(DataFrame_Intersection[score1],DataFrame_Intersection[score2])[1]
        else:
            return stats.pearsonr(DataFrame_Intersection[score1],DataFrame_Intersection[score2])
    
    def get_correlation_matrix(self,matrix='coefficient'):
        """
        The correlation between all interactomes concerning edge attributes are calculated and kept in numpy matrix
        """
        self.correlation_matrix=np.identity(len(self.Graphs))
        if self.Intersection_matrix_check==False:
            self.get_intersections()
        i1=-1
        for graph1 in self.names:  
            i1+=1            
            i2=-1
            for graph2 in self.names:
                # print (graph1,graph2)
                i2+=1
                if matrix=='coefficient':
                    self.correlation_matrix[i1,i2]=self.get_correlation_between(graph1,graph2,output='coefficient')
                    self.correlation_matrix[i2,i1]=self.correlation_matrix[i1,i2]
                    self.correlation_matrix_check=True
                elif matrix=="p-value":
                    self.correlation_matrix[i1,i2]=self.get_correlation_between(graph1,graph2,output="p-value")
                    self.correlation_matrix[i2,i1]=self.correlation_matrix[i1,i2]
                if i2>i1:
                    break
        return self.correlation_matrix
    
        
    def get_correlation_matrix_to_pandas_DataFrame(self,matrix='coefficient'):
        """
        The correlation between all interactomes concerning edge attributes are calculated and kept in pandas.DataFrame
        """
        self.get_correlation_matrix(matrix=matrix)
        return pd.DataFrame(self.correlation_matrix,columns=self.names,index=self.names)
            
    
    def visualize_confidence_scores_correlation_matrix(self, name=None,cmap="Blues"):
        """
        The correlation between all interactomes concerning edge attributes are calculated and visualized in seanborn.heatmap.
        """
        if self.correlation_matrix_check:
            if name is None:
                return self._visualize_matrix(self.correlation_matrix, name="confidence_scores_correlation_matrix",cmap=cmap)
            else:
                return self._visualize_matrix(self.correlation_matrix, name=name,cmap=cmap)
            
        else:
            self.get_correlation_matrix()
            if name is None:
                name="confidence_scores_correlation_matrix"
                return self._visualize_matrix(self.correlation_matrix,name=name,cmap=cmap)
            else:
                return self._visualize_matrix(self.correlation_matrix,name=name,cmap=cmap)
            
    
    def _visualize_matrix(self,matrix,name,cmap="Blues"):     
        
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.heatmap(matrix,xticklabels=self.names,yticklabels=self.names,cmap=cmap)        
        plt.xticks(rotation=90,fontsize=20,fontweight='bold')
        #ax.set_yticklabels(ax.get_yticklabels(), rotation=45, fontsize =20,fontweight='bold')
        plt.yticks(rotation=45,fontsize=20,fontweight='bold')
        plt.savefig("%s.svg"%name, format='svg')
    
                
                    ########################################
                    #                                      #
                    #                DEGREE                #
                    #            DISTRUBUTION              #   
                    #                                      #
                    ########################################
    
    
    def visualize_degree_distrubution(self, graph,name=None):
        Graph=self.get_network(graph)
        Degree=Graph.degree()
        Degree_sequence = sorted([d for n, d in Degree], reverse=True)
        sns.distplot(Degree_sequence,hist=False,label=graph)
        if name is None:
            plt.savefig("%s_degree_distrubution.svg"%graph, format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
    def visualize_all_degree_distribution(self,xlim=None,name=None):
        plt.figure(figsize=(30,15), dpi= 300)
        if xlim:
            plt.xlim(0,xlim)
        for graph in self.names:
            Graph=self.get_network(graph)
            Degree=Graph.degree()
            Degree_sequence = sorted([d for n, d in Degree], reverse=True)
            sns.distplot(Degree_sequence,hist=False,label=graph)
        plt.legend(fontsize=20)
        plt.xlabel("Degree",fontsize=40)
        plt.xticks(fontsize=20)
        plt.ylabel("Relative Frequency",fontsize=40)
        plt.yticks(fontsize=20)
        if name is None:
            plt.savefig("all_degree_distribution.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
            
    def visualize_all_log_degree_distribution(self,xlim=None,ylim=None,name=None):
        plt.figure(figsize=(30,15), dpi= 300)
        for graph in self.names:
            Graph=self.get_network(graph)
            Degree=Graph.degree()
            Degree_sequence = sorted([d for n, d in Degree], reverse=True)
            DegreeCount = np.array(list(collections.Counter(Degree_sequence).items()))
            DegreeCount=DegreeCount.astype(float)
            total=sum(DegreeCount[:,1])
            DegreeCount[:,1]=np.log10(DegreeCount[:,1]/total)
            DegreeCount[:,0]=np.log10(DegreeCount[:,0])
            slope, intercept, r_value, p_value, std_err = stats.linregress(DegreeCount[:,0],DegreeCount[:,1])
            #print (idegreeCount[:,1])
            #print (len(idegreeCount))
            sns.regplot(x=DegreeCount[:,0],y=DegreeCount[:,1],label='%s - %.3f'%(graph, slope*-1))
        if ylim:
            plt.ylim(ylim,0)
        if xlim:
            plt.xlim(0,xlim)
        plt.legend(fontsize=20)
        plt.xlabel("log(Degree)",fontsize=40)
        plt.xticks(fontsize=20)
        plt.ylabel("log(Relative Frequency)",fontsize=40)
        plt.yticks(fontsize=20)
        if name is None:
            plt.savefig("all_log_degree_distribution.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
    
    
        
            
                    ########################################
                    #                                      #
                    #               STRUCTUALLY            #
                    #           KNOWN and PREDICTED        # 
                    #                NODES                 #
                    #                                      #
                    ########################################       
        
        
    def get_structurally_known_interactions(self,graph,format="networkx"):
        insider_pd=pd.read_csv("Data/Human_insider_interactome_in_PDB.txt", sep=",").dropna()
        DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(graph))
        if len(DF.columns)>2:
            attribute=DF.columns[2]
        else:
            attribute=False
        intersection1=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
        intersection2=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
        #print (intersection1.head())
        if attribute:
            if format=="networkx":     
                intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2',attribute]]
                return nx.from_pandas_edgelist(intersect,"UniprotID_1", "UniprotID_2",edge_attr=attribute,create_using=nx.Graph)
            elif format=="Pandas.DataFrame":
                return pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2',attribute]]
        else:
            if format=="networkx": 
                intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2']]
                return nx.from_pandas_edgelist(intersect,"UniprotID_1", "UniprotID_2",create_using=nx.Graph)
            elif format=="Pandas.DataFrame":
                return pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2']]
        
        
    def get_coverage_of_structurally_known_interactions(self):
        coverage=pd.DataFrame(columns=['Databases','Known_Interactions','Coverage_of_PDB'])
        
        insider_pd=pd.read_csv("Data/Human_insider_interactome_in_PDB.txt", sep=",").dropna()
        len_known_PPIs=len(insider_pd.index)
        print (f'known Human PPIs counts in Insider :  {len_known_PPIs}')
        for name in self.get_name_list():
            DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
            intersection1=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
            intersection2=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
            intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)
            len_intersect=len(intersect.index)
            cover=len_intersect/len_known_PPIs
            coverage=coverage.append({'Databases':name,'Known_Interactions':len_intersect,'Coverage_of_PDB':round(cover,2)}, ignore_index=True)
        return coverage
    
    def visualize_coverage_of_structurally_known_interactions(self,name=None):
        DF=self.get_coverage_of_structurally_known_interactions()
        #pd.p
        #print (DF)
        
        plt.figure(figsize=(30,15), dpi= 100)
        ax=sns.barplot(DF['Databases'],DF['Coverage_of_PDB'])
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel('Coverage_of_PDB',fontsize=40)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')
        if name is None:
            plt.savefig("coverage_of_structurally_known_interactions.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()

        
    def get_structurally_predicted_interactions(self,name,format="networkx"):
        insider_pd=pd.read_csv("Data/Human_insider_interactome_in_I3D.txt", sep=",").dropna()
        DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
        if len(DF.columns)>2:
            attribute=DF.columns[2]
        else:
            attribute=False
        intersection1=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
        intersection2=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
        #print (intersection1.head())
        if attribute:
            if format=="networkx":     
                intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2',attribute]]
                return nx.from_pandas_edgelist(intersect,"UniprotID_1", "UniprotID_2",edge_attr=attribute,create_using=nx.Graph)
            elif format=="Pandas.DataFrame":
                return pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2',attribute]]
        else:
            if format=="networkx": 
                intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2']]
                return nx.from_pandas_edgelist(intersect,"UniprotID_1", "UniprotID_2",create_using=nx.Graph)
            elif format=="Pandas.DataFrame":
                return pd.concat([intersection1,intersection2],ignore_index=True,sort=False)[['UniprotID_1','UniprotID_2']]
        
        
    def get_coverage_of_structurally_predicted_interactions(self):
        coverage=pd.DataFrame(columns=['Databases','Predicted_Interactions','Coverage_of_I3D'])
        
        insider_pd=pd.read_csv("Data/Human_insider_interactome_in_I3D.txt", sep=",").dropna()
        len_known_PPIs=len(insider_pd.index)
        print (f'predicted Human PPIs counts in Insider :  {len_known_PPIs}')
        for name in self.get_name_list():
            DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
            intersection1=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
            intersection2=DF.merge(insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
            intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)
            len_intersect=len(intersect.index)
            cover=len_intersect/len_known_PPIs
            coverage=coverage.append({'Databases':name,'Predicted_Interactions':len_intersect,'Coverage_of_I3D':round(cover,2)}, ignore_index=True)
        return coverage

    def visualize_coverage_of_structurally_predicted_interactions(self,name):
        DF=self.get_coverage_of_structurally_predicted_interactions()
        #pd.p
        #print (DF)
        
        plt.figure(figsize=(30,15), dpi= 100)
        ax=sns.barplot(DF['Databases'],DF['Coverage_of_I3D'])
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel('Coverage_of_I3D',fontsize=40)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')
        if name is None:
            plt.savefig("coverage_of_structurally_predicted_interactions.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
    def get_coverage_structurally_known_and_predicted_interactions(self):
        coverage=pd.DataFrame(columns=['Databases','Predicted_Interactions','Coverage_of_I3D','Known_Interactions','Coverage_of_PDB'])        
        #predicted
        Pred_insider_pd=DataFrame=pd.read_csv("Data/Human_insider_interactome_in_I3D.txt", sep=",").dropna()
        len_known_PPIs=len(Pred_insider_pd.index)
        print (f'predicted Human PPIs counts in Insider :  {len_known_PPIs}')
        for name in self.get_name_list():
            DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
            intersection1=DF.merge(Pred_insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
            intersection2=DF.merge(Pred_insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
            intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)
            len_intersect=len(intersect.index)
            cover=len_intersect/len_known_PPIs
            coverage=coverage.append({'Databases':name,'Predicted_Interactions':len_intersect,'Coverage_of_I3D':round(cover,2)}, ignore_index=True)
        #known   
        
        Known_insider_pd=pd.read_csv("Data/Human_insider_interactome_in_PDB.txt", sep=",").dropna()
        len_known_PPIs=len(Known_insider_pd.index)
        print (f'known Human PPIs counts in Insider :  {len_known_PPIs}')
        PDB_cover_ratio,Known_interaction_count=[],[]
        
        for name in self.get_name_list():
            DF=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
            intersection1=DF.merge(Known_insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_1','UniprotID_2'])
            intersection2=DF.merge(Known_insider_pd, how='inner',left_on=['source','target'],right_on=['UniprotID_2','UniprotID_1'])
            intersect=pd.concat([intersection1,intersection2],ignore_index=True,sort=False)
            len_intersect=len(intersect.index)
            cover=len_intersect/len_known_PPIs
            #coverage=coverage.append({'Databases':name,'Known_Interactions':len_intersect,'Coverage_of_PDB':round(cover,2)}, ignore_index=True)
            PDB_cover_ratio.append(round(cover,2))
            Known_interaction_count.append(len_intersect)
        coverage['Coverage_of_PDB']=PDB_cover_ratio
        coverage['Known_Interactions']=Known_interaction_count
        return coverage
    
    
    
    def visualize_coverages_both_known_and_predicted_interactions(self,percentage=False,name=None):
        if percentage:
            DF=self.get_coverage_structurally_known_and_predicted_interactions()
            DF=DF[['Databases','Coverage_of_PDB','Coverage_of_I3D']]
            DF['Coverage_of_PDB']=DF['Coverage_of_PDB'].astype(float)
            DF['Coverage_of_I3D']=DF['Coverage_of_I3D'].astype(float)
            DF_melted=DF.melt('Databases',var_name='a', value_name='The number')
            plt.figure(figsize=(30,15), dpi= 300)
            ax=sns.barplot(x='Databases', y='The number', hue='a', data=DF_melted)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90,fontsize=40)
            #plt.xlabel(fontsize="xx-large")
            plt.yticks(fontsize=20)
            plt.ylabel('The coverage',fontsize=40)
            plt.xlabel("")
            plt.legend(fontsize=30)
            for p in ax.patches:
                ax.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=20,fontweight='bold')
            if name is None:
                plt.savefig("coverages_both_known_and_predicted_interactions_percentage.svg", format='svg')
            else:
                plt.savefig("%s_percentage.svg"%name, format='svg')
            return plt.show()
            
        else:
            DF=self.get_coverage_structurally_known_and_predicted_interactions()
            DF=DF[['Databases','Known_Interactions','Predicted_Interactions']]
            DF_=pd.DataFrame()
            DF_["Databases"]=DF["Databases"]
            DF_["Structurally_Known_Interactions - PDB"]=DF['Known_Interactions'].astype(int)
            DF_['Structurally_Predicted_Interactions - I3D']=DF['Predicted_Interactions'].astype(int)
            DF_melted=DF_.melt('Databases',var_name='a', value_name='The number')
            plt.figure(figsize=(30,15), dpi= 300)
            ax=sns.barplot(x='Databases', y='The number', hue='a', data=DF_melted)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90,fontsize=40)
            #plt.xlabel(fontsize="xx-large")
            plt.yticks(fontsize=20)
            plt.ylabel('The number',fontsize=40)
            plt.xlabel("")
            plt.legend(fontsize=30)
            for p in ax.patches:
                ax.annotate(format(p.get_height(), '.0f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=20,fontweight='bold')
            if name is None:
                plt.savefig("coverages_both_known_and_predicted_interactions.svg", format='svg')
            else:
                plt.savefig("%s.svg"%name, format='svg')
            return plt.show()
            
             
                    ########################################
                    #                                      #
                    #             PUBLICATION              #
                    #                                      #
                    ########################################            
            
    def visualize_publication_distributions(self,name=None):
        with open ("Data/UniProt_proteome_PubMedIDs.txt") as f:
            lines=f.readlines()[1:]
            Publication_count_Dict={}
        for line in lines:
            items=line.split("\t")
            Publication_count_Dict[items[0]]=len(items[2].split(";"))
            
        Counts=sorted(list(Publication_count_Dict.values()), reverse=True)
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.distplot(Counts)
        plt.xlabel("Publication_Count",fontsize=40)
        plt.xticks(fontsize=20)
        plt.ylabel("Relative Frequency",fontsize=40)
        plt.yticks(fontsize=20)
        if name is None:
            plt.savefig("publication_distribution.svg", format='svg')
        else:
            plt.savefig("%s_distribution.svg"%name, format='svg')
        plt.show()
        
        plt.figure(figsize=(30,15), dpi= 300)
        Pub_Count = np.array(list(collections.Counter(Counts).items()))
        Pub_Count=Pub_Count.astype(float)
        total=sum(Pub_Count[:,1])
        Pub_Count[:,1]=np.log10(Pub_Count[:,1]/total)
        Pub_Count[:,0]=np.log10(Pub_Count[:,0])
        slope, intercept, r_value, p_value, std_err = stats.linregress(Pub_Count[:,0],Pub_Count[:,1])
        #print (idegreeCount[:,1])
        #print (len(idegreeCount))
        sns.regplot(x=Pub_Count[:,0],y=Pub_Count[:,1],label='%s - %.3f'%('publication', slope*-1))
        
        
        
        plt.legend(fontsize=20)
        plt.xlabel("log(Publication_Count)",fontsize=40)
        plt.xticks(fontsize=20)
        plt.ylabel("log(Relative Frequency)",fontsize=40)
        plt.yticks(fontsize=20)
        
        if name is None:
            plt.savefig("publication_log_distribution.svg", format='svg')
        else:
            plt.savefig("%s_log_distribution.svg"%name, format='svg')
        plt.show()
                
            
    def get_correlations_between_degree_and_publication(self):
        with open ("Data/UniProt_proteome_PubMedIDs.txt") as f:
            lines=f.readlines()[1:]
            Publication_count_Dict={}
            for line in lines:
                items=line.split("\t")
                Publication_count_Dict[items[0]]=len(items[2].split(";"))
        names=self.get_name_list()
        D_F=pd.DataFrame(columns=['Databases','Correlation_Coefficient','Node_Count','p-Value'])
        for name in names:
            network=self.get_network(name)
            nodes=list(network.nodes())  
            Degrees=dict(network.degree())
            Degree_Publications={}
            for node in nodes:
                if node in Publication_count_Dict.keys():
                    Degree_Publications[node]=(np.log10(Degrees[node]),np.log10(Publication_count_Dict[node])) #log<
                #else:
                    #Degree_Publications[node]=(np.log10(Degrees[node]),0) #log
            D_P=np.array(list(Degree_Publications.values()),dtype=object)
            coeff,p=stats.pearsonr(D_P[:,0],D_P[:,1])
            D_F=D_F.append({'Databases':name,'Correlation_Coefficient':coeff,'Node_Count':len(D_P),'p-Value':p},ignore_index=True )
        return D_F
            
            
            
        
        
        
        
        
        
    def visualize_degree_vs_publication(self,graph,name=None,color=None):
        network=self.get_network(graph)
        nodes=list(network.nodes())        
        Degrees=dict(network.degree())
        with open ("Data/UniProt_proteome_PubMedIDs.txt") as f:
            lines=f.readlines()[1:]
            Publication_count_Dict={}
            for line in lines:
                items=line.split("\t")
                Publication_count_Dict[items[0]]=len(items[2].split(";"))
        
        Degree_Publications={}
        for node in nodes:
            if node in Publication_count_Dict.keys():
                Degree_Publications[node]=(np.log10(Degrees[node]),np.log10(Publication_count_Dict[node])) 
        D_P=np.array(list(Degree_Publications.values()),dtype=object)
        DF=pd.DataFrame(D_P,columns=['Degree', 'Publication']) 
        DF["Degree"]=DF["Degree"].astype('float64')
        DF["Publication"]=DF["Publication"].astype('float64')
        fig=plt.figure(figsize=(30,30),dpi=300)
        sns.set(rc={"xtick.labelsize":30,"ytick.labelsize":30})
        ax=sns.jointplot(x=DF["Degree"],y=DF['Publication'],kind='kde',color=color,height=10,cbar=True)          
        ax.set_axis_labels("%s log(Degree)"%graph,"log(Publication_count)",fontsize=40)  
        
        #plt.title(name,fontsize=30,fontweight='bold')
        if name is None:
            plt.savefig("%s_degree_vs_publication.svg"%graph, format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
        
            
            
                    ########################################
                    #                                      #
                    #                KEGG                  #
                    #                                      #
                    ########################################        
        
        
    def read_pathway_sif(self,file,sep=","):
        self.pathway_pd=pd.read_csv(file, sep=sep).dropna()
        return self.pathway_pd


    def get_common_edges_between_database_and_pathway(self,database_name,pathway_pd=None):
        if pathway_pd==None:
            pathway_pd=self.pathway_pd
        database_pd=nx.convert_matrix.to_pandas_edgelist(self.get_network(database_name))
        database_pd=database_pd.rename(columns={"source":'UniprotID_1',"target":'UniprotID_2'})
        intersection1=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_1','UniprotID_2'])
        intersection2=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_2','UniprotID_1'])
        intersection2=intersection2.drop(columns=['UniprotID_1_y','UniprotID_2_y'])
        intersection2=intersection2.rename(columns={'UniprotID_1_x':"UniprotID_1",'UniprotID_2_x':"UniprotID_2"})     
        DataFrame_Intersection=intersection1.append(intersection2)
        DataFrame_Intersection=DataFrame_Intersection[['UniprotID_1','UniprotID_2']]
        DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()        
        return DataFrame_Intersection
    
    def get_coverage_between_database_and_pathway(self,database_name,pathway_pd=None):
        if pathway_pd==None:
            pathway_pd=self.pathway_pd
        database_pd=nx.convert_matrix.to_pandas_edgelist(self.get_network(database_name))
        database_pd=database_pd.rename(columns={"source":'UniprotID_1',"target":'UniprotID_2'})
        intersection1=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_1','UniprotID_2'])
        intersection2=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_2','UniprotID_1'])
        intersection2=intersection2.drop(columns=['UniprotID_1_y','UniprotID_2_y'])
        intersection2=intersection2.rename(columns={'UniprotID_1_x':"UniprotID_1",'UniprotID_2_x':"UniprotID_2"})     
        DataFrame_Intersection=intersection1.append(intersection2)
        DataFrame_Intersection=DataFrame_Intersection[['UniprotID_1','UniprotID_2']]
        DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()        
        return round(len(DataFrame_Intersection.index)/len(pathway_pd.index),4)
    
    
    def get_coverages_for_all_pathways(self, database_name):    
        database_pd=nx.convert_matrix.to_pandas_edgelist(self.get_network(database_name))
        database_pd=database_pd.rename(columns={"source":'UniprotID_1',"target":'UniprotID_2'})
        
        path_pds,cov_score,pathways,path_names=[],[],[],[]        
        data=pd.read_csv("Data/KEGG_hsa.txt",sep="\t",names=['KEGG_ID', "Path_name"])
        data['KEGG_ID']=data['KEGG_ID'].str.split(":",expand=True)[1]
        data['Path_name']=data['Path_name'].str.split(" - ",expand=True)[0]
        KEGG_converter=dict(zip(data['KEGG_ID'],data['Path_name']))
        
        KEGG_sifs=os.listdir("Data/SIF_UniProt")
        for sif in KEGG_sifs:
            if sif.endswith(".sif"):
                kegg_pd=pd.read_csv("Data/SIF_UniProt/%s"%sif, sep=",").dropna()
                if len(kegg_pd.index)==0:
                    #print(sif, "there is not any edge constructed by PPI")
                    continue
                path_pds.append(kegg_pd)
                s_name=sif.split("_")[1][:-4]
                pathways.append(s_name)
                path_names.append(KEGG_converter[s_name])
                #print (s_name,len(kegg_pd.index))
        for pathway_pd in path_pds:
            intersection1=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_1','UniprotID_2'])
            intersection2=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_2','UniprotID_1'])
            intersection2=intersection2.drop(columns=['UniprotID_1_y','UniprotID_2_y'])
            intersection2=intersection2.rename(columns={'UniprotID_1_x':"UniprotID_1",'UniprotID_2_x':"UniprotID_2"})     
            DataFrame_Intersection=intersection1.append(intersection2)
            DataFrame_Intersection=DataFrame_Intersection[['UniprotID_1','UniprotID_2']]
            DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()        
            cov_score.append(round(len(DataFrame_Intersection.index)/len(pathway_pd.index),4))
        Res_DataFrame=pd.DataFrame()
        Res_DataFrame['KEGG_ID']=pathways
        Res_DataFrame['Pathway']=path_names
        Res_DataFrame['%s_Coverage'%database_name]=cov_score
        print (len(path_pds),len(pathways),len(cov_score))
        return Res_DataFrame
    
    
    def get_all_coverages_for_all_pathways(self):   
        
        names=self.get_name_list()
        
        path_pds,pathways,path_names=[],[],[]      
        data=pd.read_csv("Data/KEGG_hsa.txt",sep="\t",names=['KEGG_ID', "Path_name"])
        data['KEGG_ID']=data['KEGG_ID'].str.split(":",expand=True)[1]
        data['Path_name']=data['Path_name'].str.split(" - ",expand=True)[0]
        KEGG_converter=dict(zip(data['KEGG_ID'],data['Path_name']))
        
        KEGG_sifs=os.listdir("Data/SIF_UniProt")
        for sif in KEGG_sifs:
            if sif.endswith(".sif"):
                kegg_pd=pd.read_csv("Data/SIF_UniProt/%s"%sif, sep=",").dropna()
                if len(kegg_pd.index)<=30:
                    #print(sif, "there are less than 10 edges constructed by PPI")
                    continue
                path_pds.append(kegg_pd)
                s_name=sif.split("_")[1][:-4]
                pathways.append(s_name)
                path_names.append(KEGG_converter[s_name])
                #print (s_name,len(kegg_pd.index))
            
        
        Res_DataFrame=pd.DataFrame()
        Res_DataFrame['KEGG_ID']=pathways
        Res_DataFrame['Pathway']=path_names 
        
        for name in names:
            database_pd=nx.convert_matrix.to_pandas_edgelist(self.get_network(name))
            database_pd=database_pd.rename(columns={"source":'UniprotID_1',"target":'UniprotID_2'})        
            cov_score=[]
            for pathway_pd in path_pds:
                intersection1=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_1','UniprotID_2'])
                intersection2=database_pd.merge(pathway_pd, how='inner',left_on=['UniprotID_1','UniprotID_2'],right_on=['UniprotID_2','UniprotID_1'])
                intersection2=intersection2.drop(columns=['UniprotID_1_y','UniprotID_2_y'])
                intersection2=intersection2.rename(columns={'UniprotID_1_x':"UniprotID_1",'UniprotID_2_x':"UniprotID_2"})     
                DataFrame_Intersection=intersection1.append(intersection2)
                DataFrame_Intersection=DataFrame_Intersection[['UniprotID_1','UniprotID_2']]
                DataFrame_Intersection=DataFrame_Intersection.drop_duplicates()        
                cov_score.append(round(len(DataFrame_Intersection.index)/len(pathway_pd.index),4))
            Res_DataFrame[name]=cov_score
        #print (len(path_pds),len(pathways),len(cov_score))
        return Res_DataFrame
    
    def visualize_all_coverages_for_all_pathways(self,violin_plot=False,name=None):  
        if violin_plot:
            return self.__violinplot_all_coverages_for_all_pathways(name=name)
             
            
        DataFrame=self.get_all_coverages_for_all_pathways()
        names=self.names
        Average={}
        for graph in names:
            DataFrame[graph]=DataFrame[graph].astype(float)
            Average[graph]=DataFrame[graph].mean()
        
        Sorted_names =[i for i,x in sorted(Average.items(), key=lambda x: x[1], reverse=True)]
        #print (DataFrame)
        #print ("\n"*3,Average,"\n",Sorted_names)
        DataFrame=DataFrame.sort_values(by=Sorted_names[0])
        Pathways=DataFrame['Pathway']
        #print (list(Pathways)[-1],list(Pathways)[0])
        #print (DataFrame,list(Pathways)[-1],"highest",list(Pathways)[0],"\ngraph",list(Pathways))
        DF=DataFrame[Sorted_names]
        #print ("\n",Sorted_names[0])
        
        DF=DF.sort_values(by=Sorted_names[0])
        print (DF)
        np_array=DF.reset_index().values[:,1:].astype(float)
        
        fig=plt.figure(figsize=(30,50),dpi=600)
        ax=sns.heatmap(np_array, cmap="Blues",xticklabels=True, yticklabels=True)
        ax.set_xticklabels(Sorted_names,fontsize=40,fontweight='bold',rotation=90)
        ax.set_yticklabels(Pathways,fontsize=15)
        if name is None:
            plt.savefig("all_coverages_for_all_pathways", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
        
        
    def __violinplot_all_coverages_for_all_pathways(self,name=None):  
        DataFrame=self.get_all_coverages_for_all_pathways()
        names=self.names
        Average={}
        for graph in names:
            DataFrame[graph]=DataFrame[graph].astype(float)
            Average[graph]=DataFrame[graph].mean()
        print (Average)
        Sorted_names =[i for i,x in sorted(Average.items(), key=lambda x: x[1], reverse=True)]
        #print (DataFrame)
        #print ("\n"*3,Average,"\n",Sorted_names)
        DataFrame=DataFrame.sort_values(by=Sorted_names[0])
        Pathways=DataFrame['Pathway']
        DF=DataFrame[Sorted_names]
        
        DF=DF.sort_values(by=Sorted_names[0]).melt(var_name='Databases', value_name='Scores')     
        fig=plt.figure(figsize=(30,20),dpi=300)
        ax=sns.violinplot(x='Databases',y='Scores', data=DF,scale='width',inner="point",cut=0)
        ax.set_xticklabels(Sorted_names,fontsize=40,fontweight='bold',rotation=90)
        #ax.set_yticklabels(ax.get_yticklabels(),fontsize=20,fontweight='bold')
        plt.ylabel("Coverage",fontsize=30)
        plt.yticks(fontsize=30)
        plt.xlabel('')
        if name is None:
            plt.savefig("all_coverages_for_all_pathways", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
        
        
        
            
                    ########################################
                    #                                      #
                    #                COSMIC                #
                    #                                      #
                    ########################################    
                    
    def visualize_cancer_gene_census_publication_distribution(self,name=None):
            
        with open ("Data/COSMIC_Cancer_Gene_Census_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        with open ("Data/UniProt_proteome_PubMedIDs.txt") as f:
            lines=f.readlines()[1:]
            names,counts,cgc_counts=[],[],[]
        for line in lines :
            items=line.split("\t")
            names.append(items[0])
            counts.append(len(items[2].split(";")))
            if items[0] in Census:
                cgc_counts.append(len(items[2].split(";")))
                
        print ("CGC Publications Mean   :",np.mean(cgc_counts),"stand dev", stats.tstd(cgc_counts))
       
        print("All protein Publication mean :",np.mean(counts),"stand dev:" , stats.tstd(counts))
        fig=plt.figure(figsize=(15,10),dpi=300)
        sns.distplot(np.array(counts),kde=True,hist=False, label='All Protein Publications in PubMED')    
        sns.distplot(cgc_counts,kde=True,hist=False, label="Publication count of CGC in PubMED")
        plt.xlim(xmin=-10,xmax=100)    
        plt.xlabel('Publication Counts',fontsize=30)
        plt.ylabel('Density',fontsize=30)
        ttest=stats.ttest_ind(counts,cgc_counts)
        plt.legend(['Publication count of overall proteins in PubMED',"Publication count of CGC in PubMED"],fontsize=15)
        print ("t-test results   :",ttest)
        if name is None:
            plt.savefig("intOGen_CGC_publication_distribution", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show() 
        
        

        
            
        
            
            
            
                    
    def visualize_cancer_gene_census_degree(self,name=None):
        with open ("Data/COSMIC_Cancer_Gene_Census_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        df=pd.DataFrame()
        for graph in self.names:
            network=self.get_network(graph)
            network_nodes=network.nodes()
            temp_df=pd.DataFrame()
            temp_df["Degrees"]=[network.degree(gene) for gene in Census if gene in network_nodes]
            temp_df['Intercatomes']=graph
            df=pd.concat([df,temp_df])
            # sns.boxplot(list(degrees.values()))
            # plt.xlim(0,120)
            # plt.show()
        # print (df.head())
        fig=plt.figure(figsize=(15,10),dpi=300)
        ax=sns.boxplot(x=df['Intercatomes'],y=df['Degrees'], data=df)
        ax.set_xticklabels(self.names,fontsize=40,fontweight='bold',rotation=90)
        #ax.set_yticklabels(ax.get_yticklabels(),fontsize=20,fontweight='bold')
        plt.ylim(0,150)
        plt.xlabel('Publication Counts',fontsize=30)
        plt.ylabel('Density',fontsize=30)
        plt.yticks(fontsize=30)
        if name is None:
            plt.savefig("Cancer_Gene_Census_degree", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
                    
                    
    def get_coverage_cancer_gene_census(self, name):
        with open ("Data/COSMIC_Cancer_Gene_Census_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        network=self.get_network(name)
        Nodes=network.nodes
        c=0
        for protein in Census:
            if protein in Nodes:
                c+=1
        return (c/len(Census))
    
    def get_all_coverages_for_cancer_gene_census(self):
        with open ("Data/COSMIC_Cancer_Gene_Census_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        Census_len=len(Census)
        Res={}
        names=self.get_name_list()
        for name in names:
            network=self.get_network(name)
            Nodes=network.nodes
            c=0
            for protein in Census:
                if protein in Nodes:
                    c+=1
            Res[name]=c/Census_len
        return Res
    
    def visualize_coverages_for_cancer_gene_census(self,name=None):
        DB_coverage=self.get_all_coverages_for_cancer_gene_census()
        databases,coverage=list(DB_coverage.keys()),list(DB_coverage.values())               
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.barplot(databases,coverage)
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel('Coverage_of_Cancer_Gene_Census',fontsize=30)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')
        if name is None:
            plt.savefig("coverages_for_cancer_gene_census.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
            
            
        

  


    

        
        
        
            
                    ########################################
                    #                                      #
                    #             intOGen                  #
                    #                                      #
                    ########################################     
                    
    def visualize_intOgen_publication_distribution(self,name=None):
            
        with open ("Data/intOGen_Drivers_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        with open ("Data/UniProt_proteome_PubMedIDs.txt") as f:
            lines=f.readlines()[1:]
            names,counts,idg_counts=[],[],[]
        for line in lines :
            items=line.split("\t")
            names.append(items[0])
            counts.append(len(items[2].split(";")))
            if items[0] in Census:
                idg_counts.append(len(items[2].split(";")))
                
        print ("intOGen Driver Genes Publications Mean   :",np.mean(idg_counts),"stand dev", stats.tstd(idg_counts))
       
        print("All protein Publication mean :",np.mean(counts),"stand dev:" , stats.tstd(counts))
        fig=plt.figure(figsize=(15,10),dpi=300)
        sns.distplot(np.array(counts),kde=True,hist=True, label='All Protein Publications in PubMED')    
        sns.distplot(idg_counts,kde=True,hist=True, label="Publication count of intOGen Driver Genes in PubMED")
        # plt.xlim(xmin=-10,xmax=100)    
        plt.xlabel('Publication Counts',fontsize=30)
        plt.ylabel('Density',fontsize=30)
        ttest=stats.ttest_ind(counts,idg_counts)
        plt.legend(['Publication count of overall proteins in PubMED',"Publication count of intOGen Driver Genes in PubMED"],fontsize=15)
        print ("t-test results   :",ttest)
        if name is None:
            plt.savefig("intOGen_Driver_Genes_publication_distribution.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()                   
                    
                    
                    
                    
    def visualize_intOgen_driver_genes_degree(self,name=None):
        with open ("Data/intOGen_Drivers_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        df=pd.DataFrame()
        for graph in self.names:
            network=self.get_network(graph)
            network_nodes=network.nodes()
            temp_df=pd.DataFrame()
            temp_df["Degrees"]=[network.degree(gene) for gene in Census if gene in network_nodes]
            temp_df['Intercatomes']=graph
            df=pd.concat([df,temp_df])
            # sns.boxplot(list(degrees.values()))
            # plt.xlim(0,120)
            # plt.show()
        # print (df.head())
        fig=plt.figure(figsize=(30,20),dpi=300)
        ax=sns.boxplot(x=df['Intercatomes'],y=df['Degrees'], data=df)
        ax.set_xticklabels(self.names,fontsize=40,fontweight='bold',rotation=90)
        #ax.set_yticklabels(ax.get_yticklabels(),fontsize=20,fontweight='bold')
        plt.ylim(0,1500)
        plt.ylabel('Degrees',fontsize=30)
        plt.xlabel('intOgen driver genes degree',fontsize=30)
        plt.yticks(fontsize=30)
        if name is None:
            plt.savefig("intOgen_driver_genes_degree", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()                    
                    
    def get_coverage_intOGen_Drivers(self, name):
        with open ("Data/intOGen_Drivers_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        network=self.get_network(name)
        Nodes=network.nodes
        c=0
        for protein in Census:
            if protein in Nodes:
                c+=1
        return (c/len(Census))
    
    def get_all_coverages_for_intOGen_Drivers(self):
        with open ("Data/intOGen_Drivers_with_UniProt.txt","r") as f:
            Census=f.readline().split(",")
        Census_len=len(Census)
        Res={}
        names=self.get_name_list()
        for name in names:
            network=self.get_network(name)
            Nodes=network.nodes
            c=0
            for protein in Census:
                if protein in Nodes:
                    c+=1
            Res[name]=c/Census_len
        return Res
    
    def visualize_coverages_for_intOGen_Drivers(self,name=None):
        DB_coverage=self.get_all_coverages_for_intOGen_Drivers()
        databases,coverage=list(DB_coverage.keys()),list(DB_coverage.values())               
        plt.figure(figsize=(30,15), dpi= 300)
        ax=sns.barplot(databases,coverage)
        plt.xticks(rotation=90,fontsize=40)
        #plt.xlabel(fontsize="xx-large")
        plt.yticks(fontsize=20)
        plt.xlabel('')
        plt.ylabel('Coverage_of_intOGen_Drivers',fontsize=30)
        for p in ax.patches:
            ax.annotate(format(p.get_height(), '.2f'), (p.get_x() + p.get_width() / 2., p.get_height()), ha = 'center', va = 'baseline',fontsize=30,fontweight='bold')
        if name is None:
            plt.savefig("coverages_for_intOGen_Drivers.svg", format='svg')
        else:
            plt.savefig("%s.svg"%name, format='svg')
        return plt.show()
            
            
        
# =============================================================================
# 
# =============================================================================
  
if __name__ == "__main__":
    print ("main")
    GC=GraphCom()
    GC.add_graph_from_CSV("../Interactomes/PathwayCommons_with_UniProt.txt",'PathwayCommons')  
    
    GC.add_graph_from_CSV("../Interactomes/ConsensusPATH_with_UniProt.txt",'ConsensusPath',sep=",",edge_attr='confidence')
    
    GC.add_graph_from_CSV("../Interactomes/Hippie_with_UniProt.txt",'HIPPIE',sep=",",edge_attr='confidence')
    
#    GC=GraphCom()  
#    GC.visualize_intOgen_publication_distribution()
    # GC.add_graph_from_CSV("../Data/Interactomes/iREF_with_UniProt.txt",'iREF',sep=",",edge_attr='confidence')
    
    # GC.add_graph_from_CSV("../Data/Interactomes/Hippie_with_UniProt.txt",'HIPPIE',sep=",",edge_attr='confidence')
    
    # GC.add_graph_from_CSV("../Data/Interactomes/ConsensusPATH_with_UniProt.txt",'ConsensusPath',sep=",",edge_attr='confidence')
    
    # GC.add_graph_from_CSV("../Data/Interactomes/STRING_with_UniProt.txt","STRING",sep=",",edge_attr="calculated_combined_score")
    
    # GC.add_graph_from_CSV("../Data/Interactomes/OmniPath_with_UniProt.txt",'OmniPath', sep=",")
    
    # GC.add_graph_from_CSV("../Data/Interactomes/PathwayCommons_with_UniProt.txt",'PathwayCommons')
        
    # GC.visualize_confidence_scores_correlation_matrix()
    
    # GC.visualize_driver_genes_degree_distributions()
    # GC.visualize_driver_genes_degree_distributions('HIPPIE')
    # GC.visualize_driver_genes_degree_distributions('ConsensusPath')
    # GC.visualize_intOgen_driver_genes_degree_distributions()
    