

import networkx as nx
import logging
import pickle

from .load_utils import *


class KGraph:
    
    def __init__(self) -> None:
        self.kg = nx.DiGraph()
        
        self.logger = logging.getLogger("KGraph")
        self.logger.setLevel(logging.INFO)
        
        
    def print_kg_info(self):
        self.logger.info("Current KG {}".format(str(self.kg)))
        
    def load_kgraph_base(self, data_dir, go=True, omnipath=True, opentargets=True, reactome=True, STRING=True):
        
        self.data_dir = data_dir
        
        if go:
            self.logger.info("Loading GeneOntology Graph")
            load_go(self.kg, self.data_dir)
            self.print_kg_info()
            
        if reactome:
            self.logger.info("Loading Reactome Graph")
            load_reactome(self.kg, self.data_dir)
            self.print_kg_info()
            
        if omnipath:
            self.logger.info("Loading OmniPath Graph")
            load_omnipath(self.kg, self.data_dir)
            self.print_kg_info()
            
        if opentargets:
            self.logger.info("Loading OpenTargets Graph")
            load_opentargets(self.kg, self.data_dir)
            self.print_kg_info()
            
        if STRING:
            self.logger.info("Loading STRING Graph")
            load_STRING(self.kg, self.data_dir)
            self.print_kg_info()

    def load_kgraph(self, infile):
        self.kg = nx.read_gpickle(infile)
    
    
    def save_kgraph(self, outfile):
        self.logger.info("Saving KG {}".format(outfile))
        with open(outfile, 'wb') as f:
            pickle.dump(self.kg, f)
    
    
    def get_node_data(self, node):
        
        if not node in self.kg.nodes:
            return None
        
        return self.kg.nodes[node]
    
    def get_node_edges(self, node, in_types=None, out_types=None):
        
        if not node in self.kg.nodes:
            return None
        
        inEdges = [x for x in self.kg.in_edges(node)]
        outEdges = [x for x in self.kg.out_edges(node)]
        
        if not in_types is None:
            inEdges = [x for x in inEdges if self.kg.nodes[x[0]]["type"] in in_types]
        if not out_types is None:
            outEdges = [x for x in outEdges if self.kg.nodes[x[1]]["type"] in out_types]
        
        
        return inEdges + outEdges
    
    
    
    def add_gene_expression(self, exprDF):
        
        gene2expression = {}
        
        for ri, row in exprDF.iterrows():           
            gene2expression[row["gene"]] = { "mean": row["mean"], "perc_expr": row["anum"] / row["num"], "median": row["median"]}
            
        print(len(gene2expression))
        
        foundGenes = set()
                
        for node in self.kg.nodes():           
            if node in gene2expression:
                self.kg.nodes[node]["expression"] = gene2expression[node]
                foundGenes.add(node)
                
        print("Found Genes", len(foundGenes))
        
        
        
    def score_gene_gene_edges(self):
        
        
        self.score_edges()
        
        
    def score_edges(self, edge_accessor, value_accessor, scorers, in_types=None, out_types=None):
        
        
        for edge in kg.get_node_edges("CCL2", in_types=["gene"], out_types=["gene"]):
    
            inExpr = self.get_node_data(edge[0]).get("expression", {}).get("mean", 0)
            outExpr = self.get_node_data(edge[1]).get("expression", {}).get("mean", 0)
            edgeType = self.kg.edges[edge].get("type", "-")
            
            print(edge, edgeType, scoring_gene_gene_expression[edgeType](inExpr, outExpr), inExpr, outExpr)