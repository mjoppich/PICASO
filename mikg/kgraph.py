import networkx as nx


from load_utils import *




class KGraph:
    
    def __init__(self, data_dir) -> None:
        
        self.data_dir = data_dir        
        self.kg = nx.DiGraph()
        
        
    def load_kgraph_base(self, go=True, omnipath=True, opentargets=True, reactome=True, STRING=True):
        
        if go:
            load_go(self.kg, self.data_dir)
            
        if omnipath:
            load_omnipath(self.kg, self.data_dir)
            
        if opentargets:
            load_opentargets(self.kg, self.data_dir)
            
        if STRING:
            load_STRING(self.kg, self.data_dir)
    
    
    