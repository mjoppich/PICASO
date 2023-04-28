

import networkx as nx
import logging
import pickle

from .load_utils import *
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import seaborn as sns
from itertools import chain

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
    
    
    
    def get_node_types(self):
        
        nodeTypes = Counter()
        for node in self.kg.nodes:
            nodeData = self.kg.nodes[node]
            
            nodeType = nodeData.get("type", None)
            nodeTypes[nodeType] += 1
            
        return nodeTypes
       
    def get_edge_types(self, field="type"):
        
        edgeTypes = Counter()
        for edge in self.kg.edges:
           
            edgeType = self.kg.edges[edge].get(field, None)
            
            edgeTypes[ edgeType ] += 1
            
        return edgeTypes
    
    def get_edge_node_types(self, edge, field="type"):
        srcType = self.kg.nodes[edge[0]].get(field, None)
        tgtType = self.kg.nodes[edge[1]].get(field, None)
        
        return srcType, tgtType
    
    def get_edge_between_type(self):
        
        edgeTypes = Counter()
        for edge in self.kg.edges:           
            srcType, tgtType = self.get_edge_node_types(edge)
            edgeTypes[(srcType, tgtType)] += 1
            
        return edgeTypes
    
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
        
        
        
    def get_edges_between_ntypes(self, in_types, out_types):
        
        returnEdges=set()
        for edge in self.kg.edges:
            
            
            if self.kg.nodes[edge[0]].get("type", None) is None:
                print(edge, self.kg.edges[edge])
                print(edge[0], "is Nonetype")
            if self.kg.nodes[edge[1]].get("type", None) is None:
                print(edge, self.kg.edges[edge])
                print(edge[1], "is Nonetype")      
                
                
            if not in_types is None:
                srcAccept = self.kg.nodes[edge[0]].get("type", None) in in_types
            else:
                srcAccept = True
                
            if not out_types is None:
                tgtAccept = self.kg.nodes[edge[1]].get("type", None) in out_types 
            else:
                tgtAccept=True

            if srcAccept and tgtAccept:
                returnEdges.add(edge)
                
        return returnEdges
        
    def score_gene_gene_edges(self, scoring_gene_gene_expression):       
        self.score_edges(lambda x: x.get("expression", {}).get("mean", 0), "type", scoring_gene_gene_expression, in_types=["gene"], out_types=["gene"])
        
        
    def score_edges(self, value_accessor, edge_accessor, scorers, in_types=None, out_types=None, ignore_edge_types=None):
        
        for edge in self.get_edges_between_ntypes(in_types=in_types, out_types=out_types):
            
            srcType, tgtType = self.get_edge_node_types(edge)
            
            if not ignore_edge_types is None:
                if (srcType, tgtType) in ignore_edge_types:
                    continue
    
            inExpr = value_accessor(self.get_node_data(edge[0]))
            outExpr = value_accessor(self.get_node_data(edge[1]))
            
            edgeType = self.kg.edges[edge].get( edge_accessor , "-")
            
            edgeScore = scorers[edgeType](inExpr, outExpr)
            
            self.kg.edges[edge]["score"] = edgeScore
            
    def score_nodes(self, ntype="geneset", consider_edge_type=[("gene", "geneset")], scoring_function=None):
        
        assert(not scoring_function is None)
        
        inTypes = [x[0] for x in consider_edge_type]
        outTypes = [x[1] for x in consider_edge_type]
        
        for node in self.kg.nodes:
            
            if self.get_node_data(node).get("type", None) == ntype:
                
                edges = self.get_node_edges(node, in_types=inTypes, out_types=outTypes)
                
                nodeScore = scoring_function(self.kg, node, edges)
                
                self.kg.nodes[node]["score"] = nodeScore
                
                
    def get_edge_scores(self, score_accessor=lambda x: x.get("score", 0), edge_types=None):
        
        allScores = []
        
        for iedge, edge in enumerate(self.kg.edges):
            
            if not edge_types is None:
                entype = self.get_edge_node_types(edge)
                if not entype in edge_types:
                    continue
                                
            
            edgeScore = score_accessor(self.kg.edges[edge])
            allScores.append(edgeScore)
            
        return allScores
        

    def get_edge_scores_per_type(self, score_accessor=lambda x: x.get("score", 0), edge_types=None):
        
        allScores = defaultdict(list)
        
        for iedge, edge in enumerate(self.kg.edges):
            
            entype = self.get_edge_node_types(edge)
            if not edge_types is None:
                if not entype in edge_types:
                    continue
                                
            
            edgeScore = score_accessor(self.kg.edges[edge])
            allScores[entype].append(edgeScore)
            
        return allScores        
            
        
    def plot_score_histogram(self, edge_types=None):
        
        scores = self.get_edge_scores(edge_types=edge_types)
        
        fig, ax = plt.subplots(figsize=(8, 4))

        # plot the cumulative histogram
        n, bins, patches = ax.hist(scores, len(scores), density=True, histtype='step',
                                cumulative=True, label='Empirical Score Distribution')

        # tidy up the figure
        ax.grid(True)
        ax.legend(loc='right')
        ax.set_title('Cumulative Histogram of Edge Scores')
        ax.set_xlabel('Score')
        ax.set_ylabel('Likelihood of Score')

        plt.show()
        
    def plot_score_violin(self, per_edge_type=False, edge_types=None):
        
        scores = dict()
        if per_edge_type:
            scores = self.get_edge_scores_per_type(edge_types=edge_types)
        else:
            scores[("all",)] = self.get_edge_scores(edge_types=edge_types)
            
        oldkeys = [x for x in scores]
        for okey in oldkeys:
            values = scores[okey]
            del scores[okey]
            scores[ " -> ".join(okey) ] = values
            
        print([x for x in scores])
            
        keys, values = map(chain.from_iterable, zip(*(([k]*len(v), v) for k, v in scores.items())))
        scoresDF = pd.DataFrame({'class': list(keys), 'value': list(values)})
        
        print(scoresDF.head())
                
        chart = sns.violinplot(data=scoresDF, x="class", y="value")
        chart.set_xticklabels(chart.get_xticklabels(), rotation=45, horizontalalignment='right')

        plt.show()
        
    def score_subgraphs(self, Gs):
        
        scores = defaultdict(list)
        for gname in Gs:
            
            G = Gs[gname]
            
            for edge in G.edges:
                
                escore = G.edges[edge].get("score", 0)
                scores[gname].append(escore)
                
        return scores            
            
            
    def plot_subgraph_scores(self, scores):
        
        import matplotlib.gridspec as grid_spec
        from sklearn.neighbors import KernelDensity
        
        
        df = pd.DataFrame(columns=["group", "score"])
        for x in scores:
            print(x, np.mean(scores[x]))
            for y in scores[x]:
                
                df.loc[len(df)] = (x,y)
        
        
        xlim = df.score.max()
        
        
        groups = [x for x in np.unique(df.group)]
        
        gs = grid_spec.GridSpec(len(groups),1)
        fig = plt.figure(figsize=(16,2*len(groups)))

        colors = plt.cm.viridis(np.linspace(0, 1, len(groups)))

        ax_objs = []
        for i, group in enumerate(groups):

            x = np.array(df[df.group == group].score)
            x_d = np.linspace(0,xlim, 1000)

            kde = KernelDensity(bandwidth=0.2, kernel='gaussian')
            kde.fit(x[:, None])

            logprob = kde.score_samples(x_d[:, None])

            # creating new axes object
            ax_objs.append(fig.add_subplot(gs[i:i+1, 0:]))

            # plotting the distribution
            ax_objs[-1].plot(x_d, np.exp(logprob),color=colors[i],lw=1)
            ax_objs[-1].fill_between(x_d, np.exp(logprob), alpha=1, color=colors[i])


            # setting uniform x and y lims
            ax_objs[-1].set_xlim(0,xlim)
            ax_objs[-1].set_ylim(0,2)

            # make background transparent
            rect = ax_objs[-1].patch
            rect.set_alpha(0)

            # remove borders, axis ticks, and labels
            ax_objs[-1].set_yticklabels([])

            if i == len(groups)-1:
                ax_objs[-1].set_xlabel("Score", fontsize=16,fontweight="bold")
            else:
                ax_objs[-1].set_xticklabels([])

            spines = ["top","right","left","bottom"]
            for s in spines:
                ax_objs[-1].spines[s].set_visible(False)

            adj_country = group.replace(" ","\n")
            ax_objs[-1].text(-0.02,0,adj_country,fontweight="bold",fontsize=14,ha="right")


        gs.update(hspace=-0.5)

        fig.text(0.07,0.85,"RidgePlot of Scores",fontsize=20)

        plt.tight_layout()
        plt.show()
        
    def score_subgraphs_for_subnet(self, KGs, subnet):
        
        kgSubgraphs = {}        
        
        for group in KGs:
            kgSubgraphs[group] = KGs[group].kg.subgraph(subnet)
        
        kgScores = self.score_subgraphs(kgSubgraphs)
        
        return kgScores