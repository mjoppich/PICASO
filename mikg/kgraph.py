

import networkx as nx
import logging
import pickle

from .load_utils import *
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import seaborn as sns
from itertools import chain
import tempfile
import community
import glob
import random
import json

from itertools import product
import leidenalg as la
import igraph as ig

import markov_clustering as mc
import scipy

from adjustText import adjust_text
from sklearn import decomposition
from sklearn import preprocessing
import infomap
from natsort import natsorted

import networkx as nx
import progressbar

# from https://stackoverflow.com/questions/25500541/matplotlib-bwr-colormap-always-centered-on-zero
import matplotlib

class MidpointNormalize(matplotlib.colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        matplotlib.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        
        assert(clip is None)
        
        wasScalar = False
        if isinstance(value, (float, int)):
            value = np.array([value])
            wasScalar = True

        allvalues = []
        for x in value:
            if x < self.midpoint:
                allvalues.append( np.interp([x], [self.vmin, self.midpoint], (0, 0.5))[0] )
            else:
                allvalues.append( 0.5 + np.interp([x], [self.midpoint, self.vmax], (0, 0.5))[0] )

        if not wasScalar:
            return np.array(allvalues)
        
        return allvalues[0]
        
#####

def sjoined(inlist):
    return ";".join([str(x) for x in sorted(set(inlist))])

class KGraph:
    
    def __init__(self, random_state=42, kgraph_name="KGraph") -> None:
        self.kg = nx.DiGraph()
        
        self.random_state=random_state
        random.seed(random_state)
        
        self.kgraph_name = kgraph_name
        
        self.logger = logging.getLogger(kgraph_name)
        self.logger.setLevel(logging.INFO)
               

    def is_score_field(self, x):
        return "score" in x
    
       
    def get_node_attributes(self):
        
        node_attribs = set()
        for n in self.kg.nodes:
            for attrib in self.kg.nodes[n]:
                node_attribs.add(attrib)
                
        return node_attribs
       
    def get_edge_attributes(self):
        
        edge_attribs = set()
        for e in self.kg.edges:
            for attrib in self.kg.edges[e]:
                edge_attribs.add(attrib)
                
        return edge_attribs   
       
    def __repr__(self):
        return "KGraph {} with {} nodes and {} edges".format(self.kgraph_name, len(self.kg.nodes), len(self.kg.edges))
       
    def copy(self, suffix=None):
        
        kgname=self.kgraph_name
        if not suffix is None:
            kgname="{}_{}".format(self.kgraph_name, suffix)
        
        ckg = KGraph(random_state=self.random_state, kgraph_name=kgname)       
        ckg.kg = self.kg.copy()
        
        return ckg
        
    def print_kg_info(self):
        self.logger.info("Current KG {}".format(str(self.kg)))
        
    def load_kgraph_base(self, data_dir, go=True, TFs=True, omnipath=True, opentargets=True, reactome=True, kegg=True, STRING=True, NPINTER=True, ot_min_disease_assoc_score=0.8, hallmark_genesets="kegg_gmts/human/c1.all.v2023.2.Hs.symbols.gmt", curated_genesets="kegg_gmts/human/c2.all.v2023.2.Hs.symbols.gmt"):
        
        self.data_dir = data_dir

        if STRING:
            self.logger.info("Loading STRING Graph")
            load_STRING(self.kg, self.data_dir)
            self.print_kg_info()
        
        if TFs:
            self.logger.info("Loading Human Transcription Factors Graph")
            load_human_transcription_factors(self.kg, self.data_dir)
            self.print_kg_info()
        
        if go:
            self.logger.info("Loading GeneOntology Graph")
            load_go(self.kg, self.data_dir)
            self.print_kg_info()
            
        if reactome:
            self.logger.info("Loading Reactome Graph")
            load_reactome(self.kg, self.data_dir)
            self.print_kg_info()

        if kegg:
            self.logger.info("Loading KEGG Graph")
            load_kegg(self.kg, self.data_dir, hallmark_genesets=hallmark_genesets, curated_genesets=curated_genesets)
            self.print_kg_info()
            
        if omnipath:
            self.logger.info("Loading OmniPath Graph")
            load_omnipath(self.kg, self.data_dir)
            self.print_kg_info()
            
        if opentargets:
            self.logger.info("Loading OpenTargets Graph")
            load_opentargets(self.kg, self.data_dir, min_disease_association_score=ot_min_disease_assoc_score)
            self.print_kg_info()

        if NPINTER:
            #NPINTER should always be loaded last, because then all genes are already added to the graph!
            self.logger.info("Loading NPINTER5 Graph")
            load_npinter(self.kg, self.data_dir)
            self.print_kg_info()
            
        #remove singletons
        remove = [node for node,degree in dict(self.kg.degree()).items() if degree < 1]
        print("Removing {} singletons".format(len(remove)))
        self.kg.remove_nodes_from(remove)

    def load_kgraph(self, infile):
        self.logger.info("Loading KG {}".format(infile))
        with open(infile, 'rb') as f:
            self.kg = pickle.load(f)
    
    
    def save_kgraph(self, outfile):
        self.logger.info("Saving KG {}".format(outfile))
        with open(outfile, 'wb') as f:
            pickle.dump(self.kg, f)
    
    
    def get_node_data(self, node):
        
        if not node in self.kg.nodes:
            return None
        
        return self.kg.nodes[node]
    
    def get_node_edges(self, node, src_types=None, tgt_types=None):
        
        if not node in self.kg.nodes:
            return None
               
        inEdges = [x for x in self.kg.in_edges(node)]
        outEdges = [x for x in self.kg.out_edges(node)]
        
        if not src_types is None:
            src_types = set(src_types)
            inEdges = [x for x in inEdges if len(src_types.intersection(self.kg.nodes[x[0]]["type"])) > 0]
        if not tgt_types is None:
            tgt_types = set(tgt_types)
            outEdges = [x for x in outEdges if len(tgt_types.intersection(self.kg.nodes[x[0]]["type"])) > 0]
        
        return inEdges + outEdges
    
    
    def get_nodes(self, nodetype):
        if nodetype is None:
            return [x for x in self.kg.nodes]
        else:
            if isinstance(nodetype, str):
                nodetype = [nodetype]
            
            return [ x for x in self.kg.nodes if self.node_type_overlap(x, nodetype)]
    
    def get_nodes_with_highest_scores(self, n=3000, nodetype=None, score_accessor=lambda x: x["score"], node_filter=None):
        
        scored = Counter()
        
        for node in self.kg.nodes:
            if not nodetype is None:
                if not nodetype in self.kg.nodes[node].get("type", []):
                    continue
                
            if not node_filter is None:
                if not node_filter(node, self):
                    continue
                                
            scored[node] = score_accessor(self.kg.nodes[node])
            
        return scored.most_common(n)
        
    def get_edges_with_highest_scores(self, n=3000, src_types=None, tgt_types=None, score_accessor=lambda x: x["score"]):
        
        scored = Counter()
        
        if not src_types is None:
            src_types = set(src_types)
        if not tgt_types is None:
            tgt_types = set(tgt_types)
        
        for edge in self.kg.edges:
            
            src = edge[0]
            tgt = edge[1]
            
            if not src_types is None:
                if src_types.intersection(self.kg.nodes[src].get("type", [])) == 0:
                    continue
            if not tgt_types is None:
                if tgt_types.intersection(self.kg.nodes[tgt].get("type", [])) == 0:
                    continue
                                
            scored[edge] = score_accessor(self.kg.edges[edge])
            
        return scored.most_common(n)
        
    def node_types(self, node):
        nodeNodeType = self.kg.nodes[node].get("type", [])
        return nodeNodeType
    
    def get_node_types(self, single=False):
        
        nodeTypes = Counter()
        if not single:
            
            for node in self.kg.nodes:
                nodeData = self.kg.nodes[node]
                
                nodeType = sjoined(nodeData.get("type", []))
                nodeTypes[nodeType] += 1
                
            
        else:
            
            for node in self.kg.nodes:
                nodeData = self.kg.nodes[node]
                nodeNodeType = nodeData.get("type", [])
                
                for nodeType in nodeNodeType:
                    nodeTypes[nodeType] += 1
            
        return nodeTypes
       
    def plot_node_types(self):
        
        ntCounter = self.get_node_types()
        
        counteritems = [("{}\nn={}".format(x, ntCounter[x]), ntCounter[x]) for x in ntCounter]
        counteritems = sorted(counteritems, key=lambda x: x[1])
                   
        self._plot_pie(counteritems)
       
    def plot_node_children(self, node_types):
        
        scores = []
        for x in self.get_nodes(node_types):
            num_Children = len(self.kg.in_edges(x))
            scores.append(num_Children)
        
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
       
    #
    def plot_edge_types(self):
        
        ntCounter = self.get_edge_types()
        
        counteritems = [("{}\nn={}".format(x, ntCounter[x]), ntCounter[x]) for x in ntCounter]
        counteritems = sorted(counteritems, key=lambda x: x[1])
                   
        self._plot_pie(counteritems)
        
    def plot_edge_sources(self):
        
        ntCounter = self.get_edge_types(field="source")
        
        counteritems = [("{}\nn={}".format(x, ntCounter[x]), ntCounter[x]) for x in ntCounter]
        counteritems = sorted(counteritems, key=lambda x: x[1])
                   
        self._plot_pie(counteritems)
       
    def plot_edge_between_types(self, show_threshold=0.02):
        
        ntCounter = self.get_edge_between_type()
        
        counteritems = [("{} →\n{}, n={}".format(x[0], x[1], ntCounter[x]), ntCounter[x]) for x in ntCounter]
        counteritems = sorted(counteritems, key=lambda x: x[1])
        
        totalcount = sum([x[1] for x in counteritems])
        
        filtereditems = []
        other_count = 0
        for x in counteritems:
            
            if (x[1]/totalcount) < show_threshold:
                other_count += x[1]
            else:
                filtereditems.append(x)
                
        if other_count > 0:
            filtereditems.append(("Other, n={}".format(other_count), other_count))                
                   
        self._plot_pie(filtereditems)
        
            
    def _plot_pie(self, counteritems):
        
        counteritems = sorted(counteritems, key=lambda x: x[1])

        finallist = []
        for i in range(0, len(counteritems)):
            if i % 2 == 0:
                finallist.append( counteritems.pop() )
            else:
                finallist.append(counteritems.pop(0))
        
        data = []
        labels = []        
        explode = []
        
        for x,y in finallist:
            data.append(y)
            labels.append(x)
            
            if len(explode) % 2 == 0:
                explode.append(0.01)
            else:
                explode.append(0.1)
        
        # Creating plot
        fig, ax1 = plt.subplots()

        wedges, texts, autotexts = plt.pie(data, labeldistance=1.2,
                                          explode = explode, autopct = '%1.1f%%',
                                          wedgeprops = {"edgecolor" : "black", 'linewidth' : 2,
                                          'antialiased': True}, startangle=25)
        
        for autotext in autotexts:
            autotext.set_color('white')
        
        kw=dict(xycoords='data',textcoords='data',arrowprops=dict(arrowstyle='-'),zorder=0,va='center')

        labels2pos = {}
        for i,p in enumerate(wedges):
            ang=(p.theta2-p.theta1)/2. +p.theta1
            y=np.sin(np.deg2rad(ang))
            x=np.cos(np.deg2rad(ang))
            
            offset = 1.5 + (i%2)*0.5
                        
            horizontalalignment={-1:"right",1:"left"}[int(np.sign(x))]
            connectionstyle="angle,angleA=0,angleB={}".format(ang)
            kw["arrowprops"].update({"connectionstyle":connectionstyle})
            ax1.annotate(labels[i],xy=(x, y),xytext=(offset*np.sign(x),1.3*y),
                        horizontalalignment=horizontalalignment,**kw)


    def get_edge_types(self, field="type"):
        
        edgeTypes = Counter()
        for edge in self.kg.edges:
           
            edgeType = self.kg.edges[edge].get(field, None)
            
            edgeTypes[ edgeType ] += 1
            
        return edgeTypes
    
    def get_edge_node_types(self, edge, field="type"):
        srcTypes = self.kg.nodes[edge[0]].get(field, None)
        tgtTypes = self.kg.nodes[edge[1]].get(field, None)
        
        return srcTypes, tgtTypes


    def get_edge_edge_types(self, src_types=None, tgt_types=None, type_accessor="type", ignore_edge_types=None):

        eetypes = Counter()
        
        for edge in self.get_edges_between_ntypes(src_types=src_types, tgt_types=tgt_types):
            
            srcTypes, tgtTypes = self.get_edge_node_types(edge)
            
            srcType = sjoined(srcTypes)
            tgtType = sjoined(tgtTypes)
                    
            if not ignore_edge_types is None:
                if (srcType, tgtType) in ignore_edge_types:
                    continue
                    
            edgeType = self.kg.edges[edge].get( type_accessor , "-")

            eetypes[(srcType, tgtType, edgeType)] += 1


        return eetypes
    
    
    def get_edge_between_type(self):
        
        edgeTypes = Counter()
        for edge in self.kg.edges:           
            srcTypes, tgtTypes = self.get_edge_node_types(edge)
            
            srcType = sjoined(srcTypes)
            tgtType = sjoined(tgtTypes)
            
            edgeTypes[(srcType, tgtType)] += 1
            
        return edgeTypes
    
    def add_gene_expression(self, exprDF, mean_column="mean", sd_column="sd", perc_expr_column = "perc_expr", median_column="median", num_column="num", allnum_column="group_cells"):
        
        gene2expression = {}
        
        for ri, row in exprDF.iterrows():           
            gene2expression[row["gene"]] = { "mean": row[mean_column],
                                             "sd": row[sd_column],
                                             "perc_expr": row[perc_expr_column],
                                             "median": row[median_column],
                                             "num_measured": row[num_column],
                                             "num_all": row[allnum_column]}
                        
        print("Measured Genes", len(gene2expression))
        
        foundGenes = set()
                
        for node in self.kg.nodes():           
            if node in gene2expression:
                self.kg.nodes[node]["type"].add("measured_expression")
                self.kg.nodes[node]["expression"] = gene2expression[node]
                foundGenes.add(node)
                
        print("Found Genes", len(foundGenes))
        
        
    def get_node_type(self, node):
        
        assert node in self.kg.nodes
        
        return self.kg.nodes[node].get("type", None)
        
    def get_edges_to_type(self, node, otype):
        
        allEdges = []
        for inEdge in self.kg.in_edges(node):
            oNode = inEdge[0]
            
            if otype in self.get_node_type(oNode):
                allEdges.append(oNode)
                
        for outEdge in self.kg.out_edges(node):
            oNode = outEdge[1]
            
            if otype in self.get_node_type(oNode):
                allEdges.append(oNode)
                        
        return allEdges
        
    def get_edges_between_ntypes(self, src_types, tgt_types):
        
        returnEdges=set()
        for edge in self.kg.edges:
            
            
            if self.kg.nodes[edge[0]].get("type", None) is None:
                print(edge, self.kg.edges[edge])
                print(edge[0], "is Nonetype")
            if self.kg.nodes[edge[1]].get("type", None) is None:
                print(edge, self.kg.edges[edge])
                print(edge[1], "is Nonetype")      
                
                
            if not src_types is None:
                src_types = set(src_types)
                srcAccept = len(src_types.intersection(self.kg.nodes[edge[0]].get("type", []))) > 0
            else:
                srcAccept = True
                
            if not tgt_types is None:
                tgt_types = set(tgt_types)
                tgtAccept = len(tgt_types.intersection(self.kg.nodes[edge[1]].get("type", []))) > 0
            else:
                tgtAccept=True

            if srcAccept and tgtAccept:
                returnEdges.add(edge)
                
        return returnEdges
                        
                
    def _get_predecessors(self, start_node, ntype, n=10):
        
        if not isinstance(ntype, (list, set, tuple)):
            ntype = set([ntype])
        
        if self.node_type_overlap(start_node, ntype):
            return set([start_node])
    
        if n == 0:
            return set()
        
        inEdges = self.kg.in_edges(start_node)
        
        returnChildren = set()
        
        for child, curnode in inEdges:
                        
            if self.node_type_overlap(child, ntype):
                returnChildren.add(child)
                
            else:
                if None in ntype:
                    returnChildren.add(child)
                
                recChildren = self._get_predecessors(child, ntype, n-1)
                returnChildren = returnChildren.union(recChildren)
                
        return returnChildren
                
                
    def score_nodes_hierarchically(self, ntype="geneset", target_ntype="gene", relevance_threshold=0, child_score_accessor=lambda x: x.get("score", 0)):
        
        assert(not child_score_accessor is None)
        
        geneset2score = {}

        for node in self.kg.nodes:
            
            if ntype in self.get_node_data(node).get("type", []):
                
                #get all children of node of type target_ntype
                geneChildren = self._get_predecessors(node, target_ntype)
                
                #get all scores for children
                childrenScores = []
                for child in geneChildren:
                    childrenScores.append( child_score_accessor(self.kg.nodes[child]) )
                
                    
                if len(childrenScores) == 0:
                    nodeMedian = 0
                    nodeScore = 0    
                else:
                    nodeScore = np.mean(childrenScores)
                    nodeSD = np.std(childrenScores)
                    nodeMedian = np.median(childrenScores)
                
                if nodeMedian < relevance_threshold:
                    #assume geneset not relevant
                    nodeScore = 0
                    
                self.kg.nodes[node]["score"] = nodeScore
                geneset2score[node] = (nodeScore, nodeMedian)
                    
        return geneset2score
        
    def node_type_overlap(self, node, types):
        if not isinstance(types, (list, tuple, set)):
            types = set([types])
            
        return len(set(types).intersection(self.kg.nodes[node]["type"])) > 0
                
    def get_edge_scores(self, score_accessor=lambda x: x.get("score", 0), edge_types=None, nodes=None):
        
        allScores = []
        
        for iedge, edge in enumerate(self.kg.edges):
            
            if not edge_types is None:
                srcTypes, tgtTypes = self.get_edge_node_types(edge)
                
                etypes = set([x for x in product(srcTypes, tgtTypes)])
                if not len(etypes.intersection(edge_types)) > 0:
                    continue
                  
            if not nodes is None:
                if not (edge[0] in nodes and edge[1] in nodes):
                    continue
            
            edgeScore = score_accessor(self.kg.edges[edge])
            allScores.append(edgeScore)
            
        return allScores
        

    def get_edge_scores_per_type(self, score_accessor=lambda x: x.get("score", 0), edge_types=None, single=False):
        
        allScores = defaultdict(list)
        
        for iedge, edge in enumerate(self.kg.edges):
            
            srcTypes, tgtTypes = self.get_edge_node_types(edge)
            if not edge_types is None:
                
                etypes = set([x for x in product(srcTypes, tgtTypes)])
                if not len(etypes.intersection(edge_types)) > 0:
                    continue
                                
            
            edgeScore = score_accessor(self.kg.edges[edge])
            
            if not single:
                edgeRep = (sjoined(srcTypes), sjoined(tgtTypes))
                allScores[edgeRep].append(edgeScore)
            else:
                etypes = set([x for x in product(srcTypes, tgtTypes)])
                for etype in etypes:
                    allScores[etype].append(edgeScore)
            
        return allScores        
            
        
    def plot_score_histogram(self, edge_types=None, score_accessor=lambda x: x.get("score", 0)):
        
        scores = self.get_edge_scores(edge_types=edge_types, score_accessor=score_accessor)
        
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
        

        
        
    def plot_edge_attribute_distribution(self, edge_types=None, score_accessor=lambda x: x.get("score", 0), ax=None, title=None):
        
        
        scores = self.get_edge_scores(edge_types=edge_types, score_accessor=score_accessor)
     
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))

        if title is None:
            title = 'Score Distribution'
            
        # plot the cumulative histogram
        ax.boxplot( [scores], labels=["Score"] )

        # tidy up the figure
        ax.grid(True)
        ax.legend(loc='right')
        ax.set_title('Boxplot of edge scores')
        ax.set_xlabel('Node type')
        ax.set_ylabel('Score')
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, horizontalalignment='right')
        
        
        
    def plot_node_attribute_distribution(self, attribute_accessor, node_types=None, ax=None, title=None):
        
        
        if node_types == None:
            node_types = [x for x in self.get_node_types(single=True)]
        
        scoreByNodetype = defaultdict(list)
        
        for node in self.kg.nodes:
            if self.node_type_overlap(node, node_types):
                for ntype in self.get_node_type(node):
                    scoreByNodetype[ntype].append( attribute_accessor(self.kg.nodes[node]))
                    
                
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))

        if title is None:
            title = 'Score Distribution'
            
        allNodeTypes = natsorted([x for x in scoreByNodetype])

        # plot the cumulative histogram
        ax.boxplot( [scoreByNodetype[x] for x in allNodeTypes], labels=allNodeTypes )

        # tidy up the figure
        ax.grid(True)
        ax.legend(loc='right')
        ax.set_title('Boxplot of node scores by node type')
        ax.set_xlabel('Node type')
        ax.set_ylabel('Score')
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 45, horizontalalignment='right')
        
        
    def plot_node_attribute_histogram(self, attribute_accessor, nodes=None, node_type=None, ax=None, title=None):
        
        if nodes is None:
            nodes = self.get_nodes(node_type)
        else:
            nodes = set(nodes).intersection(self.kg.nodes)
        
        scores = []
        for node in nodes:
            scores.append( attribute_accessor( self.kg.nodes[node] ) )
            
        if len(scores) == 0:
            self.logger.info("No nodes selected.")
            return
                
        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 4))

        if title is None:
            title = 'Empirical Score Distribution'

        # plot the cumulative histogram
        n, bins, patches = ax.hist(scores, len(scores), density=True, histtype='step',
                                cumulative=True, label=title)

        # tidy up the figure
        ax.grid(True)
        ax.legend(loc='right')
        ax.set_title('Cumulative Histogram of Node Scores')
        ax.set_xlabel('Score')
        ax.set_ylabel('Likelihood of Score')
        
    def plot_score_violin(self, per_edge_type=False, single_edge_types=False, edge_types=None, score_accessor=lambda x: x.get("score", 0)):
        
        scores = dict()
        if per_edge_type:
            scores = self.get_edge_scores_per_type(edge_types=edge_types, score_accessor=score_accessor, single=single_edge_types)
            
            if single_edge_types:
                self.logger.warning("Caution: if an edge's nodes have multiple node types, its score is attributed to the product of all node-types.")
        else:
            scores[("all",)] = self.get_edge_scores(edge_types=edge_types, score_accessor=score_accessor)
            
        oldkeys = [x for x in scores]
        for okey in oldkeys:
            values = scores[okey]
            del scores[okey]
            scores[ " -> ".join(okey) ] = values
                       
        keys, values = map(chain.from_iterable, zip(*(([k]*len(v), v) for k, v in scores.items())))
        scoresDF = pd.DataFrame({'class': list(keys), 'value': list(values)})
                        
        chart = sns.violinplot(data=scoresDF, x="class", y="value")
        chart.set_xticklabels(chart.get_xticklabels(), rotation=45, horizontalalignment='right')

        plt.show()
        
    
            
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
        

    def filter_nodes(self, filter_function):
        assert(not filter_function is None)
        
        retainNodes = list()

        for node in self.kg.nodes:
            if not filter_function(node, self):
                continue
        
            retainNodes.append(node)
            
        return self.subset_kg(retainNodes, suffix="filter")

    def filter_edges(self, filter_function):
        assert(not filter_function is None)
        
        removeEdges = list()

        for edge in self.kg.edges:
            if not filter_function(edge, self):
                removeEdges.append(edge)
                    
        retKG = self.copy(suffix="filter")
        retKG.kg.remove_edges_from(removeEdges)
            
        return retKG


    def subset_kg(self, retainedNodes, suffix="subset"):
        
        ret = KGraph(random_state=self.random_state, kgraph_name="{}_{}".format(self.kgraph_name, suffix))
        ret.kg = nx.subgraph(self.kg, retainedNodes).copy()
    
        return ret
    
    def induce_kg(self, retainedNodes, suffix="subset", radius=1):
        
        allNodes = set()
        for en in retainedNodes:
            allNodes.add(en)
            sg = nx.ego_graph(self.kg, en, radius=radius)
            
            for extnode in sg.nodes:
                allNodes.add(extnode)
                
        return self.subset_kg(allNodes, suffix=suffix)
    
    def to_gene_kgraph(self):
        
        ret = KGraph(random_state=self.random_state, kgraph_name="{}_gene".format(self.kgraph_name))
        ret.kg = nx.subgraph(self.kg, self.get_nodes("gene")).copy()
    
        return ret
    
    def _filter_edge_score(self, score_field, minEdgeScore=3.0, undirected=True):
        
        sub_kg = self.kg.edge_subgraph([x for x in self.kg.edges if ((minEdgeScore is None) or (self.kg.edges[x].get(score_field, 0) > minEdgeScore))]).copy()
        
        if undirected:
            return sub_kg.to_undirected()
        
        return sub_kg

    def get_communities_connectedcomponent(self, minEdgeScore = 3.0, resolution=5, prefix="Module", sep="_", score_field="score"):
        
        sub_kg_ud = self._filter_edge_score(score_field, minEdgeScore=minEdgeScore)
        
        stage_comms = {}
        prefix="cc"
        for ci, cc in enumerate(nx.connected_components(sub_kg_ud)):
            stage_comms["{}{}{}".format(prefix, sep, ci)] = cc
            
        return stage_comms
    
    def get_communities_greedymodularity(self, minEdgeScore = 3.0, resolution=0.5, prefix="Module", sep="_", score_field="score"):
               
        sub_kg_ud = self._filter_edge_score(score_field, minEdgeScore=minEdgeScore)
               
        stage_comms = {}
        prefix="cc"
        for ci, cc in enumerate(nx.community.greedy_modularity_communities(sub_kg_ud, weight=score_field, resolution=resolution)):
            stage_comms["{}{}{}".format(prefix, sep, ci)] = cc
            
        return stage_comms
    
    
    def get_communities_asyn_label_propagation(self, minEdgeScore = 3.0, prefix="Module", sep="_", score_field="score", seed=None):
               
        sub_kg_ud = self._filter_edge_score(score_field, minEdgeScore=minEdgeScore)
        
        if seed is None:
            seed = self.random_state
               
        stage_comms = {}
        prefix="cc"
        for ci, cc in enumerate(nx.community.asyn_lpa_communities(sub_kg_ud, weight=score_field, seed=seed)):
            stage_comms["{}{}{}".format(prefix, sep, ci)] = cc
            
        return stage_comms
    
    
    def get_communities(self, minEdgeScore = 3.0, resolution=5, prefix="Module", sep="_", score_field="score"):
        
        sub_kg_ud = self._filter_edge_score(score_field, minEdgeScore=minEdgeScore)
        partition = community.best_partition(sub_kg_ud, weight=score_field, resolution=resolution, random_state=self.random_state)
        
        rev_partition = defaultdict(set)
        for x in partition:
            rev_partition[ "{}{}{}".format(prefix, sep, partition[x]) ].add(x)
            
        return rev_partition
    
    def get_communities_ecg(self, minEdgeScore = 3.0, resolution=5, ens_size=16, prefix="Module", sep="_", score_field="score"):
        import partition_networkx
        import community
        
        sub_kg_ud = self._filter_edge_score(score_field, minEdgeScore=minEdgeScore)
        partition = community.ecg(sub_kg_ud, weight=score_field, resolution=resolution, ens_size=ens_size)
        
        rev_partition = defaultdict(set)
        for x in partition:
            rev_partition[ "{}{}{}".format(prefix, sep, partition[x]) ].add(x)
            
        return rev_partition
      
    
    
    def get_communities_negpos(self, max_comm_size=50, prefix="Module", sep="_", score_field="score"):
        
        #sub_kg = self.kg.edge_subgraph([x for x in self.kg.edges if ((minEdgeScore is None) or (self.kg.edges[x].get(score_field, 0) > minEdgeScore))]).copy()
        sub_kg = self.kg.copy()
        sub_kg_ud = ig.Graph.from_networkx(sub_kg.to_undirected())
        
        G_pos = sub_kg_ud.subgraph_edges(sub_kg_ud.es.select(lambda x: x[score_field] > 0), delete_vertices=False)
        G_neg = sub_kg_ud.subgraph_edges(sub_kg_ud.es.select(lambda x: x[score_field] <= 0), delete_vertices=False)
        G_neg.es[score_field] = [-w for w in G_neg.es[score_field]]
        
        memberships, diff = la.find_partition_multiplex(graphs=[G_pos, G_neg],
                                                        partition_type=la.ModularityVertexPartition,
                                                        max_comm_size=max_comm_size,
                                                        weights=score_field,
                                                        layer_weights=[1,-1],)
        
        node_mods = [(x["_nx_name"], y) for x,y in zip(G_pos.vs, memberships)]
        
        mod2nodes = defaultdict(set)
        for x in node_mods:
            node = x[0]
            mod = x[1]
            mod = "{}{}{}".format(prefix, sep, mod)
            mod2nodes[mod].add(node)
            
        return mod2nodes
    
    def get_communities_infomap(self, prefix="Module", sep=","):
        """
        Partition network with the Infomap algorithm.
        Annotates nodes with 'community' id and return number of communities found.
        """
        
        infomapWrapper = infomap.Infomap("--two-level --silent --flow-model directed")

        print("Building Infomap network from a NetworkX graph...")
        infomapWrapper.add_networkx_graph(self.kg)

        print("Find communities with Infomap...")
        infomapWrapper.run();

        print("Found %d modules with codelength: %f" % (infomapWrapper.num_leaf_modules,  infomapWrapper.codelength))
        
        communities = defaultdict(set)
        nodeModules = infomapWrapper.get_modules()
        allnodes = [x for x in self.kg.nodes]
        
        for node in nodeModules:
            cluster = nodeModules[node]
            communities["{}{}{}".format(prefix, sep, cluster)].add(allnodes[node])

        return communities
    
    
    def get_communities_markovclustering(self, prefix="cluster", inflation=1.4):

        matrix = scipy.sparse.csr_matrix(nx.to_scipy_sparse_array(self.kg))
        result = mc.run_mcl(matrix, inflation=inflation)
        clusters = mc.get_clusters(result)
        nodelist = [x for x in self.kg.nodes]

        stage_comms = {"{}_{}".format(prefix, xi): [nodelist[y] for y in x] for xi, x in enumerate(clusters)}
        
        return stage_comms
    
    
    def get_communities_link(self, minEdgeScore=3.0, threshold=0.15, score_field="score", prefix="Module", sep="_"):
        
        link_comm_file = os.path.join(os.path.dirname(__file__), "link_clustering.py")
        interpreter = sys.executable
        
        edgefile = tempfile.NamedTemporaryFile(mode="w")
        print(edgefile.name)

        printedEdges = 0
        for edge in self.kg.edges:
            src = edge[0]
            tgt = edge[1]
            
            score = self.kg.edges[edge][score_field]
            
            if not minEdgeScore is None:
                if score >= minEdgeScore:
                    print(src, tgt, score, sep="\t", file=edgefile)
                    printedEdges += 1
            else:
                print(src, tgt, score, sep="\t", file=edgefile)
                printedEdges += 1
                
        print("Printed Edges", printedEdges)           
        edgefile.seek(0)
        
        syscall = "{} {} {} --threshold {}".format(interpreter, link_comm_file, edgefile.name, threshold)
        print(syscall)
        os.system(syscall)
        
        modulefile = glob.glob("{}_*.comm2nodes.txt".format(edgefile.name))[0]
        
        communities = {}        
        for line in open(modulefile):
            line = line.strip().split("\t")
            
            #7	EFO:0006919	EFO:0008344	MONDO:0003008	FHIT	WWOX	CTNNB1
            moduleID = line[0]
            moduleNodes = line[1:]
            communities["{}{}{}".format(prefix, sep, moduleID)] = moduleNodes
               
        edgefile.close()       
        
        outfiles = glob.glob("{}*".format(edgefile.name))
        for ofile in outfiles:
            os.remove(ofile)     
                        
        return communities
    
    def describe_communities(self, comms):
        
        print("Number of communities:", len(comms))
        
        commSizes = [len(comms[x]) for x in comms]
        quants = [0, 0.25, 0.5, 0.75, 1]
        
        if len(commSizes) > 0:
            commMean = np.mean(commSizes)
            commMedian = np.median(commSizes)
            commQuantiles = np.quantile(commSizes, quants)
        else:
            commMean = "-"
            commMedian = "-"
            commQuantiles = "-"
        
        print("Average community size", commMean)
        print("Median community size", commMedian)
        print("Quantile ({}) community size".format(",".join([str(x) for x in quants])), commQuantiles)
       

    def get_kg_subgraph(self, genes):
        
        ret = KGraph(random_state=self.random_state, kgraph_name="{}_sub".format(self.kgraph_name))
        ret.kg = self.get_nx_subgraph(genes)
        
        return ret


    def get_nx_subgraph(self, genes):
        return self.kg.subgraph(genes).copy() # copy avoids edge view problematics ...

    def plot_graph(self, ax=None, figsize=(6,6), title="", pos=None, close=True, font_size=8,
                   edge_score_normalizer=None,
                   node_score_normalizer=None,
                   edge_cmap = plt.cm.Reds,
                   max_node_size = 200,
                   nodetype2color={"gene": "#239756", "geneset": "#3fc37e", "disease": "#5047ee", "drug": "#3026c1", "NA": "#f37855" },
                   nodecolors = {"gene": "#239756", "geneset": "#3fc37e", "disease": "#5047ee", "drug": "#e600e6", "NA": "#f37855" },
                   nodeshapes = {"gene": "o", "geneset": "s", "disease": "^", "drug": "p", "NA": "o" },
                   edge_score_accessor=lambda x: x.get("score", 0),
                   node_score_accessor=lambda x: x.get("score", 0)
                   ):   
                
        if ax is None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
            
        G = self.kg

        #pos = nx.kamada_kawai_layout(G, pos=nx.spring_layout(G, k=0.15, iterations=20))  # For better example looking
        if pos is None:
            pos = nx.spring_layout(G, k=0.3, iterations=50)
            
        
        nodecolor = None
        if not nodetype2color is None:
            nodecolor = {}
            for node in G.nodes:
                nodecolor[node] = nodecolors.get(sorted(G.nodes[node].get("type", ["NA"]))[0], "#f37855")
                
        if node_score_normalizer is None:         
            allNodeScores = [node_score_accessor(G.nodes[node]) for node in G.nodes]
            node_score_normalizer = plt.Normalize(vmin = min(allNodeScores), vmax=max(allNodeScores))
            
        node_size_func = lambda x: 50 + max_node_size*node_score_normalizer(x)
        
                                        
        nodeshape = None
        if not nodeshapes is None:
            nodeshape = []
            
            nodelists = defaultdict(list)
            nodesizes = defaultdict(list)
            
            for i,node in enumerate(G.nodes):
                
                nodeNodeTypes = G.nodes[node].get("type", ["NA"])
                
                if len(set(nodeNodeTypes).intersection(set(nodeshapes.keys()))) > 0:
                    for ntype in nodeNodeTypes:
                        if ntype in nodeshapes:
                            ns = nodeshapes.get(ntype, "o")               
                            nodelists[ns].append(node)
                            nodesizes[ns].append( node_size_func(node_score_accessor(G.nodes[node])) )
                            continue
                else:
                    ntype = "NA"
                    ns = nodeshapes.get(ntype, "o")               
                    nodelists[ns].append(node)
                    nodesizes[ns].append( node_size_func(node_score_accessor(G.nodes[node])) )
            
            for ns in nodelists:
                nx.draw_networkx_nodes(G, pos, nodelist=nodelists[ns], node_size=nodesizes[ns], ax=ax, node_color=[nodecolor[n] for n in nodelists[ns]], node_shape=ns)
            
            
        else:
            nodesizes = []
            for node in G.nodes:
                nodesizes.append( node_size_func(node_score_accessor(G.nodes[node])) )
                        
            nx.draw_networkx_nodes(G, pos, node_size=nodesizes, ax=ax, node_color=[nodecolor[n] for n in G.nodes], node_shape=nodeshape)
            
            
        posnodes = {}
        for x in pos:
            posnodes[x] = list(pos[x])
            posnodes[x][1] += 0.02
            
        nodelabels = {}
        for x in G.nodes:
            
            if "name" in G.nodes[x]:
                nodeName = G.nodes[x].get("name", x)
                
                if nodeName == x:
                    nodelabels[x] = "{}".format(x)
                else:
                    nodelabels[x] = "{}\n({})".format(nodeName, x)
            else:
                nodelabels[x] = x
            
          
                           
        nx.draw_networkx_labels(G, posnodes, labels=nodelabels, font_size=font_size, ax=ax)
        
        edge_scores = [edge_score_accessor(G.edges[e]) for e in G.edges]
        
        if not edge_score_normalizer is None:
            edge_scores = [edge_score_normalizer(x) for x in edge_scores]
            edge_min = 0.0
            edge_max = 1.0
            
        else:
            if len(edge_scores) == 0:
                edge_min = -0.25
            else:
                edge_min = np.floor(np.min(edge_scores))
                
            if len(edge_scores) == 0:
                edge_max = 1
            else:
                edge_max = np.ceil(np.max(edge_scores))        
                
                            
        nx.draw_networkx_edges(G, pos, width=2, edge_vmin=edge_min, edge_vmax=edge_max, edge_cmap = edge_cmap, edge_color=edge_scores, ax=ax)
        ax.set_title(title, loc='left')
        
        if close:
            
            if edge_score_normalizer is None:
                edge_score_normalizer = plt.Normalize(vmin = edge_min, vmax=edge_max)
                
            sm = plt.cm.ScalarMappable(cmap=edge_cmap, norm=edge_score_normalizer)
            cbar = plt.gcf().colorbar(sm, orientation="horizontal", pad=0.1, ax=ax, shrink=1.0)
            cbar.ax.set_xlabel("Edge Score".format())
            
            l = matplotlib.ticker.MaxNLocator(nbins=5)
            l.create_dummy_axis()
            newticks = l.tick_values(edge_min, edge_max)
            cbar.ax.set_xticks(newticks)
            
            
            def get_legend(norm, ax, num=5):
                
                sizes = np.linspace(norm.vmin, norm.vmax, num)
                
                plotsizes = [node_size_func(x) for x in sizes]

                xsizes = [i for i in range(0, len(sizes))]
                ysizes = [0]*len(sizes)
                sc = ax.scatter(xsizes, ysizes, s=plotsizes, color='#239756')
                
                for i, (xi, yi) in enumerate(zip(xsizes, ysizes)):
                    plt.text(xi, yi-0.05, "{:.3f}".format(sizes[i]), va='bottom', ha='center', fontsize="small")

                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['left'].set_visible(False)
            
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(cbar.ax)
            sax = divider.append_axes('right', size='100%', pad=0.2)
            get_legend(node_score_normalizer, sax)

            plt.show()
            plt.close()
            
        return pos
                
    def get_node_degrees(self, intype, connecttype=None):
        
        relNodes = self.get_nodes(intype)
        
        if connecttype is None:
            all_degrees = [x[1] for x in nx.degree(self.kg, relNodes)]        
        else:
            all_degrees = []
            for node in relNodes:
                num_edges = len(self.get_edges_to_type(node, connecttype))
                all_degrees.append(num_edges)
                
        return all_degrees
                
        
                
class NetworkScorer:
    """_summary_ Interface class for any network scoring class.
    
    Hint: All scores must have their minimum at 0 ! This is required such that differential KGs can be successfully generated!    
    
    """
    
    def __init__(self, random_state=42) -> None:
        
        self.random_state=random_state
        random.seed(random_state)
        
        self.logger = logging.getLogger("NetworkScorer")
        self.logger.setLevel(logging.INFO)

    def score(self):
        pass

    
class MeanNetworkScorer(NetworkScorer):
    
    def __init__(self, random_state=42) -> None:
        super().__init__(random_state)
        
        
        def scoring_represses(x, y):
            return max(x - y, 0)

        def scoring_activates(x, y):
            return x*y

        def scoring_interacts(x,y):
            return x*y

        def scoring_null(x,y):
            return 0

        self.scoring_gene_gene_expression = {"represses": scoring_represses,
                                             "activates": scoring_activates,
                                             "interacts": scoring_interacts,
                                             "-": scoring_null}

        self.scoring_interactions = {"relevant_in": scoring_interacts,
                                "activates": scoring_activates,
                                "part_of": scoring_interacts,
                                "represses": scoring_represses,
                                "interacts": scoring_interacts,
                                "target_of": scoring_interacts,
                                "affected_by": scoring_interacts,
                                "-": scoring_null
                                }





    def score_edges(self, kgraph, value_accessor, edge_accessor, scorers, src_types=None, tgt_types=None, ignore_edge_types=None):
        
        for edge in kgraph.get_edges_between_ntypes(src_types=src_types, tgt_types=tgt_types):
            
            srcTypes, tgtTypes = kgraph.get_edge_node_types(edge)           
            etypes = set([x for x in product(srcTypes, tgtTypes)])
            
            if not ignore_edge_types is None:               
                if len(etypes.intersection(ignore_edge_types)) > 0:
                    #ignore
                    continue
    
            inExpr = value_accessor(kgraph.get_node_data(edge[0]))
            outExpr = value_accessor(kgraph.get_node_data(edge[1]))
            
            edgeType = kgraph.kg.edges[edge].get( edge_accessor , "-")
            
            edgeScore = scorers[edgeType](inExpr, outExpr)
            
            kgraph.kg.edges[edge]["score"] = edgeScore


    def score_nodes_from_properties(self, kgraph, ntypes=["measured_expression"],
                                    scoring_function=lambda x: x["mean"]*x["perc_expr"],
                                    property_accessor=lambda x: x.get("expression", None)):
        
        
        
        for node in kgraph.kg.nodes:
            
            typeIntersect = set(ntypes).intersection( kgraph.get_node_data(node).get("type", set()) )
            if len(typeIntersect) > 0:
                nodeData = property_accessor(kgraph.kg.nodes[node])
                
                if nodeData is None:
                    continue
                
                kgraph.kg.nodes[node]["score"] = scoring_function(nodeData)
        
    
            
    def score_nodes_from_edges(self, kgraph, ntype="geneset", consider_edge_type=[("gene", "geneset")], scoring_function=None, overwrite_score=False):
        
        assert(not scoring_function is None)
        
        src_types = [x[0] for x in consider_edge_type]
        tgt_types = [x[1] for x in consider_edge_type]
        
        for node in kgraph.kg.nodes:
            
            if kgraph.get_node_data(node).get("type", None) == ntype:
                
                edges = kgraph.get_node_edges(node, src_types=src_types, tgt_types=tgt_types)
                
                nodeScore = scoring_function(kgraph.kg, node, edges)
                
                if overwrite_score or (not "score" in kgraph.kg.nodes[node]):
                    kgraph.kg.nodes[node]["score"] = nodeScore
    

    def score_nx(self, kg:nx.DiGraph):
        
        okg = KGraph(kgraph_name="nxsub_kg")
        okg.kg = kg
        
        self.score(okg)
        
        return okg.kg
        
        

    def score(self, kgraph:KGraph):
        
        
        def gene_geneset(kg, node, edges):
            
            allExpr = list()
            for edge in edges:
                
                otherEnd = edge[0] if node==edge[1] else edge[1]
                otherExpr = kg.nodes[otherEnd].get("score", 0)
                
                allExpr.append(otherExpr)
            
            return np.mean(allExpr)
                
                

        def geneset_geneset(kg, node, edges):
            
            allExpr = list()
            for edge in edges:
                
                otherEnd = edge[0] if node==edge[1] else edge[1]
                otherExpr = kg.nodes[otherEnd].get("score", 0)
                
                allExpr.append(otherExpr)
            
            return np.mean(allExpr)

        def get_score(x):
            if "gene" in x.get("type", ["-"]):
                return x.get("score", 0)
            else:
                return x.get("score", 0)
        
        
        self.score_edges(kgraph, get_score, "type", self.scoring_gene_gene_expression, src_types=["gene"], tgt_types=["gene"])
    
        _ = kgraph.score_nodes_hierarchically(ntype="geneset", target_ntype="gene")
        _ = kgraph.score_nodes_hierarchically(ntype="disease", target_ntype="gene")
        _ = kgraph.score_nodes_hierarchically(ntype="ncRNA", target_ntype="gene")
        _ = kgraph.score_nodes_hierarchically(ntype="drug", target_ntype=["gene", "disease", "geneset", "ncRNA"])

        #
        ## are these calls actually needed?
        #

        # drug interactions
        self.score_nodes_from_edges(kgraph, ntype="drug", consider_edge_type=[("gene", "drug")], scoring_function=gene_geneset)
        self.score_nodes_from_edges(kgraph, ntype="drug", consider_edge_type=[("disease", "drug")], scoring_function=geneset_geneset)

        # genesets
        self.score_nodes_from_edges(kgraph, ntype="geneset", consider_edge_type=[("geneset", "geneset")], scoring_function=geneset_geneset)

        # non-coding RNA scoring
        self.score_nodes_from_edges(kgraph, ntype="ncRNA", consider_edge_type=[('ncRNA', 'gene'), ('gene', 'ncRNA')], scoring_function=geneset_geneset)
        self.score_nodes_from_edges(kgraph, ntype="ncRNA", consider_edge_type=[('ncRNA', 'ncRNA')], scoring_function=geneset_geneset)

        #
        ## all nodes should have a score by now!
        #

        self.score_edges(kgraph, get_score, "type", self.scoring_interactions, src_types=None, tgt_types=None, ignore_edge_types=[("gene", "gene")])
        self.calculate_edge_zscores(kgraph, score_accessor="score", edge_accessor="type", src_types=None, tgt_types=None, ignore_edge_types=None)

    def calculate_edge_zscores(self, kgraph, score_accessor, edge_accessor, src_types=None, tgt_types=None, ignore_edge_types=None):

        etype2scores = defaultdict(list)

        for edge in kgraph.get_edges_between_ntypes(src_types=src_types, tgt_types=tgt_types):
            
            srcTypes, tgtTypes = kgraph.get_edge_node_types(edge)
            edgeType = kgraph.kg.edges[edge].get( edge_accessor , "-")
            
            for srcType in srcTypes:
                for tgtType in tgtTypes:
                    
                    if not ignore_edge_types is None:
                        if (srcType, tgtType) in ignore_edge_types:
                            continue
                    
                    etype = (srcType, tgtType, edgeType)

                    edge_score = kgraph.kg.edges[edge][ score_accessor ]
                    etype2scores[etype].append(edge_score)

        etype2std = {x: np.std(etype2scores[x]) for x in etype2scores}
        etype2mean = {x: np.mean(etype2scores[x]) for x in etype2scores}

        #print("Identified the following edge types:", *[x for x in etype2scores])

        for edge in kgraph.get_edges_between_ntypes(src_types=src_types, tgt_types=tgt_types):
            
            srcTypes, tgtTypes = kgraph.get_edge_node_types(edge)
            
            if not ignore_edge_types is None:
                if (srcType, tgtType) in ignore_edge_types:
                    continue

            edgeType = kgraph.kg.edges[edge].get( edge_accessor , "-")
            
            for srcType in srcTypes:
                for tgtType in tgtTypes:
                    etype = (srcType, tgtType, edgeType)
            
                    edge_score = kgraph.kg.edges[edge][score_accessor]
                    edge_zscore = (edge_score-etype2mean[etype])/etype2std[etype]
                    kgraph.kg.edges[edge]["{}_zscore".format(score_accessor)] = edge_zscore




class NetworkExtender:
    
    def __init__(self) -> None:
        pass
    
    
    def extend_network(self, nodes, fullKG: KGraph, radius=1, scorer:NetworkScorer=None,
                       min_children_gs=2,
                       max_size_gs=100,
                       minFraction_small = 0.4,
                       minFraction_large = 0.5,
                       node_types = ["geneset", "disease", "ncRNA"],
                       minGeneSpec={"geneset": 0.8, "disease": 0.6},
                       min_edge_score=1.0,
                       score_field="score",
                       verbose=False):
        
        orig_kg = fullKG.kg
        
        
        
        relNodes = set()
        for en in nodes:
            sg = nx.ego_graph(orig_kg, en, radius=radius)
            
            # take all non-gene nodes, but no genesets/diseases
            geneNodes = [x for x in sg.nodes if fullKG.node_type_overlap(x, ["gene", "TF"])]
            
            acceptNodes = []
            #acceptNodes = [x for x in sg.nodes if not x in geneNodes]
            #acceptNodes = [x for x in acceptNodes if not fullKG.node_type_overlap(x, ["geneset", "disease", "ncRNA", "drug"])]
            
            for n in sg.nodes:
                nodeTypes = orig_kg.nodes[n].get("type", "")
                if len(set(node_types).intersection(nodeTypes)) > 0:
                                
                    geneNeighborsS = [x for x in orig_kg.successors(n) if "gene" in orig_kg.nodes[x].get("type", [])]
                    geneNeighborsP = [x for x in orig_kg.predecessors(n) if "gene" in orig_kg.nodes[x].get("type", [])]
                    
                    # TODO add spec-measure for genesets or diseases!
                    accGeneNeighborsS = []
                    for x in geneNeighborsS:
                        acceptX = False
                        for nodeType in nodeTypes:
                            acceptX = acceptX or (orig_kg.edges[(n, x)].get("{}_spec".format(nodeType), 1.0) > minGeneSpec.get(nodeType, 0))
                            if acceptX:
                                break
                        if acceptX:
                            accGeneNeighborsS.append(x)
                        
                    accGeneNeighborsP = []
                    for x in geneNeighborsP:
                        acceptX = False
                        for nodeType in nodeTypes:
                            acceptX = acceptX or (orig_kg.edges[(x,n)].get("{}_spec".format(nodeType), 1.0) > minGeneSpec.get(nodeType, 0))
                            if acceptX:
                                break
                        if acceptX:
                            accGeneNeighborsP.append(x)
                                        
                    geneNeighbors=list(set(accGeneNeighborsS+accGeneNeighborsP))
                    
                    #if len(geneNeighbors) == 0:
                    #    print("zero", n, geneNeighbors)
                    
                    if len(geneNeighbors) > max_size_gs:
                        continue

                    accEdges = []
                    #get all edges of n in origkg
                    for u, v, data in [x for x in orig_kg.in_edges(n, data=True)] + [x for x in orig_kg.out_edges(n, data=True)]:
                        #if edge could also be in extended nx
                        if u in nodes or v in nodes:                           
                            #if not edge score > minscore: continue
                            if data.get(score_field, 0) > min_edge_score:
                                accEdges.append((u,v))

                    if len(accEdges) == 0:
                        continue
                                        
                    if len(geneNeighbors) == 0:
                        containedNeighbours = 0
                        fractionNeighbours = 0
                    else:
                        containedNeighbours = len(set(geneNeighbors).intersection( nodes ))
                        fractionNeighbours = containedNeighbours/len(geneNeighbors)
                    
                    if len(geneNeighbors) < min_children_gs:
                        continue
                    
                    if len(geneNeighbors) >= 5:
                        
                        if fractionNeighbours > minFraction_large and containedNeighbours >= min_children_gs:
                            if verbose:
                                print(n, "large", "\"{}\"".format(fullKG.kg.nodes[n].get("name", n)), len(geneNeighbors), containedNeighbours, fractionNeighbours)
                            acceptNodes.append(n)
                            
                    else:
                        if fractionNeighbours > minFraction_small:
                            if verbose:
                                print(n, "small", "\"{}\"".format(fullKG.kg.nodes[n].get("name", n)), len(geneNeighbors), containedNeighbours, fractionNeighbours)
                            acceptNodes.append(n)
            
            relNodes.update(acceptNodes)
            
        print("Input Graph Nodes", len(nodes))
        print("Extended Graph Nodes", len(relNodes))
        
        enodes = set(list(relNodes) + list(nodes))
                
        okg = KGraph(fullKG.random_state, kgraph_name="nwe_sub")
        okg.kg = nx.subgraph(fullKG.kg, enodes).copy()
        
        print("Extended Graph Edges", len(okg.kg.edges))
        
        if not scorer is None:
            scorer.score(okg)
            
        return okg
    
    def extend_network_force(self, eKG:KGraph, fullKG:KGraph, nodetype, acceptor=None, edge_acceptor=None, minSpec=0.5):
        
        graphnodes = [x for x in eKG.kg.nodes]
        for n in graphnodes:
            
            possibleEdges = fullKG.get_edges_to_type(n, nodetype)
            
            for oNode in possibleEdges:
                
                if not acceptor is None:
                    if acceptor(oNode, fullKG) == False:
                        continue
                                        
                if (n, oNode) in fullKG.kg.edges:
                    edge = (n, oNode)
                else:
                    edge = (oNode, n)
                    
                    
                if not edge_acceptor is None:
                    if edge_acceptor(edge, fullKG) == False:
                        continue
                    
                ntypes = fullKG.get_node_types(n)
                
                ignore=False
                for ntype in ntypes:
                    specName = "{}_spec".format(ntype)
                    if specName in fullKG.kg.nodes[n]:
                        if fullKG.kg.nodes[n][specName] < minSpec:
                            ignore = True

                if ignore:
                    continue
                
                print(n, "-->", "({})".format(fullKG.kg.edges[edge]), oNode, fullKG.kg.nodes[oNode])

                if not oNode in eKG.kg.nodes:
                    eKG.kg.add_node(oNode)
                    
                    for x in fullKG.kg.nodes[oNode]:
                        eKG.kg.nodes[oNode][x] = fullKG.kg.nodes[oNode][x]
                    
                eKG.kg.add_edge( edge[0], edge[1] )
                    
                for x in fullKG.kg.edges[(edge[0], edge[1])]:
                    eKG.kg.edges[(edge[0], edge[1])][x] = fullKG.kg.edges[(edge[0], edge[1])][x]
            
    
    def extend_nodetypes(self, eKG:KGraph, fullKG:KGraph, nodetype, node_score_accessor=None, edge_score_accessor=None, min_node_score=0.5, verbose=False):
        
        exEdgeScores = [edge_score_accessor(eKG.kg.edges[e]) for e in eKG.kg.edges]
        edgeScoreQuantiles = np.quantile(exEdgeScores, [0.25, 0.75])

        graphnodes = [x for x in eKG.kg.nodes]
        for n in graphnodes:
            
            possibleEdges = fullKG.get_edges_to_type(n, nodetype)

            if len(possibleEdges) == 0:
                continue

            if node_score_accessor(eKG.kg.nodes[n]) < min_node_score:
                if verbose:
                    print(n, "node_score", node_score_accessor(eKG.kg.nodes[n]), "<", min_node_score)
                continue
            
            #print(n, possibleEdges)
            #print(n)
            #print(eKG.kg.nodes[n])

            for oNode in possibleEdges:
                                            
                if (n, oNode) in fullKG.kg.edges:
                    edge = (n, oNode)
                else:
                    edge = (oNode, n)

                edge_score = edge_score_accessor(fullKG.kg.edges[edge])
                if edgeScoreQuantiles[0] <= edge_score <= edgeScoreQuantiles[1]:
               
                    if verbose:
                        print(n, "-->", "({})".format(fullKG.kg.edges[edge]), oNode, fullKG.kg.nodes[oNode])

                    if not oNode in eKG.kg.nodes:
                        eKG.kg.add_node(oNode)
                        
                        for x in fullKG.kg.nodes[oNode]:
                            eKG.kg.nodes[oNode][x] = fullKG.kg.nodes[oNode][x]
                        
                    eKG.kg.add_edge( edge[0], edge[1] )
                        
                    for x in fullKG.kg.edges[(edge[0], edge[1])]:
                        eKG.kg.edges[(edge[0], edge[1])][x] = fullKG.kg.edges[(edge[0], edge[1])][x]
                        
                else:
                    if verbose:
                        print(edge, "score not", edgeScoreQuantiles[0], "<", edge_score, "<", edgeScoreQuantiles[1])
       
       

class CommunityTool:
    
    def __init__(self) -> None:
        pass
    
    def scatter_community_scores(self, module_detail):
        
        for xi, x in enumerate(module_detail):
            diffScores = module_detail[x]
            if xi == 0:
                baseline = np.array(diffScores["own_scores"])

            plt.scatter(baseline, diffScores["other_scores"], label=x, s=0.1)
            
            m, b = np.polyfit(baseline, diffScores["other_scores"], 1)

            xs = np.linspace(min(baseline), max(baseline), 1000)
            
            #use red as color for regression line
            plt.plot(xs, m*xs+b)

        plt.xlabel("Own Scores")
        plt.ylabel("Other Scores")

        plt.yscale('log')
        plt.xscale('log')

        plt.legend()
    
    
    def visualize_communities(self, details, title, subsetOrderFunc=None, field="median"):
        
        
        allSubsets = set()

        for mod in details:
            for subset in details[mod]:
                allSubsets.add(subset)

        if not subsetOrderFunc is None:
            allSubsets = sorted([x for x in allSubsets], key=lambda x: subsetOrderFunc(x))

        modScores = []
        modNames = []
        for mod in details:

            scores = []
            for subset in allSubsets:
                score = details[mod].get(subset, {}).get("{}_other".format(field), 0)
                scores.append(score)

            modScores.append( tuple(scores) ) #medians
            modNames.append(mod)
            
        modDF = pd.DataFrame(np.array(modScores).T)
        modDF.columns = modNames
        modDF.index = ["{} ({})".format(x, field) for x in allSubsets] #["{}_mean".format(x) for x in allSubsets]

        g = sns.clustermap(modDF, figsize=(12, 4), row_cluster=False, xticklabels=True, yticklabels=True)    
        g.ax_heatmap.set_title("{} modules".format(title))
            
         
    def sort_communities(self, comm_details, details=False):
        
        avgCohens = {}
        for c in comm_details:
            
            cohens = []
            for reg in comm_details[c]:
                cohens.append( comm_details[c][reg].get("cohend", 0) )
            
            avgCohen = sum(cohens)    
            avgCohen = avgCohen / (len(cohens))
            
            avgCohens[c] = avgCohen
         
        sortedComms = sorted([x for x in comm_details], key=lambda x: abs(avgCohens[x]), reverse=True)
        
        if details:
            return {x: avgCohens[x] for x in sortedComms}
        
        return sortedComms
       
    def plot_community(self, KGs, communityNodes, own, main_net=None, num_columns=4,
                         title=None, nodecolors = {"gene": "#239756", "geneset": "#3fc37e", "disease": "#5047ee", "drug": "#3026c1", "NA": "#f37855" },
                         outfile=None, dpi=500,
                         edge_score_accessor=lambda x: x.get("score", 0),
                         node_score_accessor=lambda x: x.get("score", 0),
                         verbose=False):

        assert( own in KGs)
            
        import math
        
        numCols = min(len(KGs)+1, num_columns)
        numRows = math.ceil((len(KGs))/numCols)
        
        fig, axs = plt.subplots(numRows, numCols, figsize=(6*numCols, 6*numRows+1), constrained_layout=True)
        flataxs = axs.flat
            
        # which genes to use for scoring
        scoreGenes = communityNodes
        if not main_net is None:
            scoreGenes = main_net
            
        ownKG = KGs[own].get_kg_subgraph(communityNodes)           
        ownEdgeScore = ownKG.get_edge_scores(nodes=scoreGenes, score_accessor=edge_score_accessor)
                    


        edgeScores = {}      
        allEdgeScores = []
        allNodeScores = []    
        for kgname in KGs:
            
            plotKG = KGs[kgname].get_kg_subgraph(communityNodes)
            kgScores = plotKG.get_edge_scores(nodes=scoreGenes, score_accessor=edge_score_accessor)
            edgeScores[kgname] = kgScores
            allEdgeScores = allEdgeScores + list(kgScores)
            allNodeScores = [node_score_accessor(plotKG.kg.nodes[node]) for node in plotKG.kg.nodes] + allNodeScores
            
        edgeScore_max = 0
        edgeScore_min = 0
        
        if len(allEdgeScores) > 0:
            edgeScore_max = np.ceil(np.max(allEdgeScores))
            edgeScore_min = np.floor(np.min(allEdgeScores))
            
        if verbose:
            print("Min Edge Score: ", min(allEdgeScores), edgeScore_min)
            print("Max Edge Score: ", max(allEdgeScores), edgeScore_max)
            
           
        
        
        nodeScore_max = 0
        nodeScore_min = 0   
        if len(allNodeScores) > 0:
            nodeScore_max = np.ceil(np.max(allNodeScores))
            nodeScore_min = np.floor(np.min(allNodeScores))
          
        node_score_normalizer = plt.Normalize(vmin = nodeScore_min, vmax=nodeScore_max)
            
        if verbose:
            print("Min Node Score: ", min(allNodeScores), nodeScore_min)
            print("Max Node Score: ", max(allNodeScores), nodeScore_max) 
        
        if edgeScore_min < 0 and edgeScore_max > 0:
            edge_cmap = plt.cm.coolwarm
            #edge_score_normalizer = MidpointNormalize(vmin = edgeScore_min, midpoint=0, vmax=edgeScore_max)
            edge_score_normalizer = matplotlib.colors.TwoSlopeNorm(vmin=edgeScore_min, vcenter=0, vmax=edgeScore_max)
            
        else:
            
            if edgeScore_min < 0 and edgeScore_max < 0:
                edge_cmap = matplotlib.cm.get_cmap('Blues_r')
            else:
                edge_cmap = plt.cm.Reds
            
            edge_score_normalizer = plt.Normalize(vmin = edgeScore_min, vmax=edgeScore_max)
            
        sm = plt.cm.ScalarMappable(cmap=edge_cmap, norm=edge_score_normalizer)
        
        # plot once in order to get node positions!
        pos=ownKG.plot_graph(ax=flataxs[0], title="--", close=False,
                             nodetype2color=nodecolors, edge_score_accessor=edge_score_accessor)
        flataxs[0].clear()
        max_node_size=200
        
        if not KGs is None:
            for kgi, kgname in enumerate(KGs):

                plotKG = KGs[kgname].get_kg_subgraph(communityNodes)
                
                kgScores = edgeScores[kgname]
                kgMedian = np.median(kgScores)
                kgMean = np.mean(kgScores)
                
                
                _=plotKG.plot_graph( ax=flataxs[kgi], pos=pos, title="{} (median score: {:.3f}, mean {:.3f})".format(kgname, kgMedian, kgMean),
                                    close=False, nodetype2color=nodecolors,
                                    edge_score_normalizer=edge_score_normalizer,
                                    node_score_normalizer=node_score_normalizer,
                                    edge_cmap=edge_cmap,
                                    edge_score_accessor=edge_score_accessor,
                                    node_score_accessor=node_score_accessor,
                                    max_node_size=max_node_size)
        
        
        cbar = fig.colorbar(sm, orientation="horizontal", ax=axs, pad=0.1, shrink=0.3)
        cbar.ax.set_xlabel("Edge Score".format())
        
        if edgeScore_min == 0:
            edgeScore_min = -1
        
        #cbarxticks = list(cbar.ax.get_xticks())
        #print("CbarTicks Orig", cbarxticks)
        #if not edgeScore_min in cbarxticks:
        #    cbarxticks.append(edgeScore_min)
        #if not edgeScore_max in cbarxticks:
        #    cbarxticks.append(edgeScore_max)
        #cbarxticks = sorted(cbarxticks)
        #print("CbarTicks", cbarxticks)
        #cbar.ax.set_xticks(cbarxticks)

        node_size_func = lambda x: 50 + max_node_size*node_score_normalizer(x)


        def get_legend(norm, ax, num=5):
            
            sizes = np.linspace(norm.vmin, norm.vmax, num)
            
            plotsizes = [node_size_func(x) for x in sizes]

            xsizes = [i for i in range(0, len(sizes))]
            ysizes = [0]*len(sizes)
            sc = ax.scatter(xsizes, ysizes, s=plotsizes, color='#239756')
            
            for i, (xi, yi) in enumerate(zip(xsizes, ysizes)):
                plt.text(xi, yi-0.05, "{:.3f}".format(sizes[i]), va='bottom', ha='center', fontsize="small")

            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(cbar.ax)
        sax = divider.append_axes('right', size='100%', pad=0.2)
        get_legend(node_score_normalizer, sax)


        if not title is None:           
            #fig.subplots_adjust(top=0.90)

            plt.gcf().suptitle(title, x=0.025, fontweight="bold")
                    
        if not outfile is None:
            plt.savefig(outfile + ".png", bbox_inches='tight', dpi=dpi)
            plt.savefig(outfile + ".pdf", bbox_inches='tight', dpi=dpi)
                        
                
            with open(outfile + ".tsv", "w") as fout:
                
                print("kgname", "node1", "node1_name", "node1_type", "node1_score", "node2", "node2_name", "node2_type", "node2_score", "edge_score", "edge_type", "node1_dict", "node2_dict", "edge_dict", sep="\t", file=fout)
                
                for kgname in KGs:
                    p_kg = KGs[kgname].get_kg_subgraph(communityNodes)
                    pkg = p_kg.kg
                    for edge in pkg.edges:
                        
                        n1types = pkg.nodes[edge[0]].get("type", [edge[0]])
                        n1score = pkg.nodes[edge[0]].get("score", 0)
                        
                        n2types = pkg.nodes[edge[1]].get("type", [edge[1]])
                        n2score = pkg.nodes[edge[1]].get("score", 0)
                        
                        print(kgname,
                              edge[0], pkg.nodes[edge[0]].get("name", edge[0]), ";".join(n1types), n1score,
                              edge[1], pkg.nodes[edge[1]].get("name", edge[1]), ";".join(n2types), n2score,
                              pkg.edges[edge].get("score", 0), pkg.edges[edge].get("type", 0),
                              pkg.nodes[edge[0]], pkg.nodes[edge[1]], pkg.edges[edge],
                              sep="\t", file=fout)
                                        
        plt.show()
        plt.close()
        
    def compare_modules(self, comms, figsize=(16, 12)):
        
        modules = [x for x in comms]

        simMatrix = np.zeros( (len(modules), len(modules)) )

        for mi1,  mod1 in enumerate(modules):
            for mi2, mod2 in enumerate(modules):

                overlap = len(set(comms[mod1]).intersection(comms[mod2]))
                overlap = overlap / len(comms[mod1])

                simMatrix[ mi1, mi2 ] = overlap


        fig, ax = plt.subplots(figsize=figsize)

        for row in range(0, simMatrix.shape[0]):
            x = [i for i in range(0, simMatrix.shape[1])]
            y = [row] * simMatrix.shape[1]
            size = [simMatrix[row, i]*50 for i in x]
            
            plt.scatter( x, y,  s=size, c="blue")

        
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(1))

        plt.xlim(-1, len(modules))
        plt.ylim(-1, len(modules))

        ax.set_xticklabels(["", ""] + modules + [""], rotation=45, ha='right')
        ax.set_yticklabels(["", ""] + modules + [""], rotation=0, ha='right')

        #print(ax.get_xticklabels())
        
        plt.show()
        
        
    
class DifferentialModuleIdentifier:
    
    def __init__(self) -> None:
        pass
    
    

            
    # function to calculate Cohen's d for independent samples
    # from: https://machinelearningmastery.com/effect-size-measures-in-python/
    @classmethod
    def cohend(cls, d1, d2):
        # calculate the size of samples
        n1, n2 = len(d1), len(d2)
        # calculate the variance of the samples
        s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
        # calculate the pooled standard deviation
        s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
        # calculate the means of the samples
        u1, u2 = np.mean(d1), np.mean(d2)
        # calculate the effect size
        return (u1 - u2) / s
            
    
    def identify_differential_communities(self, communities, ref_kg, KGs, sort_function=None, score_field="score", min_nodes=10, min_enriched=0.5, min_effect_size=0.2, all_verbose=False, verbose=False):
        
        
        if sort_function is None:
            
            def ct_sort_function(rel_details):
                ct = CommunityTool()
                return ct.sort_communities(rel_details, details=True)
            
            sort_function = ct_sort_function

        relcomms = []
        reldetails = {}
        
        if type(ref_kg) == str:
            ref_kg = [ref_kg]
        
        for cID in communities:
            
            if len(communities[cID]) < min_nodes:
                continue
            
            
            otherKGs = {x: KGs[x] for x in KGs if not x in ref_kg}
            refKGs = {x: KGs[x] for x in KGs if x in ref_kg}
            
            commScores = self.score_subgraphs_for_subnet(otherKGs, communities[cID], score_field=score_field)           
            ownScoresDict = self.score_subgraphs_for_subnet(refKGs, communities[cID], score_field=score_field)
            
            # list of edge scores!
            ownScores = []
            for x in ownScoresDict:
                ownScores += ownScoresDict[x]
            
            diffScores = defaultdict(dict)
            for x in commScores:
                diffScores[x] = {}
                diffScores[x]["median_other"] = np.median(commScores[x])
                diffScores[x]["mean_other"] = np.mean(commScores[x])
                
                diffScores[x]["median_own"] = np.median(ownScores)
                
                diffScores[x]["cohend"] = self.cohend(commScores[x], ownScores)
                diffScores[x]["subset_scores"] = commScores[x]
                diffScores[x]["base_scores"] = ownScores
                                
                #if np.median(commScores[x]) == 0:
                #    logFC = 0
                #else:
                #    logFC = np.log2(np.median(commScores[x]) / np.median(ownScores))                               
                #diffScores[x]["logFC"] = logFC
                
                diffScores[x]["ks"] = scipy.stats.ks_2samp(commScores[x], ownScores)
                                                              
            enrichedModule = 0
            for x in diffScores:
                if abs(diffScores[x]["cohend"]) >= min_effect_size:
                    enrichedModule+=1

            accepted=False
            if (enrichedModule / len(diffScores)) >= min_enriched:
                accepted=True                    
                relcomms.append(cID)
                reldetails[cID] = diffScores
                
            
            if all_verbose or verbose:
                if all_verbose or accepted:
                    print("community", cID, len(communities[cID]), accepted)
                    
                    # add reference info for accepted
                    for refName in refKGs:
                        diffScores[refName]["median_other"] = np.median(ownScoresDict[refName])
                        diffScores[refName]["mean_other"] = np.mean(ownScoresDict[refName])
                    
                    for x in diffScores:
                        print(x, diffScores[x])
                
                
        relcomms = sort_function(reldetails) 
        reldetails = {x: reldetails[x] for x in relcomms}                
        
        return relcomms, reldetails

    def score_subgraphs_for_subnet(self, KGs, subnet, score_field="score"):
        
        kgSubgraphs = {}        
        
        for group in KGs:
            kgSubgraphs[group] = KGs[group].kg.subgraph(subnet)
        
        kgScores = self.score_subgraphs(kgSubgraphs, score_field=score_field)
        
        return kgScores
    
    def score_subgraphs(self, Gs, score_field="score", default_score=0):
        
        scores = defaultdict(list)
        for gname in Gs:
            
            G = Gs[gname]
            
            for edge in G.edges:
                
                escore = G.edges[edge].get(score_field, default_score)
                scores[gname].append(escore)
                
        return scores            

class GenesetAnnotator:
    
    def __init__(self) -> None:
        pass
    
    
    
    def annotate_genesets(self, kg:KGraph, settype="disease", targettype="gene"):
        
        targetNodes = kg.get_nodes(targettype)
        
        nodeTypeCounter = kg.get_node_types(single=True)
        
        assert settype in nodeTypeCounter
        
        setTypeCount = nodeTypeCounter[settype]
        
        targetNodesWithEdges = 0
        attr = "{}_spec".format(settype)
        attr_dist = []
        
        for node in targetNodes:
            
            targetSetNodes = kg.get_edges_to_type(node, settype)
            
            if len(targetSetNodes) == 0:
                continue
            
            # score from: Disease Specificity Index; https://www.disgenet.org/dbinfo
            nodeSpec = np.log2(len(targetSetNodes) / setTypeCount) / np.log2(1.0/setTypeCount)
            
            targetNodesWithEdges = targetNodesWithEdges + 1
            
            kg.kg.nodes[node][ attr ] = nodeSpec
            attr_dist.append(nodeSpec)

        attrMean = np.mean(attr_dist)
        attrSD = np.std(attr_dist)

        attr_z = "{}_zscore".format(attr)

        for node in targetNodes:
            if attr in kg.kg.nodes[node]:
                kg.kg.nodes[node][attr_z] = (kg.kg.nodes[node][attr]-attrMean)/attrSD      
            
        print("Processed {} of {} target nodes for settype={} and targettype={}".format(targetNodesWithEdges, len(targetNodes), settype, targettype))
            
            
class DifferentialKG:
    
    def __init__(self) -> None:
        pass
    
    
    def _get_edge_fold_changes(self, kg1:KGraph, kg2: KGraph):

        edgeDiffs = {}
        
        for edge in kg1.kg.edges:

            edge1Score = kg1.kg.edges[edge].get("score", 0)
            edge2Score = kg2.kg.edges[edge].get("score", 0)

            #if edge2Score == 0 and edge1Score != 0:
            #    scoreDiff = 512/1
            #elif edge1Score == 0 and edge2Score != 0:
            #    scoreDiff = 1/512
            #elif edge1Score == 0 and edge2Score == 0:
            #    scoreDiff = 1
            #else:
            #    try:
            #        scoreDiff = edge1Score / edge2Score
            #    except:
            #        print("err", edge, edge1Score, edge2Score)
            #        raise ValueError()
            #scoreDiff = np.log2(scoreDiff)
            scoreDiff=self._get_log_foldchange(edge1Score, edge2Score, log2=True)

            edgeDiffs[edge] = scoreDiff

        return edgeDiffs

    def _get_log_foldchange(self, node1Score, node2Score, log2=True):
        if node2Score == 0 and node1Score != 0:
            scoreDiff = 512/1
        elif node1Score == 0 and node2Score != 0:
            scoreDiff = 1/512
        elif node1Score == 0 and node2Score == 0:
            scoreDiff = 1
        else:
            scoreDiff = node1Score / node2Score

        if log2:
            scoreDiff = np.log2(scoreDiff)

        return scoreDiff

    def _get_node_fold_changes(self, kg1:KGraph, kg2: KGraph):

        nodeDiffs = {}
        numDefault = 100
        sdDefault = 1
        
        for node in kg1.kg.nodes:

            node1Num = numDefault
            node1SD = sdDefault
            node1Score = kg1.kg.nodes[node].get("score", 0)

            node2Num = numDefault
            node2SD = sdDefault
            node2Score = kg2.kg.nodes[node].get("score", 0)

            scoreDiff = self._get_log_foldchange(node1Score, node2Score)


            if "gene" in kg1.kg.nodes[node]["type"] and "gene" in kg2.kg.nodes[node]["type"]:
                node1num = kg1.kg.nodes[node].get("expression", {}).get("num_measured", numDefault)
                node1SD = kg1.kg.nodes[node].get("expression", {}).get("sd", sdDefault)

                node2num = kg2.kg.nodes[node].get("expression", {}).get("num_measured", numDefault)
                node2SD = kg2.kg.nodes[node].get("expression", {}).get("sd", sdDefault)
            
                ttestRes = scipy.stats.ttest_ind_from_stats(
                            mean1=node1Score, std1=node1SD, nobs1=node1num,
                            mean2=node2Score, std2=node2SD, nobs2=node2num)
                scoreSig = ttestRes.pvalue

            
            nodeDiffs[node] = (scoreDiff, scoreSig)
        
        return nodeDiffs
    
    
    def get_differential_graph(self, kg1:KGraph, kg2: KGraph, field="fc", rescoring=["geneset", "disease"]):
        
        ediffs = self._get_edge_fold_changes(kg1, kg2)
        ndiffs = self._get_node_fold_changes(kg1, kg2)

        #one should only keep non-score fields!!
        dkg = kg1.copy()
        dkg.kgraph_name = "{}_vs_{}".format(kg1.kgraph_name, kg2.kgraph_name)

        score_field = "{}_score".format(field)
        
        for node in dkg.kg.nodes:
            diffScore, diffSig = ndiffs[node]
            dkg.kg.nodes[node][field] = {"score": diffScore, "sig": diffSig}
            dkg.kg.nodes[node][score_field] = diffScore
        
        for edge in dkg.kg.edges:
            dkg.kg.edges[edge][field] = {"score": ediffs[edge]}
            dkg.kg.edges[edge][score_field] = ediffs[edge]


        # rescore for genesets, diseases, ...
        if not rescoring is None and len(rescoring) > 0:
            for node in dkg.kg.nodes:
                if len(set(rescoring).intersection(dkg.kg.nodes[node].get("type", []))) > 0:
                    scored_children = dkg._get_predecessors(node, ntype="gene") # this should be hierarchical?

                    
                    
                    node1ChildScores = [ kg1.kg.nodes[x].get("score", 0) for x in scored_children ]
                    node2ChildScores = [ kg2.kg.nodes[x].get("score", 0) for x in scored_children ]

                    child_scores = [ self._get_log_foldchange(n1s, n2s, log2=False) for n1s, n2s in zip(node1ChildScores, node2ChildScores) ]
                    new_score = np.log2(np.mean(child_scores))
                    
                    if min(len(node1ChildScores), len(node2ChildScores)) != 0:
                        scoreSig = scipy.stats.ttest_rel(node1ChildScores, node2ChildScores).pvalue
                    
                    else:
                        scoreSig = 1
                    
                    dkg.kg.nodes[node][field] = {"score": new_score, "sig": scoreSig}
                    dkg.kg.nodes[node][score_field] = new_score
                    
        return dkg
    
    def calculate_diffkg_list(self, exprKGs, base_case, cases):
        
        resdkgs = {}
        for case in cases:
            assert(case in exprKGs)
            
            if case == base_case:
                continue
            
            resdkgs[case] = self.get_differential_graph(exprKGs[case], exprKGs[base_case])
            
        return resdkgs
    
    
class NETSIM:

    def __init__(self, kg: KGraph, conf_score_accessor=lambda x: 1):

        def filter_func_nodes(node, k:KGraph):
            if k.node_type_overlap(node, ["gene"]):
                return True

            if k.node_type_overlap(node, ["geneset"]):
                if k.kg.nodes[node].get("source", "") == "GeneOntology":
                    return True

            return False
        
        self.kg = kg.filter_nodes(filter_function = filter_func_nodes)

        def filter_func_edges(edge, k:KGraph):
            e0t = k.node_types(edge[0])
            e1t = k.node_types(edge[1])

            if "gene" in e0t and "geneset" in e1t:
                return True
            if "gene" in e1t and "geneset" in e0t:
                return True
            if "geneset" in e0t and "geneset" in e1t:
                return True
            return False           
        

        self.orig_kg = kg
        self.conf_score_accessor=conf_score_accessor
        
        self.kg = self.kg.filter_edges(filter_function = filter_func_edges)
        self.kg_rev = self.kg.copy("reversed")

        self.kg_rev.kg = nx.DiGraph.reverse(self.kg_rev.kg)

        self.undirected_x = self.kg.kg.to_undirected()

        if not nx.is_directed_acyclic_graph(self.kg.kg):
            print("WARNING: the KG is not a DAG!")
        
        self.max_go_level = 10

        for n in self.kg.kg.nodes:
            self.max_go_level = max(self.max_go_level, max(self.kg.kg.nodes[n].get("go_level",0 ),self.kg.kg.nodes[n].get("go_depth",0)))

        self.precalculated_term_genes = {}

    def get_annotated_genes(self, ta):

        if ta in self.precalculated_term_genes:
            return self.precalculated_term_genes[ta]
            
        if ta is None:
            self.precalculated_term_genes[None] = self.kg.get_nodes("gene")
            return self.precalculated_term_genes[None]

        if self.kg.kg.nodes[ta].get("go_level", 10) <= 4:
            self.precalculated_term_genes[ta] = self.kg._get_predecessors(ta, ntype="gene", n=self.max_go_level)
            return self.precalculated_term_genes[ta]
        else:
            return self.kg._get_predecessors(ta, ntype="gene", n=self.max_go_level)
        


    def confidence(self, gi, gj):

        edges = []
        if (gi, gj) in self.orig_kg.kg.edges:
            edges.append( (gi, gj) )
        if (gj, gi) in self.orig_kg.kg.edges:
            edges.append( (gj, gi) )

        confs = []
        for edge in edges:
            confs.append( max(0, min(self.conf_score_accessor(self.orig_kg.kg.edges[edge]), 1)) )

        return sum(confs) / len(confs)

    
    def dij(self, gi, gj):
        if gi == gj:
            return 0

        if not ((gi, gj) in self.orig_kg.kg.edges) or not ((gj, gi) in self.orig_kg.kg.edges):
            return 1

        #gi -> gj or gj -> gi is in orig_kg !
        return 1.0-self.confidence(gi, gj)

    def functional_distance_genesets(self, ta, tb, ga, gb):

        # first sum
        dij_a = 0
        for gi in ga:
            pgb_dij = 1

            for gj in gb:
                pgb_dij = pgb_dij * self.dij(gi, gj)

                if pgb_dij == 0:
                    break

            dij_a += pgb_dij

        # second sum
        dij_b = 0
        for gi in gb:
            pga_dij = 1

            for gj in ga:
                pga_dij = pga_dij * self.dij(gi, gj)

                if pga_dij == 0:
                    break

            dij_b += pga_dij

        gab_elems = set(list(ga)+list(gb))
        dist = (dij_a+dij_b) / (2*len(gab_elems) - dij_a - dij_b)

        return dist

    def get_lca(self, ta, tb):
        return nx.lowest_common_ancestor(self.kg_rev.kg, ta, tb, default=None)

    def get_shortest_path(self, ta, tb):
        try:
            return nx.shortest_path(self.kg.kg, ta, tb)
        except:
            return None
    

    def get_path_constrained_annotation(self, ta, tb, p, verbose=False):

        a_shortest = self.get_shortest_path(ta, p)
        b_shortest = self.get_shortest_path(tb, p)

        if verbose:
            print("a_shortest", a_shortest)
            print("b_shortest", b_shortest)

        ga = self.get_annotated_genes(ta)
        gb = self.get_annotated_genes(tb)

        if verbose:
            print("ga", len(ga))
            print("gb", len(gb))


        #rational: if i fetch gp, then all intermediate annotated genes are also contained!
        gi = set()
        #if len(ga) > 0 and len(gb) > 0 and len(a_shortest) > 2:
        #    intermediate_nodes = a_shortest[1:-1]
        #    
        #    for ti in intermediate_nodes:
        #        ni = self.get_annotated_genes(ti)
        #        gi.update(ni)
        
        gj = set()
        #if len(ga) > 0 and len(gb) > 0 and len(b_shortest) > 2:
        #    intermediate_nodes = b_shortest[1:-1]
        #    
        #    for tj in intermediate_nodes:
        #        nj = self.get_annotated_genes(tj)
        #        gj.update(nj)

        gp = []
        if len(ga) > 0 and len(gb) > 0:
            gp = self.get_annotated_genes(p)
        
        if verbose:    
            print("gp", len(gp))
            print("gi", len(gi))
            print("gj", len(gj))
        
        totalset = set()
        totalset.update(ga)
        totalset.update(gb)
        totalset.update(gp)
        #totalset.update(gi)
        #totalset.update(gj)

        return ga, gb, gp, totalset

    def get_node_name(self, node):
        if node in self.kg.kg.nodes:
            return self.kg.kg.nodes[node]["name"]
        return "-/-"
    
    def functional_similarity(self, ta, tb, verbose=False):
       
        p = self.get_lca(ta, tb)

        if verbose:
            print("ta", ta, self.get_node_name(ta))
            print("tb", tb, self.get_node_name(tb))
            print("p", p, self.get_node_name(p))

        if p is None:
            # there is absolutely no similarity!
            return 0
        
        if ta == tb:
            return 1
        
        G = self.get_annotated_genes(None)

        Ga, Gb, Gp, Utatbp = self.get_path_constrained_annotation(ta, tb, p, verbose=verbose)

        if verbose:    
            print("LCA", p, self.kg.kg.nodes[p]["name"])
            print(self.get_shortest_path(ta, tb))

        if len(Ga) == 0 or len(Gb) == 0:
            return 0
        
        dtatb = self.functional_distance_genesets(ta,tb, Ga, Gb) ** 2
               
        ftatbp = dtatb * len(Utatbp) + (1-dtatb) * np.sqrt(len(Ga)*len(Gb))
        htatb = dtatb * len(G) + (1-dtatb) * max(len(Ga), len(Gb))
    
        if verbose:

            print("Utatbp=", len(Utatbp))
            print("dtatb=", dtatb)
            print("ftatbp=", ftatbp)
            print("htatb=", htatb)

        sim = ((2 * np.log(len(G))-2*np.log(ftatbp))/(2*np.log(len(G))-(np.log(len(Ga))+np.log(len(Gb)) ))) * (1-( (htatb/len(G)) * (len(Gp)/len(G)) ))

        if verbose:
            print("sim=", sim)

        return sim


    def get_go_level(self, ta):
        return self.orig_kg.kg.nodes[ta].get("go_level", 0)

    
    def lca_similarity(self, ta, tb, verbose=False):
       
        p = self.get_lca(ta, tb)

        if verbose:
            print("ta", ta, self.get_node_name(ta), self.get_go_level(ta))
            print("tb", tb, self.get_node_name(tb), self.get_go_level(tb))
            print("p", p, self.get_node_name(p), self.get_go_level(p))

        if p is None:
            # no lca identified
            return 0.0

        adist = self.get_go_level(ta)-self.get_go_level(p)
        bdist = self.get_go_level(tb)-self.get_go_level(p)

        total_dist = (adist + bdist) / 2.0

        if total_dist == 0:
            return 1.0

        sim = 1/total_dist

        return sim
    
    def get_relevant_goterms(self, groupKG, max_terms=None, verbose=False):
    
        allGeneSets = Counter()
        for node in groupKG.kg.nodes:
            genesetNeighbours = self.orig_kg.get_edges_to_type(node, "geneset")
            genesetNeighbours = [x for x in genesetNeighbours if self.orig_kg.kg.nodes[x].get("source", "") == "GeneOntology"]
    
            maxGoLevel = 0
            for x in genesetNeighbours:
                gol = self.orig_kg.kg.nodes[x].get("go_level", 0)
                maxGoLevel = max(gol, maxGoLevel)
    
            genesetNeighbours = [(x, self.orig_kg.kg.nodes[x].get("go_level", 0)) for x in genesetNeighbours if self.orig_kg.kg.nodes[x].get("go_level", 0) == maxGoLevel]
                            
            for n,l in genesetNeighbours:
                allGeneSets[n] += 1
    
        #print(allGeneSets)
        
        threshold = np.ceil(len(groupKG.kg.nodes)*0.01)
        retGeneSets = {}
        for n, ncount in allGeneSets.most_common(max_terms):
            if not allGeneSets[n] > threshold:
                continue
    
            numGenesetGenes = len(self.orig_kg._get_predecessors(n, ntype="gene"))
            numOccurences = allGeneSets[n]
    
            retGeneSets[n] = (numOccurences, numGenesetGenes)
    
        if verbose:
            print("Threshold", threshold)
            print("Selected genesets", len(retGeneSets))
        
        return retGeneSets


    def compare_modules_lca( self, skg1, skg2, max_terms=None, verbose=True):
    
        group1_terms = self.get_relevant_goterms(skg1, max_terms=max_terms, verbose=verbose)
        group2_terms = self.get_relevant_goterms(skg2, max_terms=max_terms, verbose=verbose)
    
        total_group1_terms = sum([group1_terms[x][0] for x in group1_terms])
        total_group2_terms = sum([group2_terms[x][0] for x in group2_terms])
    
        if verbose:
            print("group1 terms", len(group1_terms))
            print("group2 terms", len(group2_terms))
    
        
        comps = {}
        for termi in group1_terms:
            for termj in group2_terms:
                comps[(termi, termj)] = self.lca_similarity(termi, termj, verbose=verbose)
        
        weighted_sim = 0
        weight_total = 0
        for comp in comps:
            weight_g1 = group1_terms[comp[0]][0] / total_group1_terms
            weight_g2 = group2_terms[comp[1]][0] / total_group2_terms
            use_weight = weight_g1 * weight_g2
        
            if comps[comp] == 0:
                # assume different concepts!
                continue
        
            weight_total += use_weight
            weighted_sim = use_weight * comps[comp]
    
            if verbose:
                print(comp[0], self.orig_kg.kg.nodes[comp[0]]["name"], "<->", comp[1], self.orig_kg.kg.nodes[comp[1]]["name"], comps[comp], weight_g1, weight_g2, use_weight)
    
    
        if weight_total == 0:
            return 0
        
        simScore = weighted_sim/weight_total
            
        return simScore

    
    def compare_modules( self, skg1, skg2, max_terms=None, verbose=True):
    
        group1_terms = self.get_relevant_goterms(skg1, max_terms=max_terms, verbose=verbose)
        group2_terms = self.get_relevant_goterms(skg2, max_terms=max_terms, verbose=verbose)
    
        total_group1_terms = sum([group1_terms[x][0] for x in group1_terms])
        total_group2_terms = sum([group2_terms[x][0] for x in group2_terms])
    
        if verbose:
            print("group1 terms", len(group1_terms))
            print("group2 terms", len(group2_terms))
    
        
        comps = {}
        for termi in group1_terms:
            for termj in group2_terms:
                comps[(termi, termj)] = self.functional_similarity(termi, termj, verbose=verbose)
        
        weighted_sim = 0
        weight_total = 0
        for comp in comps:
            weight_g1 = group1_terms[comp[0]][0] / total_group1_terms
            weight_g2 = group2_terms[comp[1]][0] / total_group2_terms
            use_weight = weight_g1 * weight_g2
        
            if comps[comp] == 0:
                # assume different concepts!
                continue
        
            weight_total += use_weight
            weighted_sim = use_weight * comps[comp]
    
            if verbose:
                print(comp[0], self.orig_kg.kg.nodes[comp[0]]["name"], "<->", comp[1], self.orig_kg.kg.nodes[comp[1]]["name"], comps[comp], weight_g1, weight_g2, use_weight)
    
    
        if weight_total == 0:
            return 0
        
        simScore = weighted_sim/weight_total
            
        return simScore
    
class ModuleCompare:

    def __init__(self):
        pass

    # applications
    @classmethod
    def makeProgressBar(cls) -> progressbar.ProgressBar:
        return progressbar.ProgressBar(widgets=[
            progressbar.Bar(), ' ', progressbar.Percentage(), ' ', progressbar.AdaptiveETA()
            ])

    def draw_network(self, G, title=None):
        
        allEdgeWeights = [d["weight"] for (u, v, d) in G.edges(data=True) if d["weight"] > 0]
        #print(np.quantile(allEdgeWeights, [0, 0.25, 0.5, 0.75, 0.8, 1]))
        borderWeight = np.quantile(allEdgeWeights, 0.8)
        
        elarge = [(u, v) for (u, v, d) in G.edges(data=True) if d["weight"] >= borderWeight]
        esmall = [(u, v) for (u, v, d) in G.edges(data=True) if 0 < d["weight"] < borderWeight]

        largeNodes = set([x[0] for x in elarge]+[x[1] for x in elarge])

        fig, ax = plt.subplots(figsize=(12, 12))
        
        pos = nx.spring_layout(G, k=2, seed=7)  # positions for all nodes - seed for reproducibility
        
        # nodes
        nx.draw_networkx_nodes(G, pos, nodelist=largeNodes, node_size=200)
        nx.draw_networkx_nodes(G, pos, nodelist=[x for x in G.nodes if not x in largeNodes], node_size=100, alpha=0.5)
        
        # edges
        nx.draw_networkx_edges(G, pos, edgelist=elarge, width=6, edge_color="r")
        nx.draw_networkx_edges(
            G, pos, edgelist=esmall, width=6, alpha=0.1, edge_color="b", style="dashed"
        )
        
        # node labels
        nx.draw_networkx_labels(G, pos, font_size=10, font_family="sans-serif")
        # edge weight labels
        edge_labels = nx.get_edge_attributes(G, "weight")

        pedge_labels = {}
        
        for x in edge_labels:
            if not x in elarge:
                continue
            
            elabel = edge_labels[x]

            if elabel <= 0:
                continue

            if int(elabel) != elabel:
                elabel = "{:.3f}".format(elabel)

            pedge_labels[x] = elabel
        
        nx.draw_networkx_edge_labels(G, pos, pedge_labels)
    
        
        ax = plt.gca()
        
        if not title is None:
            ax.set_title(title)
        ax.margins(0.08)
        plt.axis("off")
        plt.tight_layout()
        plt.show()

    

    def network_compare_modules(self, inKGs, measure="jaccard"):
        
        assert (measure in ["jaccard", "sorensen"])

        similarityNetwork = nx.Graph()

        mods = [x for x in inKGs]

        modSims = {}

        for i, inMod in enumerate(mods):

            similarityNetwork.add_node(inMod)

            for j in range(i+1, len(mods)):
                oMod = mods[j]
                if inMod == oMod:
                    continue

                if not oMod in similarityNetwork.nodes:
                    similarityNetwork.add_node(oMod)

                inSet = set(inKGs[inMod].kg.nodes)
                oSet = set(inKGs[oMod].kg.nodes)

                node_overlap = inSet.intersection(oSet)
                node_union = inSet.union(oSet)

                similarity = None
                if measure == "jaccard":
                    similarity = len(node_overlap) / len( node_union )
                elif measure == "sorensen":
                    similarity = (2.0*len(node_overlap)) / ( len(inSet) + len(oSet) )

                #allOverlaps.append(len(node_overlap))

                modSims[(inMod, oMod)] = similarity
                similarityNetwork.add_edge(inMod, oMod, weight=similarity)                

        self.draw_network(similarityNetwork, title="Module Overlaps")
        return modSims

    def network_compare_lca(self, inKGs, max_terms=None, fullKG=None, ns=None):
        
        assert(not (fullKG is None and ns is None))

        if ns is None:           
            ns = NETSIM(fullKG)

        modNames = [x for x in inKGs]
        
        modComps = {}
        bar = self.makeProgressBar()
        allComparisons = [(i,j) for i in range(0, len(modNames)) for j in range(i+1, len(modNames))]
        for i,j in bar(allComparisons):      
            
            signame1 = modNames[i]
            signame2 = modNames[j]
            
            sig1 = inKGs[signame1]
            sig2 = inKGs[signame2]
            modComps[(signame1, signame2)] = ns.compare_modules_lca(sig1, sig2, max_terms=max_terms, verbose=False)
                
        G = nx.Graph()
        for edge in modComps:
            G.add_edge( edge[0], edge[1], weight=modComps[edge])
        
        self.draw_network(G)
        
        return modComps
    
    def network_compare_netsim(self, inKGs, max_terms=None, fullKG=None, ns=None):

        assert(not (fullKG is None and ns is None))


        if ns is None:           
            ns = NETSIM(fullKG)
            
        modNames = [x for x in inKGs]
        
        modComps = {}
        bar = self.makeProgressBar()
        allComparisons = [(i,j) for i in range(0, len(modNames)) for j in range(i+1, len(modNames))]
        for i,j in bar(allComparisons):
            #print(datetime.fromtimestamp(time.time()))
        
            
            signame1 = modNames[i]
            signame2 = modNames[j]
            
            sig1 = inKGs[signame1]
            sig2 = inKGs[signame2]
            modComps[(signame1, signame2)] = ns.compare_modules(sig1, sig2, max_terms=max_terms, verbose=False)
        
                #print(datetime.fromtimestamp(time.time()))
        
        G = nx.Graph()
        for edge in modComps:
            G.add_edge( edge[0], edge[1], weight=modComps[edge])
        
        self.draw_network(G)
        
        return modComps
    

    def module_umap(self, sigKG, kg:KGraph):

        nodeList = [x for x in kg.kg.nodes]
        cluster_names = [x for x in sigKG]
        cluster_labels = [x.split("_")[0] for x in sigKG]
        
        # rows, cols
        moduleDF = pd.DataFrame(np.zeros((len(nodeList), len(sigKG)), dtype=np.uint8))
        moduleDF.index = nodeList
        moduleDF.columns = cluster_names
                

        for comm in cluster_names:
            commNodes = [x for x in sigKG[comm].kg.nodes]

            moduleDF.loc[commNodes, comm] = 1
        
        #print(moduleDF.head())
        #print(moduleDF.sum(axis=0))

        from sklearn import decomposition
        from sklearn import preprocessing

        scalar = preprocessing.StandardScaler()
        standardized_data = scalar.fit_transform(moduleDF.T)

        pca = decomposition.PCA(n_components=2)
        pca_data = pca.fit_transform(standardized_data)

        fig, ax = plt.subplots(figsize=(12,12))
        
        df = pd.DataFrame(pca_data)
        df.index = cluster_names
        df.columns = ["pca0", "pca1"]
        df["Module"] = cluster_names
        df["label"] = cluster_labels
        df['label'] = df['label'].astype('category')

        
        
        df.plot.scatter("pca0", "pca1", c="label", colormap='viridis', ax=ax)

        texts = []
        for x, y, s in zip(df["pca0"], df["pca1"], df["Module"]):
            texts.append(plt.text(x, y, s,fontproperties={'size' : 6}))
            
        adjust_text(texts, force_points=0.2, force_text=0.2,
            expand_points=(1, 1), expand_text=(1, 1),
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
        plt.show()
        plt.close()
        
        #pca = decomposition.PCA(n_components=30)
        #pca_data = pca.fit_transform(standardized_data)
        #print("pca data", pca_data.shape)
        
        #import umap
        #reducer = umap.UMAP(
        #                    n_neighbors=5,
        #                    min_dist=0.8
        #                   )
        #embedding = reducer.fit(pca_data)
        #print("umap embedding", embedding)

        #import umap.plot      
        #umap.plot.points(embedding, labels=np.array(cluster_labels), background='black')
        #umap.plot.connectivity(embedding, show_points=True, edge_bundling='hammer')

    def plot_dendrogram(self, simDict):
        allModules = natsorted(set([x[0] for x in simDict]+[x[1] for x in simDict]))
        simMatrix = np.zeros( (len(allModules), len(allModules)) )

        for mod1, mod2 in simDict:
            simMatrix[ allModules.index(mod1), allModules.index(mod2) ] = simDict[(mod1, mod2)]
            simMatrix[ allModules.index(mod2), allModules.index(mod1) ] = simDict[(mod1, mod2)]

            

        simVec = scipy.spatial.distance.squareform(simMatrix)
        linkage =  scipy.cluster.hierarchy.linkage(1 - simVec)
        dendro  =  scipy.cluster.hierarchy.dendrogram(linkage, labels=allModules)
        
        plt.gca().set_xticklabels(plt.gca().get_xticklabels(), rotation = 45, ha="right")
        plt.show()

class CRankExplorer:

    def __init__(self):
        pass


    def calculate_connectivity(self, kg:KGraph):

        nodeDegrees = [x[1] for x in nx.degree(kg.kg)]
        return np.mean(nodeDegrees)


    def calculate_density(self, communityKG:KGraph, edge_score_accessor):

        allEdgeScores = []
        for edge in communityKG.kg.edges:
            edgeScore = abs(edge_score_accessor(communityKG.kg.edges[edge]))
            allEdgeScores.append(edgeScore)

        #print("density", np.quantile(allEdgeScores, [0, 0.2, 0.5, 0.8, 1.0]))
        return np.mean(allEdgeScores)

    def calculate_boundary(self, community, fullKG:KGraph, edge_score_accessor, density):

        allOutEdgeScores = []
        for node in community:
            for edge in fullKG.kg.out_edges(node):
                edgeScore = abs(edge_score_accessor(fullKG.kg.edges[edge]))
                allOutEdgeScores.append(edgeScore)
        
        #print("boundary", np.quantile(allOutEdgeScores, [0, 0.2, 0.5, 0.8, 1.0]))
        return np.mean(allOutEdgeScores)
                

    def calculate_allegiance(self, community, fullKG:KGraph):

        prefs = []
        for node in community:

            commEdgesCount = 0
            outEdgesCount = 0

            outEdges = fullKG.kg.out_edges(node)
            
            if len(outEdges) == 0:
                continue
            
            for edge in outEdges:
                if edge[0] in community and edge[1] in community:
                    commEdgesCount += 1
                else:
                    outEdgesCount += 1

            prefWithinCommunity = commEdgesCount / (commEdgesCount+outEdgesCount)
            prefs.append(prefWithinCommunity)

        #print("allegiance", np.quantile(prefs, [0, 0.2, 0.5, 0.8, 1.0]))
        return np.mean(prefs)

    def evaluate_community(self, community, fullKG: KGraph, edge_score_accessor=lambda x: x["score"]):#, node_score_accessor=lambda x: x["score"]):

        communityKG = fullKG.subset_kg(community)
        #undirKG:nx.Graph = communityKG.kg.to_undirected()

        connectivityScore = self.calculate_connectivity(communityKG)
        densityScore = self.calculate_density(communityKG, edge_score_accessor)
        boundaryScore = self.calculate_boundary(community, fullKG, edge_score_accessor, densityScore)
        allegianceScore = self.calculate_allegiance(community, fullKG)


        return {
            "connectivity_score": connectivityScore,
            "density_score": densityScore,
            "boundary_score": boundaryScore,
            "allegiance_score": allegianceScore,
        }


    def evaluate_communities(self, sigKGs, fullKGs):

        commScores = []
        for comm in sigKGs:
        
            mainGraph = comm.split("_mod_")[0]
            #print(comm, mainGraph)

            mainKG = fullKGs.get(mainGraph, None)

            if mainKG is None:
                print("Not a valid graph", mainGraph)
                continue
          
            commScore = self.evaluate_community(set(sigKGs[comm].kg.nodes), mainKG, edge_score_accessor=lambda x: x["fc_score"])
            commScore["module"] = comm
            commScores.append(commScore)
        
        commScoreDF = pd.DataFrame(commScores)

        for score in commScoreDF.columns:
            if not score.endswith("_score"):
                continue
                
            commScoreDF["{}_rank".format(score)] = commScoreDF[score].rank(ascending=False)    
                
        rankCols = [x for x in commScoreDF.columns if x.endswith("_rank")]
        commScoreDF["rank_mean"] = commScoreDF[rankCols].mean(axis=1)
        
        commScoreDF.sort_values("rank_mean", inplace=True)
        return commScoreDF


    def plot_scores( self, scoreDF ):


        moduleDF = scoreDF[[x for x in scoreDF.columns if x.endswith("_score")]]

        scalar = preprocessing.StandardScaler()
        standardized_data = scalar.fit_transform(moduleDF)

        pca = decomposition.PCA(n_components=2)
        pca_data = pca.fit_transform(standardized_data)


        fig, ax = plt.subplots(figsize=(12,12))
        
        df = pd.DataFrame(pca_data)
        
        df.index = scoreDF["module"]
        df.columns = ["pc0", "pc1"]
        df["Module"] = list(scoreDF["module"])
        df["label"] = [x.split("_")[0] for x in df["Module"]]
        df['label'] = df['label'].astype('category')
        
        df.plot.scatter("pc0", "pc1", c="label", colormap='viridis', ax=ax)

        texts = []
        for x, y, s in zip(df["pc0"], df["pc1"], df["Module"]):
            texts.append(plt.text(x, y, s,fontproperties={'size' : 6}))
            
        adjust_text(texts, force_points=0.2, force_text=0.2,
            expand_points=(1, 1), expand_text=(1, 1),
            arrowprops=dict(arrowstyle="-", color='black', lw=0.5))
        plt.show()
        plt.close()


