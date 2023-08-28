

import os, sys
import urllib
import zipfile
import gzip
import pandas as pd
import numpy as np

from collections import defaultdict
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag

import networkx as nx


def download_and_unzip(download_url_link, dir_path, zipped_filename,destination_dir_name, unzip=True):
    #https://www.tutorialsbuddy.com/download-and-unzip-a-zipped-file-in-python
    print("Download starting")

    urllib.request.urlretrieve(
        download_url_link, os.path.join(dir_path, zipped_filename)
    )
    print("Download complete")

    if unzip:
        print("unzipping file starting")
    
        if zipped_filename.endswith(".zip"):
            with zipfile.ZipFile(os.path.join(dir_path, zipped_filename), "r") as zip_file:
                zip_file.extractall(os.path.join(dir_path, destination_dir_name))
        elif zipped_filename.endswith(".gz"):
            print("zipfile")
            with gzip.GzipFile(os.path.join(dir_path, zipped_filename), "rb") as zip_file:
                with open(os.path.join(dir_path, destination_dir_name, zipped_filename.replace(".gz", "")), "wb") as fout:
                    fout.write(zip_file.read())
        else:
            raise NotImplementedError("NO CASE")
            
    
    print("unzipping complete")
    
    
    
    

def load_go(kg: nx.DiGraph, data_dir, source="GeneOntology", interaction_harmonize = {
        'NOT': "relevant_in",
        'acts_upstream_of': "relevant_in",
        'acts_upstream_of_negative_effect': "relevant_in",
        'acts_upstream_of_or_within': "relevant_in",
        'acts_upstream_of_or_within_negative_effect': "relevant_in",
        'acts_upstream_of_or_within_positive_effect': "relevant_in",
        'acts_upstream_of_positive_effect': "relevant_in",
        'colocalizes_with': "relevant_in",
        'contributes_to': "activates",
        'enables': "activates",
        'involved_in': "relevant_in",
        'is_active_in': "activates",
        'located_in': "relevant_in",
        'part_of': "relevant_in"
        }):
    
    if not os.path.exists(os.path.join(data_dir, "goa_human.gaf")):
        download_and_unzip("http://geneontology.org/gene-associations/goa_human.gaf.gz", ".", os.path.join(data_dir, "goa_human.gaf.gz"), data_dir)
            
    if not os.path.exists(os.path.join(data_dir, "go-basic.obo")):
        download_and_unzip("http://geneontology.org/ontology/go-basic.obo", ".", os.path.join(data_dir, "go-basic.obo"), data_dir, unzip=False)

    genenamesURL = 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'
    if not os.path.exists(os.path.join(data_dir, "hgnc_annot.tsv")):
        download_and_unzip(genenamesURL, ".", os.path.join(data_dir, "hgnc_annot.tsv"), data_dir, unzip=False)

    ogaf = GafReader(os.path.join(data_dir, "goa_human.gaf"))
    obodag = GODag(os.path.join(data_dir, "go-basic.obo"))
            

    # read hgnc <-> uniprot conversion
    hgncDF = pd.read_csv(os.path.join(data_dir, "hgnc_annot.tsv"), sep="\t")
    uniprot2hgnc = defaultdict(set)
    all_genes = set()

    for ri, row in hgncDF.iterrows():
        
        status = row["Status"]
        
        if status == "Symbol Withdrawn":
            continue
        
        symbol = row["Approved symbol"]    
        uniprot= row["UniProt ID(supplied by UniProt)"]
        
        all_genes.add(symbol)
        
        if pd.isna(uniprot):
            continue

        uniprot2hgnc[uniprot].add(symbol)
        
        
    #fetch go2gene associations
    go2gene = defaultdict(set)

    for assoc in ogaf.get_associations():

        geneID = assoc.DB_ID
        goID = assoc.GO_ID
        
        if not geneID in uniprot2hgnc:
            geneID = assoc.DB_Symbol
            geneIDs = [geneID]
            
            if not geneID in all_genes:
                #print(geneID, assoc.Taxon, assoc)
                continue
        else:
            geneIDs = uniprot2hgnc[geneID]
        
        for geneSym in geneIDs:
            go2gene[goID].add((geneSym, tuple(assoc.Qualifier)))


    #
    ## Ready to fill graph
    #
    
    # add all genes
    for gene in all_genes:
        kg.add_node(gene, type="gene", score=0)
        
    #add GO nodes with attributes
    for goEntry in obodag:
        
        termID = goEntry
        termObj = obodag[goEntry]
        
        if termObj.is_obsolete:
            continue
        
        termName = termObj.name
        termNS = str(termObj.namespace)
        
        kg.add_node(termID, id=termID, name=termName, type="geneset", ns=termNS, score=0, source=source)
        
    #add GO edges with attributes
    for goEntry in obodag:
        
        termID = goEntry
        termObj = obodag[goEntry]
        
        for child in termObj.children:
            if child in kg.nodes and not (child, goEntry) in kg.edges:
                kg.add_edge(child, goEntry, type="part_of", score=0, source=source)
                
        for parent in termObj.parents:
            if parent.id in kg.nodes and not (goEntry, parent.id) in kg.edges:
                kg.add_edge(goEntry, parent.id, type="part_of", score=0, source=source)
                
                
                
    all_interactions = set()
    for goID in list(go2gene):   
        for gene, interaction in go2gene[goID]:
            for x in interaction:
                all_interactions.add(x)

    # add edges gene -> GO
    for goID in list(go2gene):
        
        for gene, interaction in go2gene[goID]:
            if goID in kg.nodes:
                kg.add_edge( gene, goID, type=interaction_harmonize[interaction[0]], go_interaction = interaction[0], source=source)

    return kg
            
def load_omnipath(kg: nx.DiGraph, data_dir, source="omnipath"):
    
    #
    ## Checking Files
    #
    omnipathDB = os.path.join(data_dir, "omnipath.tsv")
    hgncTranslationDB = os.path.join(data_dir, "omnipath_hgnc_uniprot.tsv")
    
    if not os.path.exists(hgncTranslationDB):
        hgncURL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

        hgncTranslation = pd.read_csv(hgncURL, sep="\t")
        hgncTranslation.to_csv(hgncTranslationDB, sep="\t")

    if not os.path.exists(omnipathDB):
        
        import omnipath as op
        opd = op.interactions.AllInteractions().get()
        
        opd.to_csv(omnipathDB, sep="\t")
        
    
    
    #
    ## Adding nodes
    # 
    hgncTranslation = pd.read_csv(hgncTranslationDB, sep="\t")
    opd  = pd.read_csv(omnipathDB, sep="\t")
        
    uniprot2gene = defaultdict(set)

    for ri, row in hgncTranslation.iterrows():
        
        hgncGene = row["Approved symbol"]
        uniprotID = row["UniProt ID(supplied by UniProt)"]
        
        if pd.isna(uniprotID):
            continue
        
        uniprot2gene[uniprotID].add(hgncGene)
        
            
    
    #
    ## Adding edges
    #   
    for ri, row in opd.iterrows():
        
        srcUID = row["source"]
        tgtUID = row["target"]
        
        srcGenes = uniprot2gene.get(srcUID, set())
        tgtGenes = uniprot2gene.get(tgtUID, set())
        
        if len(srcGenes) == 0 or len(tgtGenes) == 0:
            continue
            
        for src in srcGenes:
            if not src in kg.nodes:
                kg.add_node(src, type="gene", source=source)
            
        for tgt in tgtGenes:
            if not tgt in kg.nodes:
                kg.add_node(tgt, type="gene", source=source)
            
            
    interactionTypes = {} # stimulation, inhibition
    interactionTypes[(False, False)] = "interacts"
    interactionTypes[(False, True)] = "represses"
    interactionTypes[(True, False)] = "activates"
    interactionTypes[(True, True)] = "interacts"

    ignoredCount = 0

    for ri, row in opd.iterrows():
            
        srcUID = row["source"]
        tgtUID = row["target"]
        
        srcGenes = uniprot2gene.get(srcUID, set())
        tgtGenes = uniprot2gene.get(tgtUID, set())
        
        if len(srcGenes) == 0 or len(tgtGenes) == 0:
            continue
        
        # consensus_direction -> all resources have same direction
        # consensus_stimulation -> all resources show this as stimulation
        # consensus_inhibition -> all resources show this as consensus_inhibition
        
        # it must have a consensus direction
        if not (row["consensus_direction"]):
            ignoredCount += 1
            continue
        
        # it must either have consensus stimulation or inhibition
        if not (row["consensus_stimulation"] or row["consensus_inhibition"]):
            ignoredCount += 1
            continue

        interactionType = interactionTypes[ (row["consensus_stimulation"], row["consensus_inhibition"]) ]
        
        omnipathEvidences = str(row["references"]).strip().split(";")
        omnipathType = row["type"]
        
        
        for src in srcGenes:
            
            for tgt in tgtGenes:
                
                attrDict = {
                    "type": interactionType, "source": source, "omnipath_evidences": omnipathEvidences, "omnipath_type": omnipathType
                }
                
                kg.add_edge( src, tgt, **attrDict )
                
                if not row["is_directed"]:
                    kg.add_edge(tgt, src, **attrDict)
                    
    return kg
        
        
def load_STRING(kg: nx.DiGraph, data_dir, mart_file="oct2014_mart_export.txt", source="STRING", use_evidences = ['fusion', 'coexpression','experiments','database','textmining']):


    if not os.path.exists(os.path.join(data_dir, "9606.protein.links.full.txt.gz")):
        download_and_unzip("https://stringdb-static.org/download/protein.links.full.v11.5/9606.protein.links.full.v11.5.txt.gz", ".", os.path.join(data_dir, "9606.protein.links.full.txt.gz"), data_dir)


    df = pd.read_csv(os.path.join(data_dir, "9606.protein.links.full.txt"), sep=" ")
    
    subdf = df[["protein1", "protein2"]+use_evidences]
    subdf["score"] = subdf[use_evidences].max(axis=1)/1000
    subdf = subdf[subdf.score > 0]
    subdf = subdf[subdf.score >= 0.7]

    
    all_ensp_proteins = set()
    allStringProts = set(subdf.protein1)
    allStringProts.update(subdf.protein2)

    for x in allStringProts:
        all_ensp_proteins.add(x.split(".")[1])
        
        
    martDF = pd.read_csv(os.path.join(data_dir, mart_file), sep="\t")

    ensemblProt2Gene = defaultdict(set)
    
    martEmptyProt = ~pd.isna(martDF["Ensembl Protein ID"])
    martEmptyHgnc = ~pd.isna(martDF["HGNC symbol"])
    
    martFilter = np.where(martEmptyProt & martEmptyHgnc)
    
    for ri, row in martDF.loc[martFilter].iterrows():
        
        protid = "9606.{}".format(row["Ensembl Protein ID"])
        geneid = row["HGNC symbol"]

        ensemblProt2Gene[protid].add(geneid)
        
        
    for ri, row in subdf.iterrows():
        src = row["protein1"]
        tgt = row["protein2"]
        
        if not src in ensemblProt2Gene:
            continue
        if not tgt in ensemblProt2Gene:
            continue
        
        src = ensemblProt2Gene[src]
        tgt = ensemblProt2Gene[tgt]
        
        string_scores = {}
        for sc in use_evidences + ["score"]:
            string_scores[sc] = row[sc]
                    
        for s in src:
            for t in tgt:
                
                if s == "PGM3":
                    print(s, t, string_scores)
                
                # fixes empty gene symbols in mart file ...
                if s == "":
                    s = src
                if t == "":
                    t = tgt
                
                if not s in kg.nodes:
                    kg.add_node(s, type="gene", source=source)
                if not t in kg.nodes:
                    kg.add_node(t, type="gene", source=source)
                
                kg.add_edge( s, t, type="interacts", string_scores=string_scores, source=source, score=0 )
                
    return kg


def load_opentargets(kg: nx.DiGraph, data_dir, source="opentargets", min_disease_association_score=0.8):
    
    if not os.path.exists(os.path.join(data_dir, "opentargets_knowndrugs.tsv")):
        print("Missing OpenTargets data frame knowndrugs. Prepare data according to Jupyter Notebook:", "opentargets_knowndrugs")
        exit(-1)

    if not os.path.exists(os.path.join(data_dir, "opentargets_disease_associations.tsv")):
        print("Missing OpenTargets data frame knowndrugs. Prepare data according to Jupyter Notebook:", "opentargets_disease_associations")
        exit(-1)

    ot_drugs = pd.read_csv("../data/opentargets_knowndrugs.tsv", sep="\t")
    ot_disease = pd.read_csv("../data/opentargets_disease_associations.tsv", sep="\t")
    
    
    #
    ## Adding Disease Nodes
    #
    
    for ri, row in ot_disease.iterrows():
        gene = row["targetSymbol"]
        disease_id = row["diseaseId"].replace("_", ":")
        disease_label = row["diseaseLabel"]
        
        if not gene in kg.nodes:
            kg.add_node(gene, type="gene", source=source)
            
        if not disease_id in kg.nodes:
            kg.add_node(disease_id, type="disease", name=disease_label, source=source)
            
            
    #
    ## Adding Drug Nodes
    #
    
    for ri, row in ot_drugs.iterrows():
        drug = row["drugId"]
        disease_id = row["diseaseId"].replace("_", ":")
        gene = row["targetGeneSymbol"]
        
        if not gene in kg.nodes:
            kg.add_node(gene, type="gene", source=source)
            
        if not drug in kg.nodes:
            kg.add_node(drug, type="drug", source=source)
            
        if not disease_id in kg.nodes:
            kg.add_node(disease_id, type="disease", source=source)
            
    #
    ## Adding Disease Edges
    #
    for ri, row in ot_disease.iterrows():
            
        gene = row["targetSymbol"]
        disease_id = row["diseaseId"].replace("_", ":")
        
        disease_score = float(row["datatypeHarmonicScore"])
        disease_evidences = row["datatypeEvidenceCount"]
        disease_evidence_source = row["datatypeId"]
        
        if disease_score < min_disease_association_score:
            continue
        
        disease_data = {
            "disease_score": disease_score,
            "disease_evidences": disease_evidences,
            "disease_evidence_source": disease_evidence_source,
            "type": "interacts",
            "source": source
        }
        
        kg.add_edge( gene, disease_id, **disease_data )
        
        
    #
    ## Adding Drug Edges
    #
    for ri, row in ot_drugs.iterrows():
        #drugId	targetId	diseaseId	status	diseaseName	targetGeneSymbol	targetGeneName

            
        drug = row["drugId"]
        disease_id = row["diseaseId"].replace("_", ":")
        drug_target = row["targetGeneSymbol"]
        
        evidence_status = row["status"]
        
        drug_disease_data = {
            "evidence_status": evidence_status,
            "source": source,
            "type": "affects"
        }
        
        kg.add_edge( drug, disease_id, **drug_disease_data )
        kg.add_edge( drug, drug_target, type="targeted_by", source=source )
        
    return kg
        
        
def load_reactome(kg: nx.DiGraph, data_dir, source="reactome"):

    reactomeFile = os.path.join(data_dir,"ReactomePathways.gmt")
    
    if not os.path.exists(reactomeFile):
        download_and_unzip("https://reactome.org/download/current/ReactomePathways.gmt.zip", ".", os.path.join(data_dir,"ReactomePathways.gmt.zip"), data_dir)
        
    
    with open(reactomeFile) as fin:
    
        for line in fin:
            line = line.strip().split("\t")
            
            reactomeLabel = line[0]
            reactomeID = line[1]
            
            reactomeGenes = [x for x in line[2:] if x[0].isupper()]
                
            if not reactomeID in kg.nodes:
                kg.add_node(reactomeID, id=reactomeID, type="geneset", name=reactomeLabel, score=0, source=source)
            
            
            for gene in reactomeGenes:
                if not gene in kg.nodes:
                    kg.add_node(gene, type="gene", score=0)
                    
                kg.add_edge(gene, reactomeID, type="part_of", score=0, source=source)
    
    return kg