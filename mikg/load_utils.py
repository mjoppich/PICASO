

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


def download_and_unzip(download_url_link, dir_path, zipped_filename,destination_dir_name, unzip=True, force_zip=False, force_gz=False):
    #https://www.tutorialsbuddy.com/download-and-unzip-a-zipped-file-in-python
    print("Download starting")

    urllib.request.urlretrieve(
        download_url_link, os.path.join(dir_path, zipped_filename)
    )
    print("Download complete")

    if unzip:
        print("unzipping file starting")
    
        if (zipped_filename.endswith(".zip") and not force_gz) or force_zip:
            with zipfile.ZipFile(os.path.join(dir_path, zipped_filename), "r") as zip_file:
                zip_file.extractall(os.path.join(dir_path, destination_dir_name))
        elif (zipped_filename.endswith(".gz") and not force_zip) or force_gz:
            print("gzfile")
            with gzip.GzipFile(os.path.join(dir_path, zipped_filename), "rb") as zip_file:
                destname = os.path.join(dir_path, destination_dir_name, os.path.basename(zipped_filename).replace(".gz", ""))
                print(destname)
                with open(destname, "wb") as fout:
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
        if gene in kg.nodes:
            kg.nodes[gene]["type"].add("gene")
        else:
            kg.add_node(gene, type=set(["gene"]), score=0, name=gene, source=source)
        
    #add GO nodes with attributes
    for goEntry in obodag:
        
        termID = goEntry
        termObj = obodag[goEntry]
        
        if termObj.is_obsolete:
            continue
        
        termName = termObj.name
        termNS = str(termObj.namespace)
        
        kg.add_node(termID, id=termID, name=termName, type=set(["geneset"]), ns=termNS, score=0, source=source)
        
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
            if src in kg.nodes:
                kg.nodes[src]["type"].add("gene")
            else:
                kg.add_node(src, type=set(["gene"]), source=source, name=src, score=0)
            
        for tgt in tgtGenes:
            if tgt in kg.nodes:
                kg.nodes[tgt]["type"].add("gene")
            else:
                kg.add_node(tgt, type=set(["gene"]), source=source, name=tgt, score=0)
            
            
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
                
                if s in kg.nodes:
                    kg.nodes[s]["type"].add("gene")
                else:
                    kg.add_node(s, type=set(["gene"]), source=source, name=s, score=0)
                
                if t in kg.nodes:
                    kg.nodes[t]["type"].add("gene")
                else:
                    kg.add_node(t, type=set(["gene"]), source=source, name=t, score=0)
                                    
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
            kg.add_node(gene, type=set(["gene"]), source=source, name=gene, score=0)
        else:
            kg.nodes[gene]["type"].add("gene")
            
        if not disease_id in kg.nodes:
            kg.add_node(disease_id, type=set(["disease"]), name=disease_label, source=source)
        else:
            kg.nodes[disease_id]["type"].add("disease")
            
            
    #
    ## Adding Drug Nodes
    #
    
    for ri, row in ot_drugs.iterrows():
        drug = row["drugId"]
        disease_id = row["diseaseId"].replace("_", ":")
        diseaseName = row["diseaseName"]
        gene = row["targetGeneSymbol"]
        geneName = row["targetGeneName"]
        
        drugName = row["prefName"]
        drugType = row["drugType"]
        
        if not gene in kg.nodes:
            kg.add_node(gene, type=set(["gene"]), name=geneName, source=source, score=0)
        else:
            kg.nodes[gene]["type"].add("gene")
            
        if not drug in kg.nodes:
            kg.add_node(drug, type=set(["drug"]), source=source, name=drugName, drug_type=drugType)
        else:
            kg.nodes[drug]["type"].add("drug")
                
        if not disease_id in kg.nodes:
            kg.add_node(disease_id, type=set(["disease"]), name=diseaseName, source=source)
        else:
            kg.nodes[disease_id]["type"].add("disease")
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
                kg.add_node(reactomeID, id=reactomeID, type=set(["geneset"]), name=reactomeLabel, score=0, source=source)
            else:
                kg.nodes[reactomeID]["type"].add("geneset")
            
            for gene in reactomeGenes:
                if not gene in kg.nodes:
                    kg.add_node(gene, type=set(["gene"]), name=gene, source=source, score=0)
                else:
                    kg.nodes[gene]["type"].add("gene")
                    
                kg.add_edge(gene, reactomeID, type="part_of", score=0, source=source)
    
    return kg


def load_npinter(kg: nx.DiGraph, data_dir, source="npinter5"):

    npinterFile = os.path.join(data_dir,"interaction_NPInterv5.txt")
    
    if not os.path.exists(npinterFile):
        download_and_unzip("http://bigdata.ibp.ac.cn/npinter5/download/file/interaction_NPInterv5.txt.gz", ".", os.path.join(data_dir,"interaction_NPInterv5.txt.gz"), data_dir, force_zip=True)


    df = pd.read_csv(npinterFile, sep="\t")
    df = df[df.organism == "Homo sapiens"].copy()

    colMap = {
    'binding; regulatory': 'binding;regulatory',
    'binding;': "binding"
    }
    
    for x in set(df["class"]):
        if not x in colMap:
            colMap[x]=x


    df["class"]=df["class"].map(colMap)

    type2kgtype = {
                    'miRNA': ["ncRNA", "miRNA"],
                    'lncRNA': ["ncRNA", "lncRNA"],
                    'snoRNA': ["ncRNA"],
                    'mRNA': ["gene"],
                    'snRNA': ["ncRNA"],
                    'ncRNA': ["ncRNA"],
                    'pseudogene': ["gene"],
                    'Pseudogene': ["gene"],
                    'circRNA': "ncRNA",
                    'protein': ["gene"],
                    'Protein': ["gene"],
                    'sRNA': ["ncRNA"],
                    'vtRNAs': ["ncRNA"],
                    'piRNAs': ["ncRNA"],
                    "DNA": None,
                    "TF": ["gene", "TF"]
                  }

    for ri, row in df.iterrows():
    
        ncName = row["ncName"]
        ncID = row["ncID"]
        ncType = row["ncType"]
    
        tarName = row["tarName"]
        tarID = row["tarID"]
        tarType = row["tarType"]
    
        intClass = row["class"]
        intID = row["interID"]
        intSource = row["datasource"]
    
        if type(ncID) == float:
            if np.isnan(ncID):
                ncID = ncName
                
        if type(tarID) == float:
            if np.isnan(tarID):
                tarID = tarName
    
    
        if ncName.startswith("hsa-"):
            ncName = ncName[4:]
    
        kgNcTypes = type2kgtype[ncType]
        kgTarTypes = type2kgtype[tarType]
        
        if kgNcTypes is None or kgTarTypes is None:
            continue
    
        if not ("gene" in kgNcTypes or "gene" in kgTarTypes):
            continue
            
        #if None in kgNcTypes or None in kgTarTypes:
        #    continue

        if not ncName in kg.nodes:
            kg.add_node(ncName, id=ncID, name=ncName, type=kgNcTypes, biotype=ncType, score=0, source=source)
        else:
            kg.nodes[ncName]["type"].update(kgNcTypes)
    
        if not tarName in kg.nodes:
            kg.add_node(tarName, id=tarID, name=tarName, type=kgTarTypes, biotype=tarType, score=0, source=source)
        else:
            kg.nodes[tarName]["type"].update(kgTarTypes)
    
        kg.add_edge(ncName, tarName, type="interacts", score=0, source="npinter5", source_type=intClass, source_id=intID)
    
    return kg



def load_kegg(kg: nx.DiGraph, data_dir, source="kegg",
              hallmark_genesets="kegg_gmts/human/c1.all.v2023.2.Hs.symbols.gmt",
              curated_genesets="kegg_gmts/human/c2.all.v2023.2.Hs.symbols.gmt"
              ):
    
    hallmark_genesets = os.path.join(data_dir, hallmark_genesets)
    if not hallmark_genesets is None and not os.path.exists(hallmark_genesets):
        print("Missing KEGG pathways: hallmark genesets")
        exit(-1)

    curated_genesets = os.path.join(data_dir, curated_genesets)
    if not curated_genesets is None and not os.path.exists(os.path.join(data_dir, curated_genesets)):
        print("Missing KEGG pathways: curated genesets")
        exit(-1)

    def load_kegg_gmt(infile, kg):
        
        with open(infile, "r") as fin:
            
            for line in fin:
                line = line.strip().split("\t")
                
                pathway_id = line[0]
                pathway_genes = line[2:]
                
                
                if not pathway_id in kg.nodes:
                    kg.add_node(pathway_id, id=pathway_id, type=set(["geneset"]), name=pathway_id, score=0, source=source)
                else:
                    kg.nodes[pathway_id]["type"].add("geneset")
                
                for gene in pathway_genes:
                    if not gene in kg.nodes:
                        kg.add_node(gene, type=set(["gene"]), score=0, name=gene, source=source)
                    else:
                        kg.nodes[gene]["type"].add("gene")
                        
                    kg.add_edge(gene, pathway_id, type="part_of", score=0, source=source)
                    
    load_kegg_gmt(hallmark_genesets, kg)
    load_kegg_gmt(curated_genesets, kg)
    
    return kg


def load_human_transcription_factors(kg: nx.DiGraph, data_dir, source="human_transcription_factors"):

    tfURL="http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.txt"
    tfFile = os.path.basename(tfURL)

    transcriptionFactorFile = os.path.join(data_dir, tfFile)
    if not transcriptionFactorFile is None and not os.path.exists(transcriptionFactorFile):
        tfDB = pd.read_csv(tfURL, sep="\t")
        tfDB.to_csv(transcriptionFactorFile, sep="\t", index=False)
        
    tfDB = pd.read_csv(transcriptionFactorFile, sep="\t")
    tfs = tfDB[tfDB["Is TF?"] == "Yes"]

    tfgene = set(["gene", "TF"])

    for ri, row in tfs.iterrows():
        
        gene = row["HGNC symbol"]
        
        if not gene in kg.nodes:
            kg.add_node(gene, type=set(tfgene), score=0, name=gene, source=source)
        else:
            kg.nodes[gene]["type"].update(tfgene)
    











