

import os, sys
import urllib
import zipfile
import gzip
import pandas as pd
import numpy as np

from collections import defaultdict
from goatools.anno.gaf_reader import GafReader
from goatools.obo_parser import GODag

from Bio import SeqIO
import re

import networkx as nx
import logging


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
        }, organism="human"):
    
    go_obo_file = "go-basic.obo"
    if not os.path.exists(os.path.join(data_dir, go_obo_file)):
        download_and_unzip("http://geneontology.org/ontology/go-basic.obo", ".", os.path.join(data_dir, go_obo_file), data_dir, unzip=False)

    if "human" in organism:

        goa_gaf_file = "goa_human.gaf"
        if not os.path.exists(os.path.join(data_dir, goa_gaf_file)):
            download_and_unzip("http://geneontology.org/gene-associations/goa_human.gaf.gz", ".", os.path.join(data_dir, goa_gaf_file+".gz"), data_dir)
                

        genenameFile = "hgnc_annot.tsv"
        genenamesURL = 'https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit'

        statusColumn = "Status"
        statusWithdrawn = ["Symbol Withdrawn", "Entry Withdrawn"]
        symbolColumn = "Approved symbol"
        uniprotColumn = "UniProt ID(supplied by UniProt)"


        if not os.path.exists(os.path.join(data_dir, genenameFile)):
            download_and_unzip(genenamesURL, ".", os.path.join(data_dir, genenameFile), data_dir, unzip=False)

    if "mouse" in organism:
        goa_gaf_file = "goa_mouse.gaf"
        if not os.path.exists(os.path.join(data_dir, goa_gaf_file)):
            download_and_unzip("https://current.geneontology.org/annotations/mgi.gaf.gz", ".", os.path.join(data_dir, goa_gaf_file+".gz"), data_dir)
                
        genenameFile = "mgi_annot.tsv"
        genenamesURL = 'https://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt'
        
        statusColumn = "Status"
        statusWithdrawn = []
        symbolColumn = "Marker Symbol"
        uniprotColumn = "UniProt IDs"
        
        if not os.path.exists(os.path.join(data_dir, genenameFile)):
            download_and_unzip(genenamesURL, ".", os.path.join(data_dir, genenameFile), data_dir, unzip=False)



    ogaf = GafReader(os.path.join(data_dir, goa_gaf_file))
    obodag = GODag(os.path.join(data_dir, go_obo_file))

    # read hgnc <-> uniprot conversion
    hgncDF = pd.read_csv(os.path.join(data_dir, genenameFile), sep="\t")
    uniprot2hgnc = defaultdict(set)
    all_genes = set()

    for ri, row in hgncDF.iterrows():
        
        status = row[statusColumn]
        
        if status in statusWithdrawn:
            continue
        
        symbol  = row[symbolColumn]    
        uniprots= row[uniprotColumn]

        all_genes.add(symbol)
        
        if pd.isna(uniprots):
            continue

        for uniprot in uniprots.split("|"):
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
        
        termLevel = termObj.level
        termDepth = termObj.depth
        
        kg.add_node(termID, id=termID, name=termName, type=set(["geneset"]), ns=termNS, score=0, source=source, go_level=termLevel, go_depth=termDepth)
        
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


def load_uniprot_location(kg: nx.DiGraph, data_dir, source="uniprot_celloc", organism="human"):
    
    uniprotDB = os.path.join(data_dir, "uniprot_cellular_location.obo")
    if not os.path.exists(uniprotDB):
        uniprotURL = "https://rest.uniprot.org/locations/stream?format=obo&query=%28*%29"
        download_and_unzip(uniprotURL, data_dir, "uniprot_cellular_location.obo", unzip=False)
        

    if "human" in organism:
        uniprotData = os.path.join(data_dir, "uniprot_cellular_location_human.tsv")
        if not os.path.exists(uniprotData):
            uniprotDataURL = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Cgene_primary%2Ccc_subcellular_location&format=tsv&query=%28%28proteome%3AUP000005640%29%29"
            df = pd.read_csv(uniprotDataURL, sep="\t")
            df.to_csv(uniprotData, sep="\t", index=False)

    if "mouse" in organism:
        uniprotData = os.path.join(data_dir, "uniprot_cellular_location_mouse.tsv")
        if not os.path.exists(uniprotData):
            uniprotDataURL = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Creviewed%2Cid%2Cgene_names%2Corganism_name%2Cgene_primary%2Ccc_subcellular_location&format=tsv&query=(taxonomy_id%3A10090)"
            df = pd.read_csv(uniprotDataURL, sep="\t")
            df.to_csv(uniprotData, sep="\t", index=False)
        
        
    
    uniprotdag = GODag(uniprotDB, optional_attrs=['relationship', "xref"])
    df = pd.read_csv(uniprotData, sep="\t").replace({np.nan: None})

    ccName2ID = {}
    for upEnty in uniprotdag:
        
        termID = upEnty
        termObj = uniprotdag[upEnty]
        
        if termObj.is_obsolete:
            continue
        
        termName = termObj.name
        termNS = str(termObj.namespace)
        
        termLevel = termObj.level
        termDepth = termObj.depth
        
        ccName2ID[termName] = termID

        kg.add_node(termID, id=termID, name=termName, type=set(["geneset", "subcellular_location"]), ns=termNS, score=0, source=source, term_level=termLevel, term_depth=termDepth)
        
        
    for upEnty in uniprotdag:
        
        termID = upEnty
        termObj = uniprotdag[upEnty]
        

        for parID in [x.id for x in termObj.parents]:
            if parID in kg.nodes:
                kg.add_edge(termID, parID, type="relevant_in", score=0, source=source)
            
        for xref in [x for x in termObj.xref]:
            if xref in kg.nodes:
                kg.add_edge(termID, xref, type="relevant_in", score=0, source=source)
                
        for relID in [x.id for x in termObj.relationship.get("part_of", [])]:
            if relID in kg.nodes:
                kg.add_edge(termID, relID, type="part_of", score=0, source=source)
        
            
    def extract_subcellular_annotation(inStr):
        allAnnots = []
        for annotation in inStr.split("SUBCELLULAR LOCATION: "):
            if len(annotation) == 0:
                continue
            
            if annotation.startswith("["):
                annotation = annotation.split("]:")[1]
                
            if "Note=" in annotation:
                annotation = annotation.split("Note=")[0]
            
            annots = [item for sublist in [x.strip().split(" {")[0].split(", ") for x in re.split(r"(;\s)|(\.\s)", annotation) if not x is None and "{" in x] for item in sublist]
            allAnnots.extend(annots)
        return allAnnots  
            

    for ri, row in df.iterrows():
        
        geneSymbols = row["Gene Names (primary)"]
        geneLoc = row["Subcellular location [CC]"]
        
        if geneSymbols is None:
            continue
        
        if geneLoc is None:
            continue
        
        geneSymbols = geneSymbols.split("; ")
        
        #if not geneSymbol in kg.nodes:
        #    continue
        
        geneLocs = extract_subcellular_annotation(geneLoc)
                    
                    
        for geneSymbol in geneSymbols:
            
            if not geneSymbol in kg.nodes:
                kg.add_node(geneSymbol, type=set(["gene"]), source=source, name=geneSymbol, score=0)
            
            for geneLoc in geneLocs:
                
                if not geneLoc in ccName2ID:
                    geneLoc = geneLoc[0].upper() + geneLoc[1:]
                            
                if not geneLoc in ccName2ID:
                    print(geneSymbol, geneLocs)
                
                locID = ccName2ID[geneLoc]
                
                if locID is None:
                    continue
                
                kg.add_edge(geneSymbol, locID, type="part_of", score=0, source=source)
                
    return kg
        


            
def load_omnipath(kg: nx.DiGraph, data_dir, source="omnipath", organism="human"):
    
    #
    ## Checking Files
    #

    if "human" in organism:
        omnipathDB = os.path.join(data_dir, "omnipath_human.tsv")
        geneTranslationDB = os.path.join(data_dir, "omnipath_hgnc_uniprot.tsv")

        genenamesURL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
        
        geneSymColumn = "Approved symbol"
        uniprotColumn = "UniProt ID(supplied by UniProt)"

    if "mouse" in organism:
        omnipathDB = os.path.join(data_dir, "omnipath_mouse.tsv")
        geneTranslationDB = os.path.join(data_dir, "omnipath_mgi_uniprot.tsv")

        genenamesURL = 'https://www.informatics.jax.org/downloads/reports/MRK_Sequence.rpt'

        geneSymColumn = "Marker Symbol"
        uniprotColumn = "UniProt IDs"

    
    if not os.path.exists(geneTranslationDB):

        hgncTranslation = pd.read_csv(genenamesURL, sep="\t")
        hgncTranslation.to_csv(geneTranslationDB, sep="\t")

    if not os.path.exists(os.path.join(data_dir, "miRNA.dat")):
        download_and_unzip("https://mirbase.org/download/miRNA.dat", ".", os.path.join(data_dir, "miRNA.dat"), ".", unzip=False)

    if not os.path.exists(omnipathDB):
        
        import omnipath as op
        opd = op.interactions.AllInteractions().get(organisms=organism, genesymbols=True)
        
        opd.to_csv(omnipathDB, sep="\t")
        
    MI2entry = {}
    MIMAT2sym = {}
    with open(os.path.join(data_dir, "miRNA.dat"), "r") as handle:
        for record in SeqIO.parse(handle, "embl"):

            if record.annotations["data_file_division"] == "HSA":
                #print(record.id, len(record))
        
                mir_id = record.id
                mir_name = record.name
                mir_symbol = "-".join(mir_name.split("-")[1:3])
        
                mature_mirs = set()
                for feature in record.features:
                    if "accession" in feature.qualifiers:
                        for acc in feature.qualifiers["accession"]:
                            mature_mirs.add(acc)
                    #else:
                    #    mature_mirs.add("--NA--")
                
                MI2entry[mir_id] = {"ID": mir_id, "NAME": mir_name, "SYMBOL": mir_symbol, "matures": mature_mirs}

                for acc in mature_mirs:
                    MIMAT2sym[acc] = set([mir_symbol])
    
    
    #
    ## Adding nodes
    # 
    geneTranslation = pd.read_csv(geneTranslationDB, sep="\t")
    opd  = pd.read_csv(omnipathDB, sep="\t")
        
    uniprot2gene = defaultdict(set)

    for ri, row in geneTranslation.iterrows():
        
        hgncGene = row[geneSymColumn]
        uniprotIDs = row[uniprotColumn]
        if pd.isna(uniprotIDs):
            continue
        
        for uniprotID in uniprotIDs.split("|"):
            uniprot2gene[uniprotID].add(hgncGene)
        
    
    interactionTypes = {} # stimulation, inhibition
    interactionTypes[(False, False)] = "interacts"
    interactionTypes[(False, True)] = "represses"
    interactionTypes[(True, False)] = "activates"
    interactionTypes[(True, True)] = "interacts"

    ignoredCount = 0
            
    
    #
    ## Adding edges
    #   
    for ri, row in opd.iterrows():
        
        srcUID = row["source"]
        tgtUID = row["target"]
        
        rowType = row["type"]
        
        if srcUID.startswith("COMPLEX:") or tgtUID.startswith("COMPLEX:"):
            continue
        
        srcType = set(["gene"])
        srcGenes = uniprot2gene.get(srcUID, set())
        tgtType = set(["gene"])
        tgtGenes = uniprot2gene.get(tgtUID, set())
        
        if len(srcGenes) == 0:
            srcType = set(["ncRNA", "miRNA"])
            srcGenes = MIMAT2sym.get(srcUID, set())
        if len(tgtGenes) == 0:
            tgtType = set(["ncRNA", "miRNA"])
            tgtGenes = MIMAT2sym.get(tgtUID, set())
    
        if rowType in ["lncrna_post_transcriptional"]:
            if len(srcGenes) == 0:
                srcType = set(["ncRNA", "lncRNA"])
                srcGenes = set([srcUID])
            if len(tgtGenes) == 0:
                tgtType = set(["ncRNA", "lncRNA"])
                tgtGenes = set([tgtUID])
        
        
        if len(srcGenes) == 0 or len(tgtGenes) == 0:
            continue
            
        for src in srcGenes:
            if src in kg.nodes:
                kg.nodes[src]["type"].update(srcType)
            else:
                kg.add_node(src, type=srcType, source=source, name=src, score=0)
            
        for tgt in tgtGenes:
            if tgt in kg.nodes:
                kg.nodes[tgt]["type"].update(tgtType)
            else:
                kg.add_node(tgt, type=tgtType, source=source, name=tgt, score=0)
            
            
        # consensus_direction -> all resources have same direction
        # consensus_stimulation -> all resources show this as stimulation
        # consensus_inhibition -> all resources show this as consensus_inhibition
        
        # it must have a consensus direction
        if not (row["is_directed"]):
            ignoredCount += 1
            continue
        
        # it must either have consensus stimulation or inhibition
        # if not (row["consensus_stimulation"] or row["consensus_inhibition"]):
        #    ignoredCount += 1
        #    continue

        isStimulation = row["is_stimulation"]
        isInhibition = row["is_inhibition"]

        
        omnipathEvidences = str(row["references"]).strip().split(";")
        omnipathType = row["type"]
        omnipathSources = row["sources"].split(";")
        
        if "miRNA" in srcType and "gene" in tgtType:
            if len(set(omnipathSources).intersection(["ncRDeathDB", "miRDeathDB", "miR2Disease", "miRecords", "miRTarBase"])) > 0:
                if not isStimulation and not isInhibition:
                    isStimulation = False
                    isInhibition = True
        
        interactionType = interactionTypes[ (isStimulation, isInhibition) ]
        
        
        for src in srcGenes:
            
            for tgt in tgtGenes:
                
                attrDict = {
                    "type": interactionType, "source": source, "omnipath_evidences": omnipathEvidences, "omnipath_type": omnipathType
                }
                
                kg.add_edge( src, tgt, **attrDict )
                
                if not row["is_directed"]:
                    kg.add_edge(tgt, src, **attrDict)
            
       
       
                    
    return kg
        
        
    
def load_STRING(kg: nx.DiGraph, data_dir, source="STRING",
                use_evidences = ['fusion', 'coexpression','experiments','database','textmining'], organism="human"):

    if "human" in organism:
        taxid = "9606"
        ensProtIdCol = "Ensembl Protein ID"
        geneSymbolCol = "HGNC symbol"
    if "mouse" in organism:
        taxid = "10090"
        ensProtIdCol = "Ensembl Protein ID"
        geneSymbolCol = "MGI symbol"

    logger = logging.getLogger("load_STRING")       
    stringFile = "{}.protein.links.full.v11.5.txt".format(taxid)

    logger.info("Loading file " + stringFile)

    if not os.path.exists(os.path.join(data_dir, stringFile)):
        #https://stringdb-downloads.org/download/protein.links.full.v11.5/10090.protein.links.full.v11.5.txt.gz
        stringFileURL = "https://stringdb-downloads.org/download/protein.links.full.v11.5/"+stringFile + ".gz"
        logger.info("Downloading file " + stringFileURL)
        download_and_unzip(stringFileURL, ".", os.path.join(data_dir, stringFile+".gz"), data_dir)


    df = pd.read_csv(os.path.join(data_dir, stringFile), sep=" ")
    
    subdf = df[["protein1", "protein2"]+use_evidences]
    subdf["score"] = subdf[use_evidences].max(axis=1)/1000
    subdf = subdf[subdf.score > 0]
    subdf = subdf[subdf.score >= 0.7]

    
    #low confidence - 0.15 (or better),
    #medium confidence - 0.4,
    #high confidence - 0.7,
    #highest confidence - 0.9.

    
    all_ensp_proteins = set()
    allStringProts = set(subdf.protein1)
    allStringProts.update(subdf.protein2)

    for x in allStringProts:
        all_ensp_proteins.add(x.split(".")[1])
    
    mart_file=taxid + ".mart_export_oct2014.txt"
    logger.info("Loading biomart file " + mart_file)
    martDF = pd.read_csv(os.path.join(data_dir, mart_file), sep="\t")

    ensemblProt2Gene = defaultdict(set)
    
    martEmptyProt = ~pd.isna(martDF[ensProtIdCol])
    martEmptyHgnc = ~pd.isna(martDF[geneSymbolCol])
    
    martFilter = np.where(martEmptyProt & martEmptyHgnc)
    
    for ri, row in martDF.loc[martFilter].iterrows():
        
        protid = "{}.{}".format(taxid, row[ensProtIdCol])
        geneid = row[geneSymbolCol]

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

def load_mgi_disease_ontology(kg: nx.DiGraph, data_dir, source="MGI_DO", organism="human"):

    logger = logging.getLogger("load_mgi_disease_ontology")
    logger.info("Adding MGI Disease Ontology")

    mgiDO = pd.read_csv("https://www.informatics.jax.org/downloads/reports/MGI_DO.rpt",sep="\t")

    if organism == "human":
        taxID = 9606

    if organism == "mouse":
        taxID = 10090

    useDO = mouseDO = mgiDO[mgiDO["NCBI Taxon ID"] == taxID]

    for ri, row in useDO.iterrows():

        gene = row["Symbol"]
        disease_id = row["DO Disease ID"]
        disease_label = row["DO Disease Name"]

        if not gene in kg.nodes:
            kg.add_node(gene, type=set(["gene"]), source=source, name=gene, score=0)
        else:
            kg.nodes[gene]["type"].add("gene")

        if not disease_id in kg.nodes:
            kg.add_node(disease_id, type=set(["disease"]), name=disease_label, source=source)
        else:
            kg.nodes[disease_id]["type"].add("disease")

        disease_data = {
            "type": "interacts",
            "source": source
        }
        
        kg.add_edge( gene, disease_id, **disease_data )


def load_opentargets(kg: nx.DiGraph, data_dir, source="opentargets", min_disease_association_score=0.8, organism="human"):
    
    logger = logging.getLogger("load_opentargets")

    if not os.path.exists(os.path.join(data_dir, "opentargets_knowndrugs.tsv")):
        logger.error("Missing OpenTargets data frame knowndrugs. Prepare data according to Jupyter Notebook: opentargets_knowndrugs")
        exit(-1)

    if not os.path.exists(os.path.join(data_dir, "opentargets_disease_associations.tsv")):
        logger.error("Missing OpenTargets data frame knowndrugs. Prepare data according to Jupyter Notebook: opentargets_disease_associations")
        exit(-1)

    human2mouse = defaultdict(set)

    if organism == "mouse":
        # got this idea from: https://www.biostars.org/p/9567892/
        mouse_human_genes = pd.read_csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

        mouse = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 10090]
        human = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 9606]

        mouseDF = mouse.iloc[:, [0, 3]]
        humanDF = human.iloc[:, [0, 3]]
        conversion_df = pd.merge(humanDF, mouseDF, on="DB Class Key", suffixes=["_human", "_mouse"])

        for ri, row in conversion_df.iterrows():
            human2mouse[ row["Symbol_human"] ].add(row["Symbol_mouse"])


    ot_drugs = pd.read_csv("../data/opentargets_knowndrugs.tsv", sep="\t")
    ot_disease = pd.read_csv("../data/opentargets_disease_associations.tsv", sep="\t")
    
    
    #
    ## Adding Disease Nodes
    #
    logger.info("Adding disease nodes")
    if organism == "human":
        for ri, row in ot_disease.iterrows():
            gene = row["targetSymbol"]
            disease_id = row["diseaseId"].replace("_", ":")
            disease_label = row["diseaseLabel"]

            for gene in human2mouse.get(gene, [gene]):
                if not gene in kg.nodes:
                    kg.add_node(gene, type=set(["gene"]), source=source, name=gene, score=0)
                else:
                    kg.nodes[gene]["type"].add("gene")
                
            if not disease_id in kg.nodes:
                kg.add_node(disease_id, type=set(["disease"]), name=disease_label, source=source)
            else:
                kg.nodes[disease_id]["type"].add("disease")
    else:
        logger.info("Skipping adding disease nodes")
            
            
    #
    ## Adding Drug Nodes
    #
    logger.info("Adding Drug Nodes")
    for ri, row in ot_drugs.iterrows():
        drug = row["drugId"]
        disease_id = row["diseaseId"].replace("_", ":")
        diseaseName = row["diseaseName"]
        gene = row["targetGeneSymbol"]
        geneName = row["targetGeneName"]
        
        drugName = row["prefName"]
        drugType = row["drugType"]
        
        for gene in human2mouse.get(gene, [gene]):
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
    logger.info("Adding disease edges")
    if organism == "human":
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
            
            for gene in human2mouse.get(gene, [gene]):
                kg.add_edge( gene, disease_id, **disease_data )
    else:
        logger.info("Skipping adding disease edges")
        
    #
    ## Adding Drug Edges
    #
    logger.info("Adding drug edges")
    for ri, row in ot_drugs.iterrows():
        #drugId	targetId	diseaseId	status	diseaseName	targetGeneSymbol	targetGeneName

            
        drug = row["drugId"]
        disease_id = row["diseaseId"].replace("_", ":")
        drug_target = row["targetGeneSymbol"]
        
        evidence_status = row["status"]
        
        if pd.isna(evidence_status):
            evidence_status = "Unknown status"
        
        if not isinstance(evidence_status, str):
            print("evnostr", drug, disease_id, drug_target, evidence_status)
        
        if "Withdrawn" == evidence_status:
            #do not include edges
            continue
        
        drug_disease_data = {
            "evidence_status": evidence_status,
            "source": source,
            "type": "affected_by"
        }
        
        kg.add_edge( disease_id, drug, **drug_disease_data )

        for dtarget in human2mouse.get(drug_target, [drug_target]):
            kg.add_edge( dtarget, drug, type="target_of", source=source )
        
    return kg
        
        
def load_reactome(kg: nx.DiGraph, data_dir, source="reactome", organism="human"):

    logger = logging.getLogger("load_reactome")

    if organism == "human":
        logger.info("Loading REACTOME human")
        reactomeFile = os.path.join(data_dir,"ReactomePathways.gmt")
        
        if not os.path.exists(reactomeFile):
            logger.info("Downloading REACTOME human")
            download_and_unzip("https://reactome.org/download/current/ReactomePathways.gmt.zip", ".", os.path.join(data_dir,"ReactomePathways.gmt.zip"), data_dir)

        logger.info("Loading REACTOME human from " + reactomeFile)

    if organism == "mouse":
        logger.info("Loading REACTOME mouse")
        reactomeFile = os.path.join(data_dir, "kegg_gmts/mouse/m2.cp.reactome.v2023.2.Mm.symbols.gmt") 
        logger.info("Loading REACTOME mouse from " + reactomeFile)
    
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


def load_mirtarbase(kg: nx.DiGraph, data_dir, source="mirtarbase", organism="human"):

    logger = logging.getLogger("load mirtarbase")
    logger.info("Adding mirtarbase")

    if organism == "human":
        mirtarbaseFile = "{}/hsa_MTI.xls".format(data_dir)
    if organism == "mouse":
        mirtarbaseFile = "{}/mmu_MTI.xls".format(data_dir)


    mirDF = pd.read_excel(mirtarbaseFile)

    for ri, row in mirDF.iterrows():
        mirnaID = row["miRNA"].replace("mmu-", "")
        gene = row["Target Gene"].capitalize()

        if not gene in kg.nodes:
            kg.add_node(gene, type=set(["gene"]), source=source, name=gene, score=0)
        else:
            kg.nodes[gene]["type"].add("gene")

        if not mirnaID in kg.nodes:
            kg.add_node(mirnaID, type=set(["ncRNA", "miRNA"]), source=source, name=mirnaID, score=0)
        else:
            kg.nodes[mirnaID]["type"].add("ncRNA")
            kg.nodes[mirnaID]["type"].add("miRNA")


        mirna_data = {
            "type": "represses",
            "source": source
        }
        
        kg.add_edge( mirnaID, gene, **mirna_data )

    return kg




def load_npinter(kg: nx.DiGraph, data_dir, source="npinter5", organism="human"):

    logger = logging.getLogger("load NPInter")
    logger.info("Adding NPInter")
    logger.warn("NPInter may not be accurate in terms of reported interaction directions.")

    npinterFile = os.path.join(data_dir,"interaction_NPInterv5.txt")
    
    if not os.path.exists(npinterFile):
        download_and_unzip("http://bigdata.ibp.ac.cn/npinter5/download/file/interaction_NPInterv5.txt.gz", ".", os.path.join(data_dir,"interaction_NPInterv5.txt.gz"), data_dir, force_zip=True)


    if organism == "human":
        orgname = "Homo sapiens"
    if organism == "mouse":
        orgname = "Mus musculus"

    df = pd.read_csv(npinterFile, sep="\t")
    df = df[df.organism == orgname].copy()

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
                    'circRNA': ["ncRNA"],
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
    
    
        if ncName.startswith(("hsa-", "mmu-")):
            ncName = ncName[4:]


        if ncID.startswith("NON"):
            # NONHSAG....
            ncName = ncID
        if tarID.startswith("NON"):
            # NONHSAG....
            tarName = tarID
    
    
        kgNcTypes = type2kgtype[ncType]
        kgTarTypes = type2kgtype[tarType]
        
        if kgNcTypes is None or kgTarTypes is None:
            continue
    
        if not ("gene" in kgNcTypes or "gene" in kgTarTypes):
            continue
        
        kgNcTypes=set(kgNcTypes)
        kgTarTypes=set(kgTarTypes)
            
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
              hallmark_genesets="kegg_gmts/human/h.all.v2023.2.Hs.symbols.gmt",
              curated_genesets="kegg_gmts/human/c2.all.v2023.2.Hs.symbols.gmt",
              load_reactome=False, load_wp=True
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
                
                if not load_reactome and pathway_id.startswith("REACTOME_"):
                    continue
                
                if not load_wp and pathway_id.startswith("WP_"):
                    continue
                
                
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
    

def load_TFlink(kg: nx.DiGraph, data_dir, source="TFlink", organism="human"):

    logger = logging.getLogger("load_reactome")

    if organism == "human":
        logger.info("Loading TFLINK "+organism)
        tflinkFile = os.path.join(data_dir,"TFLink_human.gmt")
        
        if not os.path.exists(tflinkFile):
            logger.info("Downloading TFLINK "+organism)
            download_and_unzip("https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_LS_GMT_proteinName_v1.0.gmt", ".",
            os.path.join(data_dir,"TFLink_{}.gmt".format(organism)), data_dir, unzip=False)

        logger.info("Loading TFLINK "+organism+" from " + tflinkFile)
        

    if organism == "mouse":
        logger.info("Loading TFLINK "+organism)
        tflinkFile = os.path.join(data_dir,"TFLink_mouse.gmt")
        
        if not os.path.exists(tflinkFile):
            logger.info("Downloading TFLINK "+organism)
            download_and_unzip("https://cdn.netbiol.org/tflink/download_files/TFLink_Mus_musculus_interactions_LS_GMT_proteinName_v1.0.gmt", ".",
            os.path.join(data_dir,"TFLink_{}.gmt".format(organism)), data_dir, unzip=False)

        logger.info("Loading TFLINK "+organism+" from " + tflinkFile)
        

    with open(tflinkFile) as fin:
    
        for line in fin:
            line = line.strip().split("\t")
            
            tfSymbol = line[0]
            tfUniprot = line[1]
            
            targetGenes = [x for x in line[2:] if x[0].isupper()]
                
            if not tfSymbol in kg.nodes:
                kg.add_node(tfSymbol, id=tfSymbol, type=set(["TF", "gene"]), name=tfSymbol, score=0, source=source)
            else:
                kg.nodes[tfSymbol]["type"].add("TF")
                kg.nodes[tfSymbol]["type"].add("gene")
            
            for gene in targetGenes:
                if not gene in kg.nodes:
                    kg.add_node(gene, type=set(["gene"]), name=gene, source=source, score=0)
                else:
                    kg.nodes[gene]["type"].add("gene")
                    
                kg.add_edge(tfSymbol, gene, type="activates", score=0, source=source)
    
    return kg





def load_xdeathdb(kg: nx.DiGraph, data_dir, source="xdeathdb", organism="human"):

    inputFile = os.path.join(data_dir, "xdeathdb.csv")        
    df = pd.read_csv(inputFile, header=0, encoding="latin-1")
    df = df[~pd.isna(df.cell_death)]
    
    human2mouse = defaultdict(set)

    if organism == "mouse":
        # got this idea from: https://www.biostars.org/p/9567892/
        mouse_human_genes = pd.read_csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

        mouse = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 10090]
        human = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 9606]

        mouseDF = mouse.iloc[:, [0, 3]]
        humanDF = human.iloc[:, [0, 3]]
        conversion_df = pd.merge(humanDF, mouseDF, on="DB Class Key", suffixes=["_human", "_mouse"])

        for ri, row in conversion_df.iterrows():
            human2mouse[ row["Symbol_human"] ].add(row["Symbol_mouse"])

    uniprot_node_names = {}
    for node in kg.nodes:
        if "subcellular_location" in kg.nodes[node].get("type", []):
            uniprot_node_names[ kg.nodes[node].get("name", "") ] = node

    for ri, row in df.iterrows():
        
        cellDeathName = row["cell_death"]
        geneSymbol = row["symbol"]
        location = row["location"]
        disease = row["disease_name"]
                
        if geneSymbol is None:
            continue
        
        if pd.isna(geneSymbol) or pd.isna(location):
            continue
        
        if not pd.isna(disease):
            cellDeathName = "{} ({})".format(cellDeathName, disease)
        
        if not cellDeathName in kg.nodes:
            kg.add_node(cellDeathName, type=set(["geneset", "cell_death_pw"]), score=0, name=cellDeathName, source=source)
        else:
            kg.nodes[cellDeathName]["type"].update(["geneset", "cell_death_pw"])
        
        
        for gene in human2mouse.get(geneSymbol, [geneSymbol]):
            if not gene in kg.nodes:
                kg.add_node(gene, type=set(["gene"]), score=0, name=gene, source=source)
            else:
                kg.nodes[gene]["type"].update(["gene"])
                
            kg.add_edge(gene, cellDeathName, type="relevant_in", score=0, source=source)
            
            if location in uniprot_node_names:
                kg.add_edge(gene, uniprot_node_names[location], type="relevant_in", score=0, source=source)

    return kg

def load_smpdb(kg: nx.DiGraph, data_dir, source="smpdb", organism="human"):

    inputFile = os.path.join(data_dir, "all_smpdb.csv")        
    df = pd.read_csv(inputFile, header=0)
    
    geneTranslationDB = os.path.join(data_dir, "omnipath_hgnc_uniprot.tsv")
    genenamesURL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=md_prot_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
    
    geneSymColumn = "Approved symbol"
    uniprotColumn = "UniProt ID(supplied by UniProt)"

    if not os.path.exists(geneTranslationDB):
        hgncTranslation = pd.read_csv(genenamesURL, sep="\t")
        hgncTranslation.to_csv(geneTranslationDB, sep="\t")

    uniprot2gene = defaultdict(set)
    geneTranslation = pd.read_csv(geneTranslationDB, sep="\t")
    for ri, row in geneTranslation.iterrows():
        
        hgncGene = row[geneSymColumn]
        uniprotIDs = row[uniprotColumn]
        if pd.isna(uniprotIDs):
            continue
        
        for uniprotID in uniprotIDs.split("|"):
            uniprot2gene[uniprotID].add(hgncGene)

    human2mouse = defaultdict(set)

    if organism == "mouse":
        # got this idea from: https://www.biostars.org/p/9567892/
        mouse_human_genes = pd.read_csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

        mouse = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 10090]
        human = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 9606]

        mouseDF = mouse.iloc[:, [0, 3]]
        humanDF = human.iloc[:, [0, 3]]
        conversion_df = pd.merge(humanDF, mouseDF, on="DB Class Key", suffixes=["_human", "_mouse"])

        for ri, row in conversion_df.iterrows():
            human2mouse[ row["Symbol_human"] ].add(row["Symbol_mouse"])


    for ri, row in df.iterrows():
        
        smpdbid = row["SMPDB ID"]
        genes = uniprot2gene.get(row["Uniprot ID"], None)
        pathway = row["Pathway Name"]
        
        if "SMPDB" in smpdbid:
            continue
        if genes is None:
            continue
        
        if not smpdbid in kg.nodes:
            kg.add_node(smpdbid, type=set(["geneset", "metabolic_pathway"]), score=0, name=pathway, source=source)
        else:
            kg.nodes[smpdbid]["type"].update(["geneset", "metabolic_pathway"])
        
        for hgene in genes:    
            for gene in human2mouse.get(hgene, [hgene]):
                if not gene in kg.nodes:
                    kg.add_node(gene, type=set(["gene"]), score=0, name=gene, source=source)
                else:
                    kg.nodes[gene]["type"].update(["gene"])
                    
                kg.add_edge(gene, smpdbid, type="relevant_in", score=0, source=source)
    
    return kg

missingUniprotIDs = set()

def load_metalinks(kg: nx.DiGraph, data_dir, source="metalinks", organism="human"):

    df = pd.read_csv(os.path.join(data_dir, "metalinks_200db_300exp_700pred_900comb.csv"), header=0)
    
    human2mouse = defaultdict(set)

    if organism == "mouse":
        # got this idea from: https://www.biostars.org/p/9567892/
        mouse_human_genes = pd.read_csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

        mouse = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 10090]
        human = mouse_human_genes[mouse_human_genes["NCBI Taxon ID"] == 9606]

        mouseDF = mouse.iloc[:, [0, 3]]
        humanDF = human.iloc[:, [0, 3]]
        conversion_df = pd.merge(humanDF, mouseDF, on="DB Class Key", suffixes=["_human", "_mouse"])

        for ri, row in conversion_df.iterrows():
            human2mouse[ row["Symbol_human"] ].add(row["Symbol_mouse"])


    for ri, row in df.iterrows():
        
        metName = row["MetName"]
        protName = row["Protein"]
        
        
        if metName is None or pd.isna(metName):
            continue
        if protName is None or pd.isna(protName):
            continue
        
       
        combinedScore = row['Combined']
        cellLoc = row["Cellular Location"]
        tissLoc = row["Tissue Location"]
        biospecLoc = row["Biospecimen Location"]
        diseases = row["Diseases"]
        pathways = row["Pathways"]
        
        add_data = {
            "interaction_score": combinedScore,
            "cell_location": cellLoc,
            "tissue_location": tissLoc,
            "biospecimen_location": biospecLoc,
            "diseases": diseases,
            "pathways": pathways,
        }
        
        
        if not metName in kg.nodes:
            kg.add_node(metName, type=set(["metabolite_targets"]), score=0, name=metName, source=source)
        else:
            if kg.nodes[metName].get("source", "") != source:
                metName = metName + "(met)"
                kg.add_node(metName, type=set(["metabolite_targets"]), score=0, name=metName, source=source)
        
        for gene in human2mouse.get(protName, [protName]):
            if not gene in kg.nodes:
                kg.add_node(gene, type=set(["gene"]), score=0, name=gene, source=source)
            else:
                kg.nodes[gene]["type"].update(["gene"])
                
            kg.add_edge(gene, metName, type="target_of", score=0, source=source, **add_data)
    
    return kg
            