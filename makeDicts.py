import cPickle as pic
import collections

def loadGo():
    """
    Columns are:                 Contents:

     1) DB                      - database contributing the file (always "SGD" for this file)
     2) DB_Object_ID                - SGDID
     3) DB_Object_Symbol                - see below
     4) NOT             (optional)  - 'NOT', 'contributes_to', or 'colocalizes_with' qualifier for a GO annotation, when needed
     5) GO ID                   - unique numeric identifier for the GO term
     6) DB:Reference(|DB:Reference)         - the reference associated with the GO annotation
     7) Evidence                    - the evidence code for the GO annotation
     8) With (or) From      (optional)  - any With or From qualifier for the GO annotation
     9) Aspect                  - which ontology the GO term belongs in
    10) DB_Object_Name(|Name)   (optional)  - a name for the gene product in words, e.g. 'acid phosphatase'
    11) DB_Object_Synonym(|Synonym) (optional)  - see below
    12) DB_Object_Type              - type of object annotated, e.g. gene, protein, etc.
    13) taxon(|taxon)               - taxonomic identifier of species encoding gene product
    14) Date                    - date GO annotation was made
    15) Assigned_by                 - source of the annotation (e.g. SGD, UniProtKB, YeastFunc, bioPIXIE_MEFIT)
    for item in goDict[(lineList[10].split('|'))[0]]:
                if item == lineList[4]: 
                    inGeneToGO == True
            if inGeneToGO == False:
                goDict[(lineList[10].split('|'))[0]].append(lineList[4])
            inGeneToGO = False
        print goDict
    """
    gene_to_go_dict = collections.defaultdict(list)
    SGDID_to_go_dict = collections.defaultdict(list)
    SGDID_All = collections.defaultdict(list)
    SGDID_to_gene = {}
    gene_to_SGDID = {}
    with open("gene_association.sgd", "r") as goFile:
        for i in xrange(0,8):
            goFile.readline()
        for lines in goFile:
            lineList = lines.split('\t')
            if not SGDID_to_go_dict[lineList[1]]:
                gene_to_go_dict[(lineList[10].split('|'))[0]] = []
                SGDID_to_go_dict[lineList[1]] = []
                SGDID_All[lineList[1]] = []
            gene_to_go_dict[(lineList[10].split('|'))[0]].append(lineList[4])
            SGDID_to_go_dict[lineList[1]].append(lineList[4])
            SGDID_All[lineList[1]].append(lineList[0:])
            SGDID_to_gene[lineList[1]] = (lineList[10].split('|'))[0]
            gene_to_SGDID[(lineList[10].split('|'))[0]] = lineList[1]
    pic.dump(gene_to_go_dict, open( "gene_to_go.pkl", "wb" ))
    pic.dump(SGDID_to_go_dict, open( "SGDID_to_go.pkl", "wb" ))
    pic.dump(SGDID_All, open( "SGDID_All.pkl", "wb" ))
    pic.dump(SGDID_to_gene, open( "SGDID_to_gene.pkl", "wb" ))
    pic.dump(gene_to_SGDID, open( "gene_to_SGDID.pkl", "wb" ))


def loadLit():
    """
    gene_literature.tab :

    Columns are :           Contents:

    1) PubMed ID (optional)     - the unique PubMed identifer for a reference
    >2) citation (mandatory)        - the citation for the publication, as stored in SGD
    3) gene name (optional)     - Gene name, if one exists
    4) feature (optional)       - Systematic name, if one exists
    5) literature_topic (mandatory) - all associated Literature Topics of the SGD Literature Guide
                      relevant to this gene/feature within this paper
                      Multiple literature topics are separated by a '|' character.
    >6) SGDID (mandatory)       - the SGDID, unique database identifier, for the gene/feature
    """
    SGDID_to_lit = collections.defaultdict(list)
    with open("gene_literature.tab.txt", "r") as litFile:
        for lines in litFile:
            lineList = lines.split('\t')
            if not SGDID_to_lit[lineList[5]]:
                SGDID_to_lit[lineList[5]] = []
            SGDID_to_lit[lineList[5]] = lineList[1]
    pic.dump(SGDID_to_lit, open( "SGDID_to_lit.pkl", "wb" ))


def loadSlimTerm():
    """
    go_slim_mapping.tab This file is TAB delimited and contains the mapping of all yeast gene products (protein or RNA)
    to a GO-Slim term.

    Columns:    

    1) ORF (mandatory)      - Systematic name of the gene
    2) Gene (optional)      - Gene name, if one exists
    >3) SGDID (mandatory)       - the SGDID, unique database identifier for the gene
    4) GO_Aspect (mandatory)    - which ontology: P=Process, F=Function, C=Component
    >5) GO Slim term (mandatory)    - the name of the GO term that was selected as a GO Slim term
    6) GOID (optional)      - the unique numerical identifier of the GO term
    7) Feature type (mandatory)     - a description of the sequence feature, such as ORF or tRNA
    """
    SGDID_to_goSlim = collections.defaultdict(list)
    with open("go_slim_mapping.tab.txt", "r") as slimFile:
        for lines in slimFile:
            lineList = lines.split('\t')
            if not SGDID_to_goSlim[lineList[2]]:
                SGDID_to_goSlim[lineList[2]] = []
            SGDID_to_goSlim[lineList[2]].append(lineList[4])
    pic.dump(SGDID_to_goSlim, open( "SGDID_to_goSlim.pkl", "wb" ))
    print SGDID_to_goSlim

loadGo()
loadLit()
loadSlimTerm()
