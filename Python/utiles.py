"""
Collection of functions that can be imported by other scripts.
"""
import re, math


def parse_bed(bed_file, columns=[], labels=[], id_subset=set(), BIN_SIZE=1):
    """
    Reads a BED file and rturns a tuple with the follwing data structures:
     
    id2elem dictionary in the format {id:{chr:"chrX", start:123, end:345, id:"id", ...}, id2={...}, ...}
    
    reg2id {chrX:[set(id1,id2),set(),set(id3,id4),set(id4),...], chrY:[set(...),...]}
    
    Parses ';' separated annotation columns, given by the columns and matching labels arguments.
    """
    
    assert len(columns) == len(labels)
    
    id2elem = {}    # mapps each element ID to its features 
    reg2id = {}     # mapps each chromosomal region (indexed bins of given size) to a set of overlapping elements
    
    for line in open(bed_file):
        
        # ignore header or comment lines
        if not line.startswith('#'):
            
            sp = line.strip().split('\t')

            id = sp[3]
            
            # read only those elements with id in the given subset
            if (not id_subset) or id in id_subset : 
                
                chr = sp[0]
                start = int(sp[1])
                end = int(sp[2])
                            
                assert start <= end
                
                # bins on the chromosom for faster access to region
                bins = range(start/BIN_SIZE, end/BIN_SIZE + 1)
                
                # insert values in id2elem
                id2elem[id] = {"id":id, "chr":chr, "start":start, "end":end, "bin":bins}
                
                # read additional annotations
                for i in range(len(columns)):
                    
                    id2elem[id][labels[i]] = sp[columns[i]].split(';') if len(sp) > columns[i] else []
                
                reg2id = update_reg2id(reg2id, id, id2elem[id], BIN_SIZE=BIN_SIZE)

    return id2elem, reg2id

def any_overlap(s1, e1, s2, e2):
    """ Returns True if there is any overlap between reg1 and reg2. """
    return s1 < e2 and s2 < e1

def complete_overlap(s1, e1, s2, e2):
    """Returns True if reg1 completely overlaps reg2, that is if
    reg2 is included in reg1 """
    return s1 <= s2 and e2 <= e1

def complete_included(s1, e1, s2, e2):
    """Returns True if reg1 is completely included in reg2, that is if
    reg2 overlaps reg1 completely """
    return s2 <= s1 and e1 <= e2

def percent50(s1, e1, s2, e2):
    """returns True if reg1 overlaps reg2 to at least 50% """
    
    if not any_overlap(s1, e1, s2, e2):        
        return False

    # compute overlapping bps:
    o = min(e1, e2) - max(s1,s2)

    # get  lenghts of elem2
    m = e2 -s2

    # return elem2 overlap >= 50%
    return m == 0 or  float(o) / m >= 0.5 

def reciprocal50(s1, e1, s2, e2):
    """Returns ture, if the two have a reciprocal overlap of >= 0.5 """

    if not any_overlap(s1, e1, s2, e2):        
        return False

    # compute overlapping bps:
    o = min(e1, e2) - max(s1,s2)

    # get  lenghts of the two elements
    n = e1 - s1
    m = e2 -s2
    
    # return reciprocal overlap >= 50 %
    return float(o)/max(n,m) >= 0.5


def annotate_region(id2elem, id2annot, reg2annot, label, overlap_func=any_overlap, BIN_SIZE=1):
    """Annotates elements in id2elem with overlpping elements in id2annot"""

    # iterate over each element
    for id in id2elem.keys():
        
        elem = id2elem[id]  # dict with features
        chr = elem["chr"]

        annotation_IDs = set()  # set of annotation IDs that overlap the element
                    
        if reg2annot.has_key(chr):
                
            # handle each bin in the region:
            for b in elem["bin"]:
                                    
                # check if for that bin an annotation exists:
                if len(reg2annot[chr]) > b:
                
                    # handle each annotation element in the bin:
                    for aID in reg2annot[chr][b]:
                        
                        # check if there is a real overlap
                        if overlap_func(elem["start"], elem["end"], id2annot[aID]["start"], id2annot[aID]["end"]):
                        
                            annotation_IDs.add(aID)
        
        id2elem[id][label] = annotation_IDs
        
    return id2elem

def get_overlapping_elements(chr, start, end, id2elem, reg2id, overlap_func = any_overlap, BIN_SIZE=1):
    """ returns a set of IDs from reg2id, that overlap at the given region. 
    reg2id = {chrX:[[id1,id2],[],[id3,id4],[id4],...], chrY:[[...],...]}
    """
    overlap_ids = set()
    
    for bin in range(start/BIN_SIZE, (end/BIN_SIZE)+1):
        
        # test if reg2id has ids in that bin:
        if reg2id.has_key(chr) and len(reg2id[chr]) > bin:

            # iterate over each id found:
            for id in reg2id[chr][bin]:
                
                # test for overlap
                if overlap_func(start, end, id2elem[id]["start"], id2elem[id]["end"]):
                    
                    overlap_ids.add(id)
    
    return overlap_ids

def parse_targetgene(inFile, phenotype_terms, pt_genes=set(), BIN_SIZE=1):
    """
    Reads a tab file with genes and stores onley those annotated with a target phenotype if phenotype_terms given else all genes.
    If pt_genes is given it stores only the genes that are contained in the pt_genes set.
    
    returns: 
        
    id2elem dictionary in the format {id:{chr:"chrX", start="123", stop="345", id="id", ...}, id2:{...}, ...}
    
    reg2id {chrX:[set(id1,id2),set(),set(id3,id4),set(id4),...], chrY:[set(...),...]}
    
    """
    #columns=[5]
    #labels=["HPO"]
    
    id2elem = {}
    reg2id = {}
    
    for line in open(inFile):
        
        if not line.startswith('#'):
            
            sp = line.strip().split('\t')
            
            id = sp[3]
            # parse associated HPO terms:
            #phenotype = set(sp[columns[labels.index("HPO")]].split(';'))
            
            # check if only a subset of genes should be read:
            if (not pt_genes) or id in pt_genes : 

                # target phenotype is not given or at last one of them is in the target phenotype set:
                if ( not phenotype_terms ): # or phenotype.intersection(phenotype_terms):
                
                    chr = sp[0]
                    start = int(sp[1])
                    end = int(sp[2])
                    
                    
                    assert start <= end
                    
                    # bins on the chromosom for faster access to region
                    bins = range(start/BIN_SIZE, end/BIN_SIZE + 1)
                    
                    # insert values in id2elem
                    id2elem[id] = {"id":id, "chr":chr, "start":start, "end":end, "bin":bins}
                    
                    reg2id = update_reg2id(reg2id, id, id2elem[id], BIN_SIZE=BIN_SIZE)
                
    return id2elem, reg2id

def update_reg2id(reg2id, id, elem, BIN_SIZE):
    """ updates a reg2ID dict with an additional element."""
    
    # check if chr is already contained
    if elem["chr"] not in reg2id:
        reg2id[elem["chr"]] = []

    # test if bin list should be extended to the largeds value in bins
    while len(reg2id[elem["chr"]]) <= elem["bin"][-1]:
        reg2id[elem["chr"]].append(set())
     
    # append id to every bin set
    for b in elem["bin"]:
        reg2id[elem["chr"]][b].add(id)    

    return reg2id

def parse_term2gene(term_to_gene_file):
    """ returns the mapping from term to genes"""
    
    term2gene ={}
    
    for line in open(term_to_gene_file):
        
        sp = line.strip().split('\t')
    
        if not sp[0] in term2gene:
            term2gene[sp[0]] = set()
            
        term2gene[sp[0]].add(sp[1])
    
    return term2gene
    
def parse_term2gene_UPO(term_to_gene_file):
    """ 
    returns the mapping from term to genes.
    Parses the HPO/crossSpeciesPheno_2_HSgenes.txt file.
    """
    term2gene ={}
    
    for line in open(term_to_gene_file):
        
        sp = line.strip('\n').split(';')
        term = sp[0]
        genes = set([g.split(',')[0] for g in  sp[2:]])

        if not term in term2gene:
            term2gene[term] = set()
            
        term2gene[term].update(genes)
    
    return term2gene
    

def round_to(i, d=3):
    """ rounds i to d significant digits. Returns a string """
    return '%s' % float('%.{0}g'.format(d) % i)

def meanstdv(x): 
    """ returns mean and sd from list x """
    n = len(x)
    mean = sum([i for i in x]) / float(n) if n != 0 else 0.0
    std = math.sqrt(sum([(i - mean)**2 for i in x]) / float(n-1)) if n-1 != 0 else 0.0
    return mean, std
	
