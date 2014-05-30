"""
This script runs the barrier CNV analysis.

"""

import argparse # reads commandline arguments
import utiles   # provides some additional functions
import HPO      # provides functions and classes for the phenotype ontology 
import numpy as np # for faster computation with big datasets

def commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cnv_file", type=str, required=True, help="input CNV file in tab separated format. With columns chr, start, end, id")
    parser.add_argument("-d", "--domain_file", type=str, required=True, help="Domain file in .bed format")
    parser.add_argument("-b", "--boundary_file", type=str, required=True, help="Domainboundary file in .bed format")
    parser.add_argument("-g", "--genes_file", type=str, required=True, help="Genes file in .tab format")
    parser.add_argument("-hg", "--term_to_gene_file", type=str, required=True, help="Tab separated file, that maps each phenotype term (including decendants) to its associated genes")
    parser.add_argument("-e", "--enhancer_dir", type=str, help="path to directory with enhancer data matching the target tissues. Assume files with <tissue>.bed.id")
    parser.add_argument("-ef", "--enhancer_file", type=str, help="path to a single file with enhancers.")
    parser.add_argument("-hpo", "--hpo_file", type=str, required=True, help="Human Phenotype Ontology file in .obo format")
    parser.add_argument("-p", "--target_phenotype_file", type=str, required=True, help="Tab separated file with target phenotypes as HPO ID in first column")
    parser.add_argument("-of", "--overlap_function", type=str, choices=["complete", "any", "percent50"], default="complete", help="Function to compute the overlap of boundaries. 'complete' requires the CNV to completely overlap a boundary, 'any' requires only a partial overlap and 'percent50' requires at least 50 percent of the boundary affect.")
    parser.add_argument("-w", "--window_size", type=int, default=400000, help="Window size for testing enhancer adaption mechanism without the boundary disruption effect. That is search for matching enahncer and gene in a fixed distance window in the flanking regions of the deletion.")
    parser.add_argument("-bs", "--bin_size", type=float, default=1000000, help="bin size for faster acces to regions while computiong region overlaps. Default is 10^6")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="output file")
    parser.add_argument("-sf", "--sparse_output_format", action="store_true", help="Write sparse output format, that is only CNVs that match the target_phenotype and only some (see soruce code) columns.")
    return parser.parse_args()

def annotate_barrier_tissue(id2cnv, id2domain, reg2dID, id2target_enhancer, reg2etID, \
    id2target_genes, reg2gtID, id2boundary, window_size, BIN_SIZE) :
    """
    Overlap computation
    Annotates each CNV with lots of labels like "ptEnhancer_ptGene" 
    Now it also searches for distance effect without taking boundary elements into account.
    """
    # heandle each CNV:
    for cID, cnv in id2cnv.items():
        
        # initialize indicator variables and annotation sets:
        ptEnhacer_ptGene = False 
        
        # lists with two sets for upstream and downstream elements
        ptEnhancer = [set(), set()]
        ptGene = [set(), set()]
        
        # coordinates for region of intrest
        region = ["NA", "NA"]   
        
        # test if CNV overlpes completely a boundary
        if cnv["boundary"]:
                     
            # get domain IDs from start end end of CNV
            start_domain = utiles.get_overlapping_elements(cnv["chr"], cnv["start"], cnv["start"], id2domain, reg2dID, utiles.any_overlap, BIN_SIZE=BIN_SIZE)
            end_domain = utiles.get_overlapping_elements(cnv["chr"], cnv["end"], cnv["end"], id2domain, reg2dID, utiles.any_overlap, BIN_SIZE=BIN_SIZE)
            
            # test if start of CNV does not lie in an domain (but in a boundary):
            if not start_domain:

                # get start coordinate of overlapping boundary:
                boundary_start_min = min([id2boundary[bID]["start"] for bID in cnv["boundary"]])
                start_domain = utiles.get_overlapping_elements(cnv["chr"], boundary_start_min-1, boundary_start_min-1, id2domain, reg2dID, utiles.any_overlap, BIN_SIZE=BIN_SIZE)

            if not end_domain:
                
                # get end coordinated of overlapping boundary:
                boundary_end_max = max([id2boundary[bID]["end"] for bID in cnv["boundary"]])
                end_domain = utiles.get_overlapping_elements(cnv["chr"], boundary_end_max+1, boundary_end_max+1, id2domain, reg2dID, utiles.any_overlap, BIN_SIZE=BIN_SIZE)
            
            # test if both CNV border lies inside domains (and not in boundaries or unorganized chromatin)
            if start_domain and end_domain and len(start_domain) == len(end_domain) == 1:
                
                start_dID = start_domain.pop()
                end_dID = end_domain.pop()
                
                assert start_dID != end_dID
                
                # get intervales between CNV and next domain boundaries:
                reg_start = [id2domain[start_dID]["start"], cnv["start"]]
                reg_end = [cnv["end"], id2domain[end_dID]["end"]]

                    
                # get enhancers and target genes in domains adjacent to CNV, that are not overlapped by the CNV:
                ptEnhancer[0] = utiles.get_overlapping_elements(cnv["chr"], reg_start[0], reg_start[1], id2all_enhancer, reg2eaID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
                ptEnhancer[1] = utiles.get_overlapping_elements(cnv["chr"], reg_end[0], reg_end[1], id2all_enhancer, reg2eaID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)                

                ptGene[0] = utiles.get_overlapping_elements(cnv["chr"], reg_start[0], reg_start[1], id2target_genes, reg2gtID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
                ptGene[1] = utiles.get_overlapping_elements(cnv["chr"], reg_end[0], reg_end[1], id2target_genes, reg2gtID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
                
                # save region of interest for visualization 
                region = [reg_start[0], reg_end[1]]
            
            
        # save annotation for each CNV:
        id2cnv[cID]["ptEnhancer"] = ptEnhancer
        id2cnv[cID]["ptGene"] = ptGene
        
        # check candidate events by testing if annotation sets are not empty for bothe, forwared and reverse interaction direction
        id2cnv[cID]["ptEnhancer_ptGene"] = bool( (ptEnhancer[0] and ptGene[1]) or (ptEnhancer[1] and ptGene[0]) )
        # save region 
        id2cnv[cID]["region"] = region
        
        ################################################################
        # Analyse distance effect: ignore boundaries and search for Gene 
        # and enhancer in fixed window around CNV.
        # This event is here named 'enhancer_adaption'.
        # variables and identifier are shortened by 'dist_' 
        ################################################################
        
        # lists with two sets for upstream and downstream elements
        dist_ptEnhancer = [set(), set()]
        dist_ptGene = [set(), set()]
        
        # get intervales between CNV and next domain boundaries:
        dist_reg_start = [cnv["start"] - window_size, cnv["start"]]
        dist_reg_end = [cnv["end"], cnv["end"] + window_size]

        # get enhancers and target genes in domains adjacent to CNV, that are not overlapped by the CNV:
        dist_ptEnhancer[0] = utiles.get_overlapping_elements(cnv["chr"], dist_reg_start[0], dist_reg_start[1], id2all_enhancer, reg2eaID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
        dist_ptEnhancer[1] = utiles.get_overlapping_elements(cnv["chr"], dist_reg_end[0], dist_reg_end[1], id2all_enhancer, reg2eaID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)                

        dist_ptGene[0] = utiles.get_overlapping_elements(cnv["chr"], dist_reg_start[0], dist_reg_start[1], id2target_genes, reg2gtID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
        dist_ptGene[1] = utiles.get_overlapping_elements(cnv["chr"], dist_reg_end[0], dist_reg_end[1], id2target_genes, reg2gtID, utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
        
        # save adjacent enahncers and genes
        id2cnv[cID]["dist_ptEnhancer"] = dist_ptEnhancer
        id2cnv[cID]["dist_ptGene"] = dist_ptGene
        
        # check candidate events by testing if annotation sets are not empty for bothe, forwared and reverse interaction direction
        id2cnv[cID]["dist_ptEnhancer_ptGene"] = bool( (dist_ptEnhancer[0] and dist_ptGene[1]) or (dist_ptEnhancer[1] and dist_ptGene[0]) )

        # save region of interest for visualization 
        id2cnv[cID]["dist_region"] = [dist_reg_start[0], dist_reg_end[1]]
        
    return id2cnv


def write_analysed_cnv(id2cnv, output_file, tisue_name, phenotype_terms, id2target_genes):
    """ writes a tab separated output file with annotations """
    
    with open(output_file, 'w') as outHandle:
        
        header = ["chr", "start", "end", "name", "type", \
        "HPOterms", "target_term", \
        "overlap_score", "adjacent_score", "dist_score", \
        "target_phenotype", "overlap_PT_genes", \
        "boundary", "nr_boundary", "ptEnhancer_ptGene", \
        "tissue_name", "ptEnhancer_up", "ptEnhancer_down", \
        "ptGene_up", "ptGene_down", "common_HPO_gene_up", "common_HPO_gene_down", \
        "reg_start", "reg_end", \
        "dist_ptEnhancer_ptGene", "dist_ptEnhancer_up", "dist_ptEnhancer_down","dist_ptGene_up", "dist_ptGene_down", \
        "dist_common_HPO_gene_up", "dist_common_HPO_gene_down", \
        "dist_reg_start", "dist_reg_end"
        ]
        outHandle.write("\t".join(header) + '\n')
        
        for cID, cnv in id2cnv.items():

            #~ phenotype = "phenotype" if cnv["target_phenotype"] else "not_phenotype"
            phenotype = bool(cnv["target_phenotype"])
            #~ overlap_PT_genes = "overlap_PT_genes" if cnv["overlap_target_genes"] else "not_overlap_PT_genes"
            overlap_PT_genes = bool(cnv["overlap_target_genes"])

            reg_start = cnv["region"][0]
            reg_end = cnv["region"][1]
            
            # output helper function
            def join_set(s):
                return ";".join(s) if s else "NA"
            
            # build set of common HPO terms between target term and adjacent genes:
            common_terms_gene_up = set()
            for gID in cnv["ptGene"][0]:
                common_terms_gene_up.update( phenotype_terms.intersection(id2target_genes[gID]["HPO"]) )
            
            common_terms_gene_down = set()
            for gID in cnv["ptGene"][1]:
                common_terms_gene_down.update( phenotype_terms.intersection(id2target_genes[gID]["HPO"]) )
            
            # calculate the same for window distance genes
            dist_common_terms_gene_up = set()
            for gID in cnv["dist_ptGene"][0]:
                common_terms_gene_up.update( phenotype_terms.intersection(id2target_genes[gID]["HPO"]) )
            
            dist_common_terms_gene_down = set()
            for gID in cnv["dist_ptGene"][1]:
                common_terms_gene_down.update( phenotype_terms.intersection(id2target_genes[gID]["HPO"]) )
            
            outCols = [cnv["chr"], cnv["start"], cnv["end"], cID, \
                cnv["type"][0], \
                join_set(cnv["HPO"]), \
                join_set(cnv["target_term"]), \
                cnv["overlap_score"][0], cnv["adjacent_score"][0], cnv["dist_score"][0], \
                phenotype, \
                overlap_PT_genes, \
                bool(cnv["boundary"]), \
                len(cnv["boundary"]), \
                id2cnv[cID]["ptEnhancer_ptGene"], \
                tissue_name, \
                join_set(id2cnv[cID]["ptEnhancer"][0]), \
                join_set(id2cnv[cID]["ptEnhancer"][1]), \
                join_set(id2cnv[cID]["ptGene"][0]), \
                join_set(id2cnv[cID]["ptGene"][1]), \
                join_set(common_terms_gene_up), \
                join_set(common_terms_gene_down), \
                reg_start, reg_end, \
                id2cnv[cID]["dist_ptEnhancer_ptGene"], \
                join_set(id2cnv[cID]["dist_ptEnhancer"][0]), join_set(id2cnv[cID]["dist_ptEnhancer"][1]), \
                join_set(id2cnv[cID]["dist_ptGene"][0]), join_set(id2cnv[cID]["dist_ptGene"][1]), \
                join_set(common_terms_gene_up), join_set(common_terms_gene_down), \
                cnv["dist_region"][0], cnv["dist_region"][1] \
                ]
            
            #        "dist_ptEnhancer_ptGene", "dist_ptEnhancer_up", "dist_ptEnhancer_down","dist_ptGene_up", "dist_ptGene_down", \
            # "dist_reg_start", "dist_reg_end"
            assert len(outCols) == len(header)
            
            outHandle.write("\t".join([str(c) for c in outCols]) + "\n")

def get_rates_for_target(id2cnv, rand_n, nr_bins=6):
    """
    Compute hitrates and number of subset terms and store it in a dict
    """
    # build mapping of randomization number to cIDs : 
    r2cID = []
    if rand_n > 1:    
        for i in range(rand_n):
            r2cID.append([])
        for cID in id2cnv.keys():
            r2cID[int(cID.split('_')[2]) - 1 ].append(cID)    
    else:
        r2cID.append([cID for cID in id2cnv.keys() if  bool(id2cnv[cID]["target_phenotype"])])

    # initialize empty list for each randomization:
    values = {}
    values["N"] = np.zeros(rand_n, dtype=np.int32)
    values["boundary"] = np.zeros(rand_n, dtype=np.int32)
    values["TDBD_and_Mixed"] = np.zeros(rand_n, dtype=np.int32)
    values["TDBD"] = np.zeros(rand_n, dtype=np.int32)
    values["GDE"] = np.zeros(rand_n, dtype=np.int32)
    values["no_data"] = np.zeros(rand_n, dtype=np.int32)
    values["EA_and_Mixed"] = np.zeros(rand_n, dtype=np.int32)
    values["EA"] = np.zeros(rand_n, dtype=np.int32)
    values["EA_no_boundary"] = np.zeros(rand_n, dtype=np.int32)
    for bin in range(1, nr_bins + 1):
        values["N_bin_" + str(bin)] = np.zeros(rand_n, dtype=np.int32)
        values["TDBD_and_Mixed_bin_" + str(bin)] = np.zeros(rand_n, dtype=np.int32)
        values["TDBD_bin_" + str(bin)] = np.zeros(rand_n, dtype=np.int32)
    
    for i in range(rand_n):

        values["N"][i] = len(r2cID[i])        

        # compute rates
        for cID in r2cID[i]:
            
            ol_score = float(id2cnv[cID]["overlap_score"][0])
            adj_score = float(id2cnv[cID]["adjacent_score"][0])
            dist_score = float(id2cnv[cID]["dist_score"][0])

            values["boundary"][i] +=  bool(id2cnv[cID]["boundary"]) 
            values["TDBD_and_Mixed"][i] += id2cnv[cID]["ptEnhancer_ptGene"] 
            values["TDBD"][i] += id2cnv[cID]["ptEnhancer_ptGene"] and (adj_score > ol_score )
            values["GDE"][i] += (not id2cnv[cID]["ptEnhancer_ptGene"]) and  (ol_score > 0)
            values["no_data"][i] += (not id2cnv[cID]["ptEnhancer_ptGene"]) and  (ol_score <= 0)
            values["EA_and_Mixed"][i] += id2cnv[cID]["dist_ptEnhancer_ptGene"]
            values["EA"][i] += id2cnv[cID]["dist_ptEnhancer_ptGene"] and (dist_score > ol_score )
            values["EA_no_boundary"][i] += id2cnv[cID]["dist_ptEnhancer_ptGene"] and (dist_score > ol_score) and not bool(id2cnv[cID]["boundary"])
            
            # get bin (nr. of boundary overlaped)
            bin = len(id2cnv[cID]["boundary"])
            if bin > 0:
                if bin >= nr_bins:
                    bin = nr_bins
                
                values["N_bin_" + str(bin)][i] += 1
                values["TDBD_and_Mixed_bin_" + str(bin)][i] += id2cnv[cID]["ptEnhancer_ptGene"]
                values["TDBD_bin_" + str(bin)][i] += id2cnv[cID]["ptEnhancer_ptGene"] and (adj_score > ol_score ) 

    return values
    
def wirte_rates(rates, outFile, nr_bins=6):
    """ takes a list with values dictionaries from  get_rates_for_target for all target phenotypes and writes them to 
    an output file so that each line is one randomization
    """
    terms_n = len(rates)
    rand_n = len(rates[0]["N"])
    
    with open(outFile, 'w') as outHandle:
        
        header = ["N", "TDBD_and_Mixed", "TDBD", "GDE", "no_data", "boundary", "EA_and_Mixed", "EA", "EA_no_boundary"]
        for bin in range(1, nr_bins + 1):
            header.append("N_bin_" + str(bin))
            header.append("TDBD_and_Mixed_bin_" + str(bin))
            header.append("TDBD_bin_" + str(bin))
        outHandle.write("\t".join(header) + "\n")
        
        # iterate over each randomization:
        for j in range(rand_n):
            
            Nj =  sum([tT["N"][j] for tT in rates])
            # compute rates:
            cols = []
            cols.append( Nj )
            cols.append( sum([tT["TDBD_and_Mixed"][j] for tT in rates]) )
            cols.append( sum([tT["TDBD"][j] for tT in rates]) )
            cols.append( sum([tT["GDE"][j] for tT in rates]) )
            cols.append( sum([tT["no_data"][j] for tT in rates]) )
            cols.append( sum([tT["boundary"][j] for tT in rates]) )
            cols.append( sum([tT["EA_and_Mixed"][j] for tT in rates]) )
            cols.append( sum([tT["EA"][j] for tT in rates]) )
            cols.append( sum([tT["EA_no_boundary"][j] for tT in rates]) )

            for bin in range(1, nr_bins + 1):
                Nb = sum([tT["N_bin_" + str(bin)][j] for tT in rates])
                cols.append( Nb )
                cols.append( sum([tT["TDBD_and_Mixed_bin_" + str(bin)][j] for tT in rates]) )
                cols.append( sum([tT["TDBD_bin_" + str(bin)][j] for tT in rates]) )
            
            outHandle.write("\t".join([str(c) for c in cols]) + "\n")
    
if __name__ == "__main__":

    # read commandline argumets
    args = commandline()
    
    BIN_SIZE = int(args.bin_size)
    
    print "INFO: Begin to read HPO data..."
    ontology = HPO.Ontology(args.hpo_file)
    print "INFO: Finished parsing HPO ."
    
    print "INFO: Begin to read CNV data from file '{0}'...".format(args.cnv_file)
    plain_id2cnv, plain_reg2cID = utiles.parse_bed(args.cnv_file, \
        columns=[4, 5, 6, 7, 8, 9], \
        labels=["type", "HPO", "overlap_score", "adjacent_score", "dist_score", "target_term"], BIN_SIZE=BIN_SIZE)
    
    print "INFO: Begin to read boundaries data..."
    id2boundary, reg2bID = utiles.parse_bed(args.boundary_file, BIN_SIZE=BIN_SIZE)
    print "INFO: Finished reading boundaries."
    
    print "INFO: Begin to read domain data..."
    id2domain, reg2dID = utiles.parse_bed(args.domain_file, BIN_SIZE=BIN_SIZE)
    print "INFO: Finished reading domains."
    
    print "INFO: Read term to gene mapping... "
    term2gene = utiles.parse_term2gene(args.term_to_gene_file)
    print "INFO: Finished reading of term to gene mapping."
    
    ####################################################################
    # Annotate some regions...
    ####################################################################        
    if args.overlap_function == "complete": overlap_func = utiles.complete_overlap
    elif args.overlap_function == "any": overlap_func = utiles.any_overlap
    elif args.overlap_function == "percent50": overlap_func = utiles.percent50
    
    print "INFO: Beginn to annotate CNVs with boundaries..."
    plain_id2cnv = utiles.annotate_region(plain_id2cnv, id2boundary, reg2bID, "boundary", overlap_func=utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
    #~ plain_id2cnv = utiles.annotate_region(plain_id2cnv, id2boundary, reg2bID, "boundary", overlap_func=utiles.percent50, BIN_SIZE=BIN_SIZE)
    #~ plain_id2cnv = utiles.annotate_region(plain_id2cnv, id2boundary, reg2bID, "boundary", overlap_func=utiles.any_overlap, BIN_SIZE=BIN_SIZE)
    print "INFO: Finished annotating CNVs with boundaries."
    
    if args.sparse_output_format:
        # list of for rates of all target phenotypes:
        rates = []
        # get number of random permutation if input is random data
        random_data = len(plain_id2cnv.keys()[0].split('_')) == 3
        rand_n = max([int(cID.split('_')[2]) for cID in plain_id2cnv.keys()]) if random_data else 1
    ####################################################################
    # for each target term run the whole analysis:
    ####################################################################
    for line in open(args.target_phenotype_file):
        
        fields = line.strip().split('\t')
        phenotype = fields[0]
        phenotype_name = fields[1]
        tissue_name = fields[2]
        
        print "#=======================================================#"
        print "INFO: Run analysis for term: {0} := '{1}'. Tissue: {2}".format(phenotype, phenotype_name, tissue_name)
    
        # get alternative term IDs for the target phenotype
        alt_terms = set([phenotype])
        if phenotype in ontology.id2alt:
            alt_terms.update(ontology.id2alt[phenotype])
        if phenotype in ontology.alt2id:
            alt_terms.update(ontology.alt2id[phenotype])

        # get set of all phenotypes that are descendants of the target term:
        phenotype_terms = set()
        for tID in alt_terms:
            if tID in ontology.id2descendants:
                phenotype_terms.update(ontology.id2descendants[tID])
        
        # set of all genes associated with the target phenotype
        target_genes_subset = set()
        for tID in alt_terms:
            if tID in term2gene:
                target_genes_subset.update(term2gene[tID])
                
        #id2target_genes, reg2gtID = utiles.parse_targetgene(args.genes_file, set(), pt_genes=term2gene[phenotype], BIN_SIZE=BIN_SIZE)
        #id2target_genes, reg2gtID = utiles.parse_bed(args.genes_file, id_subset=term2gene[phenotype], BIN_SIZE=BIN_SIZE)
        id2target_genes, reg2gtID = utiles.parse_bed(args.genes_file, columns=[5], labels=["HPO"], id_subset=target_genes_subset, BIN_SIZE=BIN_SIZE)
            
        # get tissue specific enhancers:
        enhancer_file = args.enhancer_dir + '/' + tissue_name + ".tab" if not args.enhancer_file else args.enhancer_file
        id2all_enhancer, reg2eaID = utiles.parse_bed(enhancer_file, BIN_SIZE=BIN_SIZE)
        
        # remove annotations by copying the original plain cnv dicts
        #~ id2cnv = plain_id2cnv
        #~ reg2cID = plain_reg2cID
        
        # get subset of CNVs that correspond to the target phenotype:
        print "INFO: test CNVs for matching phenotype terms..."
        if args.sparse_output_format:
            id2cnv = {cID:cnv for cID, cnv in plain_id2cnv.items() if phenotype == plain_id2cnv[cID]["target_term"][0] }
        else:
            id2cnv = plain_id2cnv
        for cID in id2cnv:
            id2cnv[cID]["target_phenotype"] = phenotype == plain_id2cnv[cID]["target_term"][0]

        # annotate CNVs with overlapping target genes:
        print "INFO: annote CNVs wiht overlapping target Genes."
        id2cnv = utiles.annotate_region(id2cnv, id2target_genes, reg2gtID, "overlap_target_genes", overlap_func=utiles.complete_overlap, BIN_SIZE=BIN_SIZE)
        
        #~ if args.sparse_output_format:
            # take only a subset of CNVs that have the 
        ####################################################################    
        # annotate each CNV with the information if it hits a boundary with adjacent enhancers and genes that match the target phenotype
        ####################################################################
        print "INFO: Compute barrier hits..."
        id2cnv = annotate_barrier_tissue(\
            id2cnv, id2domain, reg2dID, \
            id2all_enhancer, reg2eaID, \
            id2target_genes, reg2gtID, \
            id2boundary, args.window_size, BIN_SIZE)        
        ####################################################################
        # write CNVs to output file
        ####################################################################
        if args.sparse_output_format:
            print "INFO: Calculate hit rates..."
            rates.append(get_rates_for_target(id2cnv, rand_n))
            print "INFO: Finished hit rates."

        else:
            print "INFO: write output CNV file ..."
            write_analysed_cnv(id2cnv, args.output_file + "." + phenotype, tissue_name, phenotype_terms, id2target_genes)
            ####################################################################
            # get values for .bg.tab file
            ####################################################################
            id2del = {cID: cnv for (cID, cnv) in id2cnv.items() if "loss" in cnv["type"] or "gain" in cnv["type"]} 
            id2delAnyPT = {cID: cnv for (cID, cnv) in id2del.items() if cnv["HPO"] != ["NA"]} 
            id2delTargetPT = {cID: cnv for (cID, cnv) in id2delAnyPT.items() if phenotype_terms.intersection(cnv["target_term"]) } 
            id2delTargetPTnotPTgene = {cID: cnv for (cID, cnv) in id2delTargetPT.items() if cnv["overlap_target_genes"] }
            
            variant_dicts = [id2cnv, id2del, id2delAnyPT, id2delTargetPT, id2delTargetPTnotPTgene]
            variant_names = ["all", "del", "del any PT", "del target PT", "del target PT not PTgene"]
            
            with open(args.output_file + "." + phenotype + ".bg.tab", 'w') as outHandle:
                
                header = ["variant_type", "n cases", "avg length (Mb)", "avg HPO terms", "avg boundaries", "avg target genes", "avg target enhancers"]
                outHandle.write("\t".join(header) + '\n')
                
                # for all types of variants (all cnv, del, del with PT...)
                for id2var, var_name in zip(variant_dicts, variant_names):
        
                    N = len(id2var)
                    
                    # sieze in mega base pairs (Mb)
                    sizes = [(var["end"] - var["start"]) / float(10**6) for var in id2var.values()]
                    avg_size, sd_size = utiles.meanstdv(sizes)
                            
                    out_cols = [var_name, str(N), "{0}($\pm${1} SD)".format(utiles.round_to(avg_size, 3), utiles.round_to(sd_size, 3)) ]
                    
                    for annot in ["HPO", "boundary", "overlap_target_genes"]: 
                        avg_annot, sd_annot = utiles.meanstdv([len(var[annot]) for var in id2var.values()])
                        out_cols.append("{0}($\pm${1} SD)".format(utiles.round_to(avg_annot, 3), utiles.round_to(sd_annot, 3)))
                
                    outHandle.write("\t".join(out_cols) + '\n')
        
            ####################################################################
            # get some count values for .stats file
            ####################################################################
            N = len(id2delAnyPT)
            N_boundary = len(id2boundary)   # number of boundaries in the genome
            N_domains = len(id2domain)   # number of domains in the genome
            mean_boundary_size = utiles.round_to(sum([b["end"] - b["start"] for b in id2boundary.values()]) / float(N_boundary), 3)
            mean_domain_size = utiles.round_to(sum([d["end"] - d["start"] for d in id2domain.values()]) / float(N_domains), 3)
            
            boundary_cIDs = set([cID for cID in id2delAnyPT if id2delAnyPT[cID]["boundary"]])    
            ptEnhancer_ptGene = set([cID for cID in id2delAnyPT if id2delAnyPT[cID]["ptEnhancer_ptGene"]])
            
            PT_cIDs = set([cID for cID in id2delAnyPT if id2delAnyPT[cID]["target_phenotype"] ])
            
            ####################################################################
            # wirte some stats to output.stats files
            ####################################################################
            with open(args.output_file + "." + phenotype + ".stats", 'w') as outHandle:
                
                outHandle.write("Input HPO term\t{0}\n".format(phenotype))
                outHandle.write("Target HPO terms\t{0}\n".format(len(phenotype_terms)))
                outHandle.write("Background enhancers\t{0}\n".format(len(id2all_enhancer)))
                outHandle.write("Target PT genes\t{0}\n".format(len(id2target_genes)))
        
                outHandle.write("Domains\t{0}\n".format(len(id2domain)))
                outHandle.write("Mean domain size\t{0}\n".format(mean_domain_size))     
                outHandle.write("Boundaries\t{0}\n".format(N_boundary))        
                outHandle.write("Mean boundary size\t{0}\n".format(mean_boundary_size))        
                
                outHandle.write("Deletions\t{0}\n".format(N))
                
                outHandle.write("\quad Boundary\t{0}({1})\n".format(len(boundary_cIDs), utiles.round_to(len(boundary_cIDs) / float(N) if N!=0 else 0.0, 3) ))        
                outHandle.write("\quad ptEnhancer_ptGene\t{0}({1})\n".format(len(ptEnhancer_ptGene), utiles.round_to(len(ptEnhancer_ptGene) / float(N) if N!=0 else 0.0, 3) ))
                
                outHandle.write("\quad Deletions with target phenotype:\t{0}\n".format(len(PT_cIDs)))
                
                outHandle.write("\quad \quad Boundary\t{0}({1})\n".format(len(boundary_cIDs.intersection(PT_cIDs)), utiles.round_to(len(boundary_cIDs.intersection(PT_cIDs)) /  float(len(PT_cIDs)) if len(PT_cIDs) > 0 else 0.0, 3) ))        
                outHandle.write("\quad \quad ptEnhancer_ptGene\t{0}({1})\n".format(len(ptEnhancer_ptGene.intersection(PT_cIDs)), utiles.round_to(len(ptEnhancer_ptGene.intersection(PT_cIDs)) / float(len(PT_cIDs)) if len(PT_cIDs) > 0 else 0.0, 3) ))
                
    
    # put all output files together:
    if not args.sparse_output_format:    
        with open(args.output_file, 'w') as outHandle:
            for pt_nr, line in enumerate(open(args.target_phenotype_file)):
                fields = line.strip().split('\t')
                phenotype = fields[0]
                for i, inLine in enumerate(open(args.output_file + "." + phenotype )):
                    if i>0 or pt_nr == 0:
                        outHandle.write(inLine)
    else:
        # write hit rates
        wirte_rates(rates, args.output_file)

    print "INFO: END OF RUN"
