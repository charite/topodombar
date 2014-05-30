"""
Reads phenomatch scores CnvStatistics.jar output and computes phenogram scores for each CNV.
Replaces the thre columns with phenomatch scores with the following three columns:
1. col: phenogram score for genes overlaped by the CNV
2. col: maximum of phenogram scores for upstrem or downstrem adjecent genes.
3. col: maximum of phenogram scores distance based adjecent genes.
"""
import argparse # reads commandline arguments

def commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, required=True, help="input file. output of CnvStatistics.jar")
    parser.add_argument("-c", "--columns", type=int, nargs="*", default=[6,7,8,9,10], help="Columns with list of ';' separated values in the input file.")
    parser.add_argument("-f", "--function", type=str, default='sum', help="Function to combine the scores of all genes. One of 'sum', 'mean'..")
    parser.add_argument("-t", "--threshold", type=float, default=0.0, help="Threshold t to consider only phenomatch scores >= t.")
    parser.add_argument("-k", "--power_parameter", type=float, default=1.0, help="Parameter k to boost large values. Each Gene score is taken to the power of k")
    parser.add_argument("-o", "--output_file", type=str, required=True, help="output file")
    return parser.parse_args()

def mean(l):
    """ returns the mean of a list with floats"""
    return sum(f) / max(1, len(f))
    
def phenogram_score(inFile, columns, func, t, k, outFile):
    """ """
    if func == "sum":
        f = sum
    elif func == "mean":
        f= mean
    elif func == "max":
        f = max
    else:
        print "ERROR: not proper function argument"
        exit()   

    with open(outFile, 'w') as outHandle:
        
        for line in open(inFile):
            
            # ignore error output with '*'
            if not line.startswith('*'):
            
                sp = line.strip().split('\t')
                    
                # pares scores, ignore NA
                overlap_scores = [parse_float(s) if s != 'NA' else 0.0 for s in sp[columns[0]].split(';')]
                up_scores = [parse_float(s)  if s != 'NA'else 0.0 for s in sp[columns[1]].split(';')]
                down_scores = [parse_float(s)  if not s == 'NA'else 0.0 for s in sp[columns[2]].split(';')]
                dist_up_scores = [parse_float(s)  if s != 'NA'else 0.0 for s in sp[columns[3]].split(';')]
                dist_down_scores = [parse_float(s)  if not s == 'NA'else 0.0 for s in sp[columns[4]].split(';')]
                
                # filter for phenomatch scores >= t:
                overlap_scores = [s**k for s in overlap_scores if s >= t]
                up_scores = [s**k for s in up_scores if s >= t]
                down_scores = [s**k for s in down_scores if s >= t]
                dist_up_scores = [s**k for s in dist_up_scores if s >= t]
                dist_down_scores = [s**k for s in dist_down_scores if s >= t]
                
                # combine scores with input function
                pheno_overlap = f(overlap_scores) if overlap_scores else 0
                
                # take max of upsteam and downstream adjacent scores
                pheno_adjacent = max(f(up_scores) if up_scores else 0, f(down_scores) if down_scores else 0)
                
                pheno_dist = max(f(dist_up_scores) if dist_up_scores else 0, f(dist_down_scores) if dist_down_scores else 0)

                # replace input columns with combined values:
                sp[columns[0]] = str(pheno_overlap)
                sp[columns[1]] = str(pheno_adjacent)
                #~ sp[8] = str(pheno_total)
                sp[columns[2]] = str(pheno_dist)
                
                # delete other columns with per gene scores:
                sp[columns[3]:(columns[4]+1)] = []
                                            
                #outHandle.write("{0}\t{1}\t{2}\n".format(pheno_overlap, pheno_adjacent, pheno_total))
                outHandle.write("\t".join(sp) + '\n')

def parse_float(s):
    """ reads a float form a string, even if it has the german format 0,123..."""
    return float(s.replace(',', '.'))

if __name__ == "__main__":

    # read commandline argumets
    args = commandline()
    
    phenogram_score(args.input_file, args.columns, args.function, args.threshold, args.power_parameter, args.output_file)
