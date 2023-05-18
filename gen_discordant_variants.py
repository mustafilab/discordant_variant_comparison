import re, numpy as np, pandas as pd
from scipy import stats
from tqdm import tqdm

from vcf_node import node
from vcf_helper import *


## declare files and parameters
ush2a_start  = 215622891
ush2a_end    = 216423448
n826_vcf     = "data/ush2a_826.nanopore.phased.filtered.vcf"
i826_vcf     = "data/ush2a_826.illumina.hg38.varfilter.vcf.txt"
n826_bed     = "data/ush2a_826.assembly.haplotagged.depth.tsv"
i826_bed     = "data/ush2a_826.hg38.aln.depth.tsv"



### IMPORT

## bed files
nbed = pd.read_csv(n826_bed, header=None, names=['chr','x','y'], index_col=None, sep="\t")
ibed = pd.read_csv(i826_bed, header=None, names=['chr','x','y'], index_col=None, sep="\t")
## nodes
n826 = node(n826_vcf,"826","nanopore")
i826 = node(i826_vcf,"826","illumina")


### COMPARISON

## run depth comparison on bed files
A, B, AB, nAnB = compare_depth(nbed,ibed,10)
## calculate intersections
n826.intersection(i826)



#### EXPORT

## header is ['sequencing platform','partner_coverage_low','type',(vcf header)...]
with open('depth_comparison.results.np.tsv','w') as f:
    header = f"sequencing platform\tpartner_coverage_low\tvariant_type\tLOC\t" + "\t".join(n826.df.columns.values.tolist()) + "\n"
    f.write(header)
    
    ### nanopore
    ## indels
    a = get_variant_keys( n826.branches['826_illumina'].overlap['indel']['loc_disagree']['826_nanopore'], A)
    b = n826.branches['826_illumina'].overlap['indel']['loc_disagree']['826_nanopore'].keys()
    for v in b:
        if v in a:
            l = f"nanopore\tTrue\tindel\t{v}\t" + "\t".join(list(n826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
        else:
            l = f"nanopore\tFalse\tindel\t{v}\t" + "\t".join(list(n826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
            
    ## snps
    a = get_variant_keys( n826.branches['826_illumina'].overlap['snp']['loc_disagree']['826_nanopore'], A)
    b = n826.branches['826_illumina'].overlap['snp']['loc_disagree']['826_nanopore'].keys()
    for v in b:
        if v in a:
            l = f"nanopore\tTrue\tsnp\t{v}\t" + "\t".join(list(n826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
        else:
            l = f"nanopore\tFalse\tsnp\t{v}\t" + "\t".join(list(n826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
    
    f.close()
    
    
## header is ['sequencing platform','partner_coverage_low','type', (vcf header)...]
with open('depth_comparison.results.il.tsv','w') as f:
    header = f"sequencing platform\tpartner_coverage_low\tvariant_type\tLOC\t" + "\t".join(list(i826.df.columns.values)) + "\n"
    f.write(header)
    
    ### illumina
    ## indels
    a = get_variant_keys( n826.branches['826_illumina'].overlap['indel']['loc_disagree']['826_illumina'], B)
    b = n826.branches['826_illumina'].overlap['indel']['loc_disagree']['826_illumina'].keys()
    for v in b:
        if v in a:
            l = f"illumina\tTrue\tindel\t{v}\t" + "\t".join(list(i826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
        else:
            l = f"illumina\tFalse\tindel\t{v}\t" + "\t".join(list(i826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
            
    ## snps
    a = get_variant_keys( n826.branches['826_illumina'].overlap['snp']['loc_disagree']['826_illumina'], B)
    b = n826.branches['826_illumina'].overlap['snp']['loc_disagree']['826_illumina'].keys()
    for v in b:
        if v in a:
            l = f"illumina\tTrue\tsnp\t{v}\t" + "\t".join(list(i826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
        else:
            l = f"illumina\tFalse\tsnp\t{v}\t" + "\t".join(list(i826.df.loc[v,:].values.astype(str))) + "\n"
            f.write(l)
            
    f.close()

