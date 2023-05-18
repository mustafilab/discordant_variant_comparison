import re, numpy as np, pandas as pd

## class for assignment of vcf data to node in tree
## WARNING parse function requires modification to operate generally

ush2a_start  = 215622891
ush2a_end    = 216423448

class node:
    def __init__(self, filename, _id, dtype):
        self.id       = _id
        self.dtype    = dtype
        self.df       = self.load(filename)
        self.branches = {}
    
    def load(self, f):
        df = self.parse(\
                      self.parse(pd.read_csv(f, \
                                             header=0, \
                                             index_col=1, \
                                             sep="\t").rename(columns = {'ush2a_826':'VALUES', 'ush2a_827':'VALUES', 'ush2a_828':'VALUES'}), \
                            encoding="FORMAT", \
                            target="VALUES", \
                            sep=":",\
                            droplist=[]), \
                      encoding="INFO", \
                      target="INFO", \
                      sep=";", \
                      split_encoding="=", \
                      split_target="=", \
                      split_encoding_index=0,  \
                      split_target_index=1)
        return df
    
    def parse(self, \
              df, \
              encoding="FORMAT", \
              target="", \
              sep=":", \
              split_encoding=None, \
              split_target=None, \
              split_encoding_index=0,  \
              split_target_index=1, \
              droplist=[]):
        
        tmp_df    = df.copy().rename(columns = {'#CHROM':'chr'})
        field_map = {\
                 "AA":"ancestral_allele",\
                 "AC":"allele count",\
                 "AN":"total_allele_number",\
                 "AF":"allele_frequency",\
                 "AC":"allele_count",\
                 "BQ":"base_quality",\
                 "DB":"dbSNP_membership",\
                 "DP":"combined_depth",\
                 "MQ":"RMS_map_quality",\
                 "NS":"number_of_samples",\
                 "GT":"genotype",\
                 "DP":"depth",\
                 "AD":"allele_depth",\
                 "GL":"genotype_likelihood",\
                 "HQ":"haplotype_qualities",\
                 "C":"caller",\
                 "PS":"phase_set",\
                 "PL":"phase_likelihood",\
                 "PQ":"phasing_quality",\
                 "QD":"QD",\
                 "GQ":"conditional_genotype_quality",\
                 "VAF":"variant_allele_frequency"\
                }

        # split value field
        def value_parser(key, value):
            if split_target is None:
                v = value
            else:
                if len(value.split(split_target)) >= split_target_index + 1:
                    v = value.split(split_target)[split_target_index]
                else:
                    v = value

            if v == "DV" and key == "C":
                return "DeepVariant"
            elif v == "P" and key == "C":
                return "Pepper"
            elif v == np.nan:
                return pd.NA
            else:
                return v

        # split encoding field
        def key_parser(key):
            if split_encoding is None:
                k = key
            else:
                if len(key.split(split_encoding)) >= split_encoding_index + 1:
                    k = key.split(split_encoding)[split_encoding_index]
                else:
                    k = key

            if k in field_map:
                return field_map[k]
            else:
                return k
        
        # safe min for empty array or list
        def safe_min(x):
            if type(x) is list or type(x) is type(np.array([])):
                if len(x) > 0:
                    return min(x)
                elif len(x) == 0:
                    return 0
            else:
                try:
                    return min(x)
                except:
                    print("unhandled safe min error: returning 0")
                    return 0
                
        def is_indel(var):
            return True if ( (str(var).strip().upper() == "-") or (len(str(var).strip()) ) > 1) else False
        
        def is_snp(var):
            return True if ( (str(var).strip().upper() in ['A','C','G','T']) and (len(var) == 1) ) else False
    
    
        # load values into tmp_df        
        for i in df.index.values:
            for key,value in zip(df.loc[i,encoding].split(sep), df.loc[i,target].split(sep)):
                # add entry
                tmp_df.loc[i, key_parser(key)] = value_parser(key, value)

                
                
        # drop parsed columns
        tmp_df.drop(labels=set(np.unique([encoding, target])).union(*droplist),axis=1, inplace=True)
        
        # filter for variants in ush2a
        tmp_df = tmp_df.loc[ush2a_start:ush2a_end, :]
        
        # filter for variants that pass
        tmp_df = tmp_df.loc[tmp_df.FILTER=="PASS", :]
        
        
        
        ### FILTERING
        ###
        ### WARNING: tuning parameters can provide higher concordance at the risk of overfitting to the specific callset
        ### Multiple statistical approaches may be used, ultimately the baseline is compared.
        ###
        
        # init selection to slice up tmp_df -> only passed and in ush2a
        sel = (tmp_df.FILTER=="PASS")
        
        # filter for variants with allele depth >= n
        min_allele_depth = 4
#        sel &= np.array([ bool(safe_min([int(a) for a in x.split(",") if int(a) != 0]) >= min_allele_depth) for x in tmp_df.allele_depth.values ])
        
        # filter for variants with depth >= n
        min_depth = 10
#        sel &= np.array([ bool(x >= min_depth) for x in tmp_df.depth.values.astype(float)])
        
        
        
        # filter for variants with conditional_genotype_quality >= 3 for nanopore, 36 for illumina
        # values from rtg vcfeval bestfit thresholding
        if self.dtype=="nanopore":      threshold = 3  ## 3
        elif self.dtype=="illumina":    threshold = 36 ## 36 
#        sel &= (tmp_df.conditional_genotype_quality.astype(float) >= threshold)

        # filter for variants with conditional_genotype_quality >= quantile
        if self.dtype=="nanopore":      quantile  = 0.01
        elif self.dtype=="illumina":    quantile  = 0.01
        threshold = tmp_df.conditional_genotype_quality.astype(float).quantile(quantile)
#        sel &= (tmp_df.conditional_genotype_quality.astype(float) >= threshold)
        
    
        # create indel mask
        indel_mask = []
        for i in tmp_df.index.values:
            var = str(tmp_df.loc[i,'ALT']).strip().split(",")
            ref = str(tmp_df.loc[i,'REF']).strip()
            mask = False
            for v in var:
                if ( (is_indel(v)) or (is_indel(ref)) ):
                    mask = True
            indel_mask.append(mask)
        
        
        # conditionally include indels
        ## quantile
        if self.dtype=="nanopore":      quantile        = 0.01
        elif self.dtype=="illumina":    quantile        = 0.01
        indel_threshold = tmp_df.loc[indel_mask, :].conditional_genotype_quality.astype(float).quantile(quantile)
        ## hard threshold
        if self.dtype=="nanopore":      indel_threshold    = 0
        elif self.dtype=="illumina":    indel_threshold    = 0
#        sel |= ( indel_mask & (tmp_df.conditional_genotype_quality.astype(float) >= indel_threshold) )
        
    
        tmp_df = tmp_df.loc[sel, :]
        
        return tmp_df
    
    
    
    def intersection(self, node):
        overlap = branch(self, node)
        self.branches[f'{node.id}_{node.dtype}'] = overlap
        node.branches[f'{self.id}_{self.dtype}'] = overlap

        
class branch:
    
    ####################################################################################################################
    # A             
    # |------|
    #        |      [loc]
    #      [snps]----|
    #        |       |  
    # |------|       |
    # B              |   
    #                |
    #                |   differ
    #                #------------------- [data = [(a_0,b_0,ref,), (a_1,...)...]
    #                |               
    #          same  |               
    #                |               
    #                |               
    #                |               
    #                |               
    #                -----------------[variant]            
    #                                     |               
    #                                     |               
    #                                     |  differ
    #                                     #------------------- [data = [(a_0,b_0,ref,), (a_1,...)...]
    #                                     |                   
    #                               same  |                   
    #                                     |               
    #                                     |
    #                                     -----------------[variant]
    #                                                          |
    #                                                          |
    #                                                          |  differ
    #                                                          #------------------- [data=(a,b,ref,phase_a,phase_b)]
    #                                                          |                   
    #                                                    same  |                   
    #                                                          |               
    #                                                          |
    ####################################################################################################################
    
    
    def is_snp(self, var):
        return True if ( (str(var).strip().upper() in ['A','C','G','T']) and (len(var) == 1) ) else False

    def is_indel(self, var):
        return True if ( (str(var).strip().upper() == "-") or (len(str(var).strip()) ) > 1) else False

    def phase(self,a):
        l = list(filter(lambda x: x!='', re.split(r"[\|/]",a)))
        if len(np.unique(l)) == 1: 
            return 'homozygous'
        else:        
            return 'heterozygous'

    def compare_strings(self,a,b):

        def string_iou(x,y):
            xas      = list(filter(lambda x: x!='', re.split("",x)))
            xbs      = list(filter(lambda x: x!='', re.split("",y)))
            n_union  = max(len(xas),len(xbs))
            n_ident  = len( list(filter(lambda x: x[0] == x[1], zip(xas, xbs))) )


            ua,na = np.unique(xas,return_counts=True)
            ub,nb = np.unique(xbs,return_counts=True)

            ha    = dict(zip(ua,na))
            hb    = dict(zip(ub,nb))

            n_noorder_ident = 0
            for k in set(ha.keys()).intersection(set(hb.keys())):
                n_noorder_ident += min(ha[k],hb[k])            

            return max(n_ident, n_noorder_ident) / n_union

        xa  = a.split(",")
        xb  = b.split(",")
        res = []
        for i,va in enumerate(xa):
            a_res = []
            for j,vb in enumerate(xb):
                a_res.append((va,vb,string_iou(va,vb)))
            res.append(sorted(a_res, key=lambda x: x[2], reverse=True))
        res = sorted(res,key=lambda x: x[0][2], reverse=True)

        n = min(len(xa),len(xb))

        ret = []
        for l in res[0:n]:
            ret.append(l[0])
        return ret
    
    def __init__(self, node_a, node_b):
        self.a                                     = f'{node_a.id}_{node_a.dtype}'
        self.b                                     = f'{node_b.id}_{node_b.dtype}'
        self.loc_union                             = {}
        self.loc_intersection                      = {}
        self.loc_unique                            = {}
        self.snp                                   = {}
        self.indel                                 = {}
        self.overlap = {\
                        "snp":{\
                             "variant_agree":{\
                                            "phase_agree":{},\
                                            "phase_disagree":{}\
                                           },\
                             "variant_disagree":{},\
                             "loc_disagree":{}\
                            },\
                        "indel":{\
                             "variant_agree":{\
                                            "phase_agree":{},\
                                            "phase_disagree":{}\
                                           },\
                             "variant_disagree":{},\
                             "loc_disagree":{}\
                            },\
                        }
        self.stats = {\
                        "snp":{\
                             "variant_agree":{\
                                            "phase_agree":0,\
                                            "phase_disagree":0\
                                           },\
                             "variant_disagree":{\
                                            "phase_agree":0,\
                                            "phase_disagree":0\
                                           },\
                             "loc_disagree":{}\
                            },\
                        "indel":{\
                             "variant_agree":{\
                                            "phase_agree":0,\
                                            "phase_disagree":0\
                                           },\
                             "variant_disagree":{\
                                            "phase_agree":0,\
                                            "phase_disagree":0\
                                           },\
                             "loc_disagree":{}\
                            },\
                        }
                
       # def calculate_overlap(node_a, node_b):            
        self.loc_union           = set(node_a.df.index.values).union(set(node_b.df.index.values))
        self.loc_intersection    = set(node_a.df.index.values).intersection(set(node_b.df.index.values))
        self.loc_unique[self.a]  = set(node_a.df.index.values)-self.loc_intersection
        self.loc_unique[self.b]  = set(node_b.df.index.values)-self.loc_intersection



        # load snps and indels a
        for i in node_a.df.index.values:
            if not i in self.snp:    self.snp[i]   = {}
            if not i in self.indel:  self.indel[i] = {}

            var = str(node_a.df.loc[i,'ALT']).strip().split(",")
            ref = str(node_a.df.loc[i,'REF']).strip()

            for v in var:
                if ( (self.is_snp(v)) and (self.is_snp(ref)) ):
                    if not self.a in self.snp[i]:
                        self.snp[i][self.a]    = {v: self.phase(node_a.df.loc[i,'genotype'])}
                    else:
                        self.snp[i][self.a][v] = self.phase(node_a.df.loc[i,'genotype'])

                elif ( (self.is_indel(v)) or (self.is_indel(ref)) ):
                    if not self.a in self.indel[i]:
                        self.indel[i][self.a]    = {v: self.phase(node_a.df.loc[i,'genotype'])}
                    else:
                        self.indel[i][self.a][v] = self.phase(node_a.df.loc[i,'genotype'])

        # load snps and indels b
        for i in node_b.df.index.values:
            if not i in self.snp:   self.snp[i]   = {}
            if not i in self.indel: self.indel[i] = {}

            var = str(node_b.df.loc[i,'ALT']).strip().split(",")
            ref = str(node_b.df.loc[i,'REF']).strip()
                      
            for v in var:
                if ( (self.is_snp(v)) and (self.is_snp(ref)) ):
                    if not self.b in self.snp[i]:
                        self.snp[i][self.b]    = {v: self.phase(node_b.df.loc[i,'genotype'])}
                    else:
                        self.snp[i][self.b][v] = self.phase(node_b.df.loc[i,'genotype'])

                elif ( (self.is_indel(v)) or (self.is_indel(ref)) ):
                    if not self.b in self.indel[i]:
                        self.indel[i][self.b]    = {v: self.phase(node_b.df.loc[i,'genotype'])}
                    else:
                        self.indel[i][self.b][v] = self.phase(node_b.df.loc[i,'genotype'])


        # calculate snp overlap
        for i,ht in self.snp.items():
            # loc agree
            if len(ht.keys()) == 2:
                vp_a = list(ht.values())[0]
                vp_b = list(ht.values())[1]

                # variants agree
                for v in set(vp_a.keys()).intersection(vp_b.keys()):
                    # phase agree
                    if vp_a[v] == vp_b[v]:
                        if not i in self.overlap['snp']['variant_agree']['phase_agree']:
                            self.overlap['snp']['variant_agree']['phase_agree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['snp']['variant_agree']['phase_agree'][i] += [(v, vp_a[v])]
                    # phase disagree
                    else:
                        if not i in self.overlap['snp']['variant_agree']['phase_disagree']:
                            self.overlap['snp']['variant_agree']['phase_disagree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['snp']['variant_agree']['phase_disagree'][i] += [(v, vp_a[v])]

                # variants disagree
                for v in set(vp_a.keys()).symmetric_difference(vp_b.keys()):
                    if v in vp_a:
                        if not i in self.overlap['snp']['variant_disagree']:
                            self.overlap['snp']['variant_disagree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['snp']['variant_disagree'][i] += [(v, vp_a[v])]
                    elif v in vp_b:
                        if not i in self.overlap['snp']['variant_disagree']:
                            self.overlap['snp']['variant_disagree'][i] = [(v, vp_b[v])]
                        else:
                            self.overlap['snp']['variant_disagree'][i] += [(v, vp_b[v])]
            # location disagree
            else:
                for k, vp in ht.items():
                    for v in vp.keys():
                        if not k in self.overlap['snp']['loc_disagree']:
                            self.overlap['snp']['loc_disagree'][k] = {i:[(v, vp[v])]}

                        elif not i in self.overlap['snp']['loc_disagree'][k]:
                            self.overlap['snp']['loc_disagree'][k][i] = [(v, vp[v])]
                        else:
                            self.overlap['snp']['loc_disagree'][k][i] += [(v, vp[v])]

        # calculate indel overlap
        for i,ht in self.indel.items():
            # loc agree
            if len(ht.keys()) == 2:
                vp_a = list(ht.values())[0]
                vp_b = list(ht.values())[1]

                # variants agree
                for v in set(vp_a.keys()).intersection(vp_b.keys()):
                    # phase agree
                    if vp_a[v] == vp_b[v]:
                        if not i in self.overlap['indel']['variant_agree']['phase_agree']:
                            self.overlap['indel']['variant_agree']['phase_agree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['indel']['variant_agree']['phase_agree'][i] += [(v, vp_a[v])]
                    # phase disagree
                    else:
                        if not i in self.overlap['indel']['variant_agree']['phase_disagree']:
                            self.overlap['indel']['variant_agree']['phase_disagree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['indel']['variant_agree']['phase_disagree'][i] += [(v, vp_a[v])]

                # variants disagree
                for v in set(vp_a.keys()).symmetric_difference(vp_b.keys()):
                    if v in vp_a:
                        if not i in self.overlap['indel']['variant_disagree']:
                            self.overlap['indel']['variant_disagree'][i] = [(v, vp_a[v])]
                        else:
                            self.overlap['indel']['variant_disagree'][i] += [(v, vp_a[v])]
                    elif v in vp_b:
                        if not i in self.overlap['indel']['variant_disagree']:
                            self.overlap['indel']['variant_disagree'][i] = [(v, vp_b[v])]
                        else:
                            self.overlap['indel']['variant_disagree'][i] += [(v, vp_b[v])]
            
            # location disagree
            else:
                for k, vp in ht.items():
                    for v in vp.keys():
                        if not k in self.overlap['indel']['loc_disagree']:
                            self.overlap['indel']['loc_disagree'][k] = {i:[(v, vp[v])]}

                        elif not i in self.overlap['indel']['loc_disagree'][k]:
                            self.overlap['indel']['loc_disagree'][k][i] = [(v, vp[v])]
                        else:
                            self.overlap['indel']['loc_disagree'][k][i] += [(v, vp[v])]