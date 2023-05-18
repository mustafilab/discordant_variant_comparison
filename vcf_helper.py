from tqdm import tqdm

## function to compare depths on dataframes and filter at threshold t
def compare_depth(a,b,t):
    A     = {}
    B     = {}
    AB    = {}
    nAnB  = {}
    for i in tqdm(a.index.values, desc="comparing depths"):
        loc = a.loc[i,'x']
        if (a.loc[i,'y'] < t) and (b.loc[i,'y'] < t):
            nAnB[loc] = 1
        elif (a.loc[i,'y'] < t) and (b.loc[i,'y'] >= t):
            B[loc] = 1
        elif (a.loc[i,'y'] >= t) and (b.loc[i,'y'] < t):
            A[loc] = 1
        elif (a.loc[i,'y'] >= t) and (b.loc[i,'y'] >= t):
            AB[loc] = 1
        else:
            print("something unexpected happened during depth comparison")
    
    return A, B, AB, nAnB


def get_variant_depth(a,d):
    return len(set(a.keys()).intersection(set(d.keys())))

def get_variant_keys(a,d):
    return set(a.keys()).intersection(set(d.keys()))

def get_low_depth_loc(a,d):
    return list(set(a.keys()).intersection(set(d.keys())))