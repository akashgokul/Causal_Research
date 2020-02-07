import itertools

#Helper functions
def V_given_set(V,X_or_W,indices):
    """
    Input: V (set of endogenous var), list corresp to X or W values, list indices to change

    Output: Updates values in V at select indices to the value of X/W list
    """
    updated_V = V.copy()
    var_seen = 0
    for i in range(len(V)):
        if i in indices:
            updated_V[i] = X_or_W[var_seen]
            var_seen += 1
    return updated_V

def prune(dct):
    """
    Input: A dictionary

    Output: Same dictionary without any repeat elements

    """
    for key in dct.keys():
        curr_lst = dct[key]
        seen = []
        for elem in curr_lst:
            if elem not in seen:
                seen.append(elem)
        dct[key] = seen
    return dct

def subsets_finder(indices):
    """
    Input: Indices of some list

    Output: A list of all subsets (cardinality >= 1) of a given list

    Reference: https://bit.ly/2YVGYxt
    """
    n = len(indices)
    if(n < 1):
        return indices
    else:
        subsets = []
        for i in range(n+1):
            subset = list(itertools.combinations(indices, i))
            subsets += [sub for sub in subset]
        return subsets

def partitions(indices):
    """
    Input: List of indices

    Output: Set of all partitions

    Note: I used this when the function above was broken
    Reference: https://bit.ly/2KfTeFk
    """
    if(len(indices) == 0):
        return [([],[])]
    subsets = [v for a in range(len(indices)) for v in itertools.combinations(indices, a)]
    comb = []
    for i in range(len(subsets)//2 + 1):
        comb.append((list(itertools.chain(subsets[i])), [e for e in indices if e not in subsets[i]]))
    ret_lst = comb + [(tup[1],tup[0]) for tup in comb]
    return ret_lst

def negate_var(lst,indices):
    """
    Input: List, list of select indices

    Updates list[i] for i in indices to be not list[i

    Output: Updated/Negated List
    """
    return_lst = lst.copy()
    for i in range(len(lst)):
        if(i in indices):
            return_lst[i] = not lst[i]
    return return_lst

#Extracts given elements from X using t
def extract_X(X,t):
    return_lst = []
    for i in t:
        return_lst.append(X[i])
    return return_lst
