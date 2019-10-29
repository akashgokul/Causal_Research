##Helper functions

#Updates endogenous variables (V) given a certain state X,W
#i.e. sets to True/False given X or W
def V_given_set(V,X_or_W,indices):
    updated_V = V.copy()
    var_seen = 0
    for i in range(len(V)):
        if i in indices:
            updated_V[i] = X_or_W[var_seen]
            var_seen += 1
    return updated_V
#Goes through dictionary makes sure no repeat elements
def prune(dct):
    for key in dct.keys():
        curr_lst = dct[key]
        seen = []
        for elem in curr_lst:
            if elem not in seen:
                seen.append(elem)
        dct[key] = seen
    return dct
#Creates a list of all subsets (cardinality >= 1) of a given lst
#Used https://bit.ly/2YVGYxt
def subsets_finder(indices):
    #indices corresp. to elem of some lst opposed to using indices corresp to V.
    n = len(indices)
    if(n < 1):
        return indices
    else:
        subsets = []
        for i in range(n+1):
            subset = list(itertools.combinations(indices, i))
            subsets += [sub for sub in subset]
        return subsets

#returns list of all possible splits(take indices finds all )
#Used: https://bit.ly/2KfTeFk
def partitions(indices):
    if(len(indices) == 0):
        return [([],[])]
    subsets = [v for a in range(len(indices)) for v in itertools.combinations(indices, a)]
    comb = []
    for i in range(len(subsets)//2 + 1):
        comb.append((list(itertools.chain(subsets[i])), [e for e in indices if e not in subsets[i]]))
    ret_lst = comb + [(tup[1],tup[0]) for tup in comb]
    return ret_lst

#Given certain indices, updates lst[i] = not lst[i] for i in indices
def negate_var(lst,indices):
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


#Created a class of original H.P defn of Causality
class Causal_Model:
    #U,V,F are assumed to be lists for simplicity sake
    #Assumed that final var = outcome (e.g. last elem of V = Forest Fire)
    #U,V are boolean values corresp. to some list of var
    #F is a list of pointers to functions (e.g. def foo())
    def __init__(self, U, V, R, F):
        self.exogenous_var = U
        self.endogenous_var = V
        self.var_range = R
        self.function = F
        self.signature = (U,V,R)
        self.model = (self.signature,F)



    # For the following, I assume that this can be used to check causality if an outcome did not occur(?)

    #AC1 of defn, checks if outcome = function for given X
    #outcome_val refers to desired outcome
    def ac1_check(self,X,X_indices,outcome_function=self.F,outcome_val):
        V = V_given_set(self.endogenous_var,X,X_indices)
        outcome_given_X = outcome_function(self.exogenous_var,V)
        return outcome_given_X == outcome_val

    #AC2(a) checks the but-for clause, i.e changing X would lead to opposite outcome
    #2005 paper says W->w' but modified defn. paper says W->w, func is using the latter
    #This function finds the correct W, and calls ac2_b
    def ac2_a_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function=self.F,outcome_val):
        V = V_given_set(self.endogenous_var,Z,Z_indices)
        V = V_given_set(V,W,W_indices)
        x_prime = [not i for i in X]
        V_given_x_prime = V_given_set(V,x_prime,X_indices)
        outcome = outcome_function(self.exogenous_var,V_given_x_prime)
        if outcome_val != outcome:
            return True
        else:
            return False


    #Checks AC2(b) of the defn
    #Checks that outcome holds for all subsets of Z (Z') if Z' is set to original value
    def ac2_b_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function=self.F, outcome_val):
        V = self.endogenous_var
        V_fixed_W = V_given_set(V,W,W_indices)
        V = V_given_set(V_fixed_W,X,X_indices)
        subsets_of_Z = subsets_finder(list(range(len(Z))))
        orig_Z = [V[i] for i in Z_indices]
        for subset in subsets_of_Z:
            curr_Z = [orig_Z[i] if i in subset else Z[i] for i in range(len(Z))]
            curr_V = V_given_set(V_fixed_W,curr_Z,subset)
            outcome = outcome_function(self.exogenous_var,curr_V)
            if(outcome_val != outcome):
                return False
        return True

    #Checks AC2(b^u) of the defn
    #Checks that outcome holds for all subsets of W (W') and Z(Z') if W',Z' is set to original value
    def ac2_b_u_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function=self.F, outcome_val):
        V = self.endogenous_var
        orig_Z = [V[i] for i in Z_indices]
        orig_W = [V[i] for i in W_indices]
        subsets_of_Z = subsets_finder(list(range(len(Z))))
        subsets_of_W = subsets_finder(list(range(len(W))))
        for sub_Z in subsets_of_Z:
            curr_Z = [orig_Z[i] if i in sub_Z else Z[i] for i in range(len(Z))]
            for sub_W in subsets_of_W:
                curr_W = [orig_W[i] if i in sub_W else W[i] for i in range(len(W))]
                curr_V = V_given_set(V,curr_Z,Z_indices)
                curr_V = V_given_set(curr_V,curr_W,W_indices)
                outcome = outcome_function(self.exogenous_var,curr_V)
                if(outcome_val != outcome):
                    return False
        return True

    #Returns true if ac2_a^m is satisfied
    #Also gives the coresp. W,W_indices that satisfy the defn
    #Output: [W, W_indices, T/F]
    def ac2_m_check(self,X,X_indices,outcome_function=self.F,outcome_val):
        V = self.endogenous_var
        V_indices_excluding_X = [i for i in range(len(V)) if i not in X_indices]
        potential_W_indices = subsets_finder(V_indices_excluding_X)
        subsets_X = subsets_finder(X_indices)
        for sub in subsets_X:
            sub_lst = [i for i in sub]
            curr_X = negate_var(X,sub_lst)
            curr_V = V_given_set(V,curr_X,X_indices)
            for w_indices in potential_W_indices:
                curr_W = [V[i] for i in w_indices]
                curr_outcome = outcome_function(self.exogenous_var,curr_V)
                if(curr_outcome != outcome_val):
                    return[curr_W,w_indices,True]

        return [None,None,False]

    def ac2_check_given_Z_W(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function=self.F,outcome_val):
        ac2_a = self.ac2_a_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        ac2_b = self.ac2_b_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        if(ac2_a and ac2_b):
            return True
        else:
            return False

    def ac2_u_check_given_Z_W(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function=self.F,outcome_val):
        ac2_a = self.ac2_a_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        ac2_b = self.ac2_b_u_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        if(ac2_a and ac2_b):
            return True
        else:
            return False


    #Goes through all partitions of V (partitions named Z,W) and returns the first Z,W to satisfy
    #defn is an int(0,1) corresp to (original def, updated defn ('05))
    #Returns a dictionary of num_changes (from model to causal scenario) -> (corresp causal scenario's W)
    # If W_dct is empty => ac2 failed
    #If fast => return after finding first W that works
    #Incorporates num_outcomes to avoid putting outcome in W
    def Z_and_W_search(self,X,X_indices,outcome_function=self.F,outcome_val,defn,num_outcomes,fast):
        V = self.endogenous_var
        V_var = len(V)
        outcome_indices = [V_var-(i+1) for i in range(num_outcomes)]
        Useable_V_indices = [i for i in range(V_var) if ((i not in X_indices) and (i not in outcome_indices))]
        splits = partitions(Useable_V_indices)
        W_dct = {}
        for partition in splits:
            curr_Z_indices_no_X = partition[0] + V[:(-1*num_outcomes)]
            curr_W_indices = partition[1]
            curr_Z_indices = X_indices + curr_Z_indices_no_X
            curr_Z = [V[i] for i in curr_Z_indices]
            curr_W = [V[i] for i in curr_W_indices]
            ac2_check = self.ac2_check_given_Z_W(curr_Z,curr_Z_indices,X,X_indices,
                                                 curr_W,curr_W_indices,outcome_function,outcome_val)
            if(fast and ac2_check):
                return {0:[(curr_W,curr_W_indices)]}
            elif(ac2_check):
                if(0 in W_dct.keys()):
                    num_change_lst = W_dct[0]
                    num_change_lst.append((curr_W,curr_W_indices))
                    W_dct[0] = num_change_lst
                else:
                    W_dct[0] = [(curr_W,curr_W_indices)]
            else:
                subsets_of_curr_Z_no_X = subsets_finder(curr_Z_indices_no_X)
                subsets_of_curr_W = subsets_finder(curr_W_indices)
                for sub_z in subsets_of_curr_Z_no_X:
                    sub_zlst = [i for i in sub_z]
                    updated_Z = X + [not V[i] if i in sub_z else V[i] for i in curr_Z_indices_no_X]
                    for sub_w in subsets_of_curr_W:
                        sub_wlst = [i for i in sub_w]
                        updated_W = [V[i] if i in sub_w else not V[i] for i in curr_W_indices]
                        ac2_check = self.ac2_check_given_Z_W(updated_Z,curr_Z_indices,X,X_indices,updated_W,
                                                             curr_W_indices,outcome_function,outcome_val)
                        ac2_u_check = self.ac2_u_check_given_Z_W(updated_Z,curr_Z_indices,X,X_indices,updated_W,
                                                             curr_W_indices,outcome_function,outcome_val)
                        if((defn == 0 and ac2_check) or (defn == 1 and ac2_u_check)):
                            orig_W = [self.endogenous_var[i] for i in curr_W_indices]
                            num_changes = sum([1 if orig_W[i] != updated_W[i] else 0 for i in range(len(orig_W))])
                            if(fast):
                                W_dct[num_changes] = [(updated_W,curr_W_indices)]
                                return W_dct
                            elif(num_changes in W_dct):
                                num_changes_lst = W_dct[num_changes]
                                num_changes_lst.append((updated_W,curr_W_indices))
                                W_dct[num_changes] = num_changes_lst
                            else:
                                W_dct[num_changes] = [(updated_W,curr_W_indices)]

        W_dct= prune(W_dct)
        return W_dct





    #Checks that X is minimal by iterating over all subsets
    def ac3_check(self,X,X_indices,outcome_function=self.F,outcome_val):
        if(len(X) <= 1):
            return True
        subsets_of_X = subsets_finder(X_indices)
        for sub in subsets_of_X:
            sub_lst = [i for i in sub]
            updated_X = extract_X(X,sub_lst)
            W_indices = [j for j in range(len(self.endogenous_var)) if j not in sub]
            W = [self.endogenous_var[k] for k in W_indices]
            ac2_check = self.ac2_a_check(updated_X,sub_lst,updated_X,sub_lst,W,W_indices,outcome_function,outcome_val)
            if(ac2_check):
                return False

        return True

    #Finding paths to the model s.t. X is pivotal
    #(I'm writing this b/c using W generated for zultan is having problems due to multitude of Z,W pairs)
    def paths_to_pivotality(self,X,X_indices,outcome_val,outcome_func=self.F,num_outcomes):
        V = self.endogenous_var
        V_no_X_indx = [i for i in range(len(V)) if i not in X_indices]
        V_no_X = [V[i] for i in V_no_X_indx]
        paths = []
        subsets_V_no_X = subsets_finder(V_no_X_indx)
        subsets_V_no_X = [tup for tup in subsets_V_no_X if len(tup) > 0]
        for setting in subsets_V_no_X:
            curr_V_no_X = [not V_no_X[i] if i in setting else V_no_X[i] for i in range(len(V_no_X_indx))]
            curr_V = []
            curr_V_prime = []
            X_seen = 0
            for i in range(len(V)):
                if(i in range(len(X_indices))):
                    X_seen += 1
                    curr_V.append(X[i])
                    curr_V_prime.append(not X[i])
                else:
                    curr_V.append(curr_V_no_X[i - X_seen])
                    curr_V_prime.append(curr_V_no_X[i-X_seen])
            outcome_one = outcome_func(self.exogenous_var,curr_V)
            outcome_prime = outcome_func(self.exogenous_var,curr_V_prime)
            outcome_check = outcome_one == outcome_val
            pivotality_check = outcome_one != outcome_prime
            if(outcome_check and pivotality_check):
                paths.append((len(setting),curr_V))
        return paths

    #Fail checking
    def wrong_check(self,ac_1,ac_2,ac_3):
        if (not ac_1):
            print("(False b/c of AC1)")
        if(not ac_2):
            print("(False b/c of AC2)")
        if(not ac_3):
            print("(False b/c of AC3)")

    #Returns true if X satisfies HP defn, False o.w.
    def causality_check(self,X,X_indices,outcome_val,outcome_func=self.F,num_outcomes,fast):
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        dict_to_bool = lambda x: False if len(x) == 0 else True
        ac_2 = dict_to_bool(self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,0,num_outcomes,fast))
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    #Returns true if X satisfied the modified defn. (using ac1,ac2(a),ac2(b^u), ac3 - Halpern & Pearl 2005)
    def updated_causality_check(self,X,X_indices,outcome_val,outcome_func=self.F,num_outcomes,fast):
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        dict_to_bool = lambda x: False if len(x) == 0 else True
        ac_2 = dict_to_bool(self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,1,num_outcomes,fast))
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    #Returns true if X satisfies modified defn of Halpern & Pearl 2014 (ac1, ac2_m, ac3)
    def modified_causality_check(self,X,X_indices,outcome_val,outcome_func=self.F):
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        ac_2 = self.ac2_m_check(X,X_indices,outcome_func,outcome_val)[2]
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    #Returns "responsibility" as per Chockler & Halpern (2004)
    # i.e. Calls Z_and_W_search and returns the min_key
    def responsibility(self,X,X_indices,outcome_val,outcome_func=self.F,num_outcomes,fast):
        if(self.updated_causality_check(X,X_indices,outcome_val,outcome_func,num_outcomes,fast)):
            not_X = [not i for i in X]
            if(outcome_func(self.exogenous_var,V_given_set(self.endogenous_var,not_X,X_indices))!=outcome_val):
                return 1
            else:
                W_dct = self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,0,num_outcomes,False)
                min_num_changes = min(W_dct.keys())
                return 1 / (min_num_changes + 1)
        else:
            return 0


    #num_outcome_var denote var in V which are not part of structural eqn. (e.g. Forest Fire)
    #Goes through all the possible scenarios of variables values (excluding X)
    # For each scenario, checks if f(X U Scenario) != f(not X U Scenario)
    # Returns sum of outcome_changes / num_scenarios
    def influence(self,X,X_indices,num_outcome_var,outcome_func=self.F):
        U = self.exogenous_var
        end_idx = (-1)*num_outcome_var
        V = self.endogenous_var
        V_indices_no_X = [i for i in range(len(V[:end_idx])) if i not in X_indices]
        subsets_of_V = subsets_finder(V_indices_no_X)
        if(len(subsets_of_V)==0):
            return 1
        outcome_change_ct = 0
        for subset in subsets_of_V:
            non_X_var_values = [True if i in subset else False for i in V_indices_no_X]
            V_scenario = V_given_set(V,non_X_var_values,V_indices_no_X)
            V_pos_X_scenario = V_given_set(V_scenario,X,X_indices)
            V_neg_X_scenario = V_given_set(V_scenario,[not i for i in X], X_indices)
            if(outcome_func(U,V_pos_X_scenario) != outcome_func(U,V_neg_X_scenario)):
                outcome_change_ct += 1
        return outcome_change_ct / len(subsets_of_V)

    #Adj responsibility = influence*responsibility
    def adj_responsibility(self,X,X_indices,num_outcome_var,outcome_val,outcome_func=self.F,fast):
        inf = self.influence(X,X_indices,num_outcome_var,outcome_func)
        res = self.responsibility(X,X_indices,outcome_val,outcome_func,num_outcome_var,fast)
        return inf*res

    #Adj responsibility_2 = sum of influence of var in W that change
    #If multiple W, computes the above ^ for each W and sums together
    #Returns 0 if X not Causal

    #Changing this to responsibility in every W generated world, using regularization (divide by num_W_generated)
    def adj_responsibility_2(self,X,X_indices,num_outcome_var,outcome_val,outcome_func=self.F,fast):
        if(self.updated_causality_check(X,X_indices,outcome_val,outcome_func,num_outcome_var,fast)):
            not_X = [not i for i in X]
            if(outcome_func(self.exogenous_var,V_given_set(self.endogenous_var,not_X,X_indices))!=outcome_val):
                return 1
            else:
                W_dct = self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,0,num_outcome_var,fast)
                W_dct_keys = W_dct.keys()
                min_key = min(W_dct_keys)
                W_indices = W_dct[min_key][0][1]
                orig_W = [self.endogenous_var[i] for i in W_indices]
                inf = self.influence(orig_W,W_indices,num_outcome_var,outcome_func)
                return 1 / (1+len(W_dct[min_key])*min_key*inf)
        else:
            return 0

    #This function returns adj_2 and zultan at once (to save time)
    def adj_responsibility_2m(self,X,X_indices,num_outcome_var,outcome_val,outcome_func=self.F,fast):
        if(not self.updated_causality_check(X,X_indices,outcome_val,outcome_func,num_outcome_var,True)):
            return (0,0)

        else:
            not_X = [not i for i in X]
            if(outcome_func(self.exogenous_var,V_given_set(self.endogenous_var,not_X,X_indices))!=outcome_val):
                return (1,1)
            else:
                W_dct = self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,1,num_outcome_var,False)
                W_dct_keys = W_dct.keys()
                min_key = min(W_dct_keys)
                W_indices = W_dct[min_key][0][1]
                orig_W = [self.endogenous_var[i] for i in W_indices]
                inf = self.influence(orig_W,W_indices,num_outcome_var,outcome_func)
                adj_2_outcome = 1 / (1+len(W_dct[min_key])*min_key)
        #zultan:
                paths = self.paths_to_pivotality(X,X_indices,outcome_val,outcome_func,num_outcome_var)
                if(len(paths) == 0):
                    return (adj_2_outcome,1)
                N = 1 / (sum([1 / p[0] for p in paths]))
                return (adj_2_outcome, 1 / (N+1))


    #Using Multiple Counterfactual Pivotality model from Zultan, Gerstenberg, Lagnado 2012
    # If there are multiple ways that a variable X can become causal (via but-for) then instead of
    # using the 1 / (1 + min(changes_to_W)) we will change the denom to (1+N)
    # N = 1 / (sum of 1/c_i) where c_i is the number of changes to W in "path" i
    #this iterates over all possible paths where X is pivotal
    #If X is pivotal : returns 1
    #Returns 0 if X is not causal
    def zultan_responsibility(self,X,X_indices,outcome_val,outcome_func= self.F,num_outcome_var):
        if(self.causality_check(X,X_indices,outcome_val,outcome_func,num_outcome_var,True)):
            not_X = [not i for i in X]
            if(outcome_func(self.exogenous_var,V_given_set(self.endogenous_var,not_X,X_indices)) != outcome_val):
                return 1
            paths = self.paths_to_pivotality(X,X_indices,outcome_val,outcome_func,num_outcome_var)
            if(len(paths) == 0):
                return 1
            N = 1 / (sum([1 / p[0] for p in paths]))
            return 1 / (N+1)
        else:
            return 0
    #Influence function but only samples n states (assuming n >= 1)
    #Saves time as inf function is exponential ; mc = Monte Carlo method to sample from uniform [0,1]
    def mc_inf_sample(self,X,X_indices,num_outcome_var,outcome_func=self.F,n):
        if(n <= 0):
            raise ValueError
        V = self.endogenous_var
        v_idx_to_change = [i for i in range(len(V[0:(-1)*num_outcome_var])) if i not in X_indices]
        outcome_change_ct = 0
        seen = []
        for i in range(n):
            random_vec = np.random.randint(2,size = len(v_idx_to_change))
            changeable_var_vec = [True if random_vec[i] == 1 else False for i in range(len(random_vec))]
            curr_V_pos = V_given_set(V,changeable_var_vec,v_idx_to_change)
            curr_V_neg = negate_var(curr_V_pos,X_indices)
            outcome_pos = outcome_func(self.exogenous_var,curr_V_pos)
            outcome_neg = outcome_func(self.exogenous_var,curr_V_neg)
            if(outcome_pos != outcome_neg):
                outcome_change_ct += 1

        return outcome_change_ct / n
