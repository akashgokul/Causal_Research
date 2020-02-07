from utils import *
import numpy as np

class Causal_Model:
    def __init__(self, U, V, R, F):
        """
        Inputs: U: list of values of exogenous variables
                V: list of values of endogenous variables
                R: list of range of each variable (we deal with (True, False))
                F: list of functions corresponding to each variable in V

        Note: All variables are boolean valued. We assume the final variable of V is the outcome.

        """
        self.exogenous_var = U
        self.endogenous_var = V
        self.var_range = R
        self.function = F
        self.signature = (U,V,R)
        self.model = (self.signature,F)


    def ac1_check(self,X,X_indices,outcome_function,outcome_val):
        """
        Input: X: list of values for the desired input
               X_indices: Indices of X within the self.V (list of endogenous variables)
               Outcome_function: outcome function
               Outcome_val: Desired outcome check

        Note: This is AC1 condition from original Halper-Pearl definition.
        It checks that the desired outcome actually occurs for a given state of X

        Returns: True/False correponding to whether or not these X variables lead to outcome_val
        when put into the outcome_funtion.

        """
        V = V_given_set(self.endogenous_var,X,X_indices)
        outcome_given_X = outcome_function(self.exogenous_var,V)
        return outcome_given_X == outcome_val

    def ac2_a_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val):
        """
        Input: Z: A list of variable settings where the variables are in V and X ⊆ Z
                Z_indices: (list) Indices of Z variables, where these indices correspond to indices in self.V
                X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                W: Variables values for the complementary set to Z
                W_indices: (list) Indices of W variables in self.V
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test

        Note: This is the AC2(a) check of the original and updated HP definitions. It essentially checks but-for,
        i.e if X is set to the opposite value would the outcome still have happened?

        Returns: True/False corresponding to whether the outcome would differ from outcome_val if each x = not x

        """
        V = V_given_set(self.endogenous_var,Z,Z_indices)
        V = V_given_set(V,W,W_indices)
        x_prime = [not i for i in X]
        V_given_x_prime = V_given_set(V,x_prime,X_indices)
        outcome = outcome_function(self.exogenous_var,V_given_x_prime)
        if outcome_val != outcome:
            return True
        else:
            return False


    def ac2_b_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function, outcome_val):
        """
        Input: Z: A list of variable settings where the variables are in V and X ⊆ Z
                Z_indices: (list) Indices of Z variables, where these indices correspond to indices in self.V
                X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                W: Variables values for the complementary set to Z
                W_indices: (list) Indices of W variables in self.V
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test

        Note: This checks AC2(b) of the HP definitions.
        This checks that the outcome occurs even if any subset of Z is set to the actual values from the experiment.

        Returns: True/False if condiion is met

        """
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

    def ac2_b_u_check(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function, outcome_val):
        """
        Input: Z: A list of variable settings where the variables are in V and X ⊆ Z
                Z_indices: (list) Indices of Z variables, where these indices correspond to indices in self.V
                X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                W: Variables values for the complementary set to Z
                W_indices: (list) Indices of W variables in self.V
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test

        Note: This checks AC2(b^u) from the updated HP definition.
        This checks that the outcome occurs even if any subset of Z is set to its values in the experiment
        AND any subset of W set to the given values W satisfy the outcome

        Returns: True/False if condition is met

        """
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
    def ac2_m_check(self,X,X_indices,outcome_function,outcome_val):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test

        Note: This checks and finds the satisfying sets for AC2(m) from modified HP definition.
        Specifically, ensures that there is a set W, as defined in func above,
        s.t. W is set to the same settings of self.V and there is a setting of X called x''
        s.t. that outcome given W and X = x' is the opposite of the outcome_val

        Returns: [W, W_indices, T/F]: Where W and W_indices is the list that satisfy the definition.
        T/F correspond to whether there is any W that satisfies AC2(m).

        """
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

    def ac2_check_given_Z_W(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val):
        ac2_a = self.ac2_a_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        ac2_b = self.ac2_b_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        if(ac2_a and ac2_b):
            return True
        else:
            return False

    def ac2_u_check_given_Z_W(self,Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val):
        ac2_a = self.ac2_a_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        ac2_b = self.ac2_b_u_check(Z,Z_indices,X,X_indices,W,W_indices,outcome_function,outcome_val)
        if(ac2_a and ac2_b):
            return True
        else:
            return False


    def Z_and_W_search(self,X,X_indices,outcome_function,outcome_val,defn,num_outcomes,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Defn: {0,1} integer referring to original or updated definition respectively.
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test
                Fast: T/F corresponding whether to return after finding the first Z,W to satisfy

        Note: This function finds the Z and W that satisfy AC2.

        Returns: A dictionary of satisfying W. Keys: num_changes from original values of V | Value: W, W_indices
        If dictionary is empty => AC2 Failed

        """
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


    def ac3_check(self,X,X_indices,outcome_function,outcome_val):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_function: Function that determines outcome given self.U and self.V
                Outcome_val: Desired outcome val/setting to test

        Note: This function checks the AC3 condition, which checks that X is minimal by iterating over all
        subsets of X.

        Returns: T/F corresponding to whether AC3 is satisfied

        """
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

    def paths_to_pivotality(self,X,X_indices,outcome_val,outcome_func,num_outcomes=1):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)

        Note: This function is part of the Zultan model of responsibility (https://www.ncbi.nlm.nih.gov/pubmed/23855451)
        It returns all the cases where X is pivotal, i.e. states where changing the value of X would change the outcome.

        Returns: A list of tuples
                tup[0] = number of variables changed
                tup[1] = endogenous state for 1 valid path

        """
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

    def wrong_check(self,ac_1,ac_2,ac_3):
        """
        Prints which condition failed
        """
        if (not ac_1):
            print("(False b/c of AC1)")
        if(not ac_2):
            print("(False b/c of AC2)")
        if(not ac_3):
            print("(False b/c of AC3)")

    def causality_check(self,X,X_indices,outcome_val,outcome_func,num_outcomes,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This checks the original definition of causality from Halpern and Pearl.

        Returns: T/F if X is a causal variable/vector of variables

        """
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        dict_to_bool = lambda x: False if len(x) == 0 else True
        ac_2 = dict_to_bool(self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,0,num_outcomes,fast))
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    #Returns true if X satisfied the modified defn. (using ac1,ac2(a),ac2(b^u), ac3 - Halpern & Pearl 2005)
    def updated_causality_check(self,X,X_indices,outcome_val,outcome_func,num_outcomes,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This checks the updated definition of causality (AC1, AC2(a), AC2(b^u), AC3) from Halpern and Pearl.

        Returns: T/F if X is a causal variable/vector of variables

        """
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        dict_to_bool = lambda x: False if len(x) == 0 else True
        ac_2 = dict_to_bool(self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,1,num_outcomes,fast))
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    def modified_causality_check(self,X,X_indices,outcome_val,outcome_func):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V

        Note: This checks the modified definition of causality (AC1, AC2_m, AC3) from Halpern and Pearl.

        Returns: T/F if X is a causal variable/vector of variables

        """
        ac_1 = self.ac1_check(X,X_indices,outcome_func,outcome_val)
        ac_2 = self.ac2_m_check(X,X_indices,outcome_func,outcome_val)[2]
        ac_3 = self.ac3_check(X,X_indices,outcome_func,outcome_val)
        self.wrong_check(ac_1,ac_2,ac_3)
        return ac_1 and ac_2 and ac_3

    def responsibility(self,X,X_indices,outcome_val,outcome_func,num_outcomes,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This function returns the Chockler & Halpern (2004) definition of responsibility

        Returns: {1 / (minimum number of changes to make X pivotal + 1), if X is causal
                  0, o.w. }

        """
        if(self.causality_check(X,X_indices,outcome_val,outcome_func,num_outcomes,fast)):
            not_X = [not i for i in X]
            if(outcome_func(self.exogenous_var,V_given_set(self.endogenous_var,not_X,X_indices))!=outcome_val):
                return 1
            else:
                W_dct = self.Z_and_W_search(X,X_indices,outcome_func,outcome_val,0,num_outcomes,False)
                min_num_changes = min(W_dct.keys())
                return 1 / (min_num_changes + 1)
        else:
            return 0


    def influence(self,X,X_indices,num_outcome_var,outcome_func):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcome_var: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_function: Function that determines outcome given self.U and self.V

        Note: This function returns the boolean influence of a variable X
        More Info: https://en.wikipedia.org/wiki/Analysis_of_Boolean_functions#Influence

        Returns: Boolean Influence of X

        """
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

    def new_responsibility(self,X,X_indices,num_outcomes,outcome_val, outcome_func):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V

        Note: This function test a prospective definition of responsibility that we are workshopping.

        Returns: The boolean influence of X conditioned on the event that outcome_function(X) occurs

        """
        U = self.exogenous_var
        end_idx = (-1)*num_outcomes
        V = self.endogenous_var
        V_indices_no_X = [i for i in range(len(V[:end_idx])) if i not in X_indices]
        subsets_of_V = subsets_finder(V_indices_no_X)
        if(len(subsets_of_V) == 0):
            print("V")
            return 1
        outcome_change_ct = 0
        for subset in subsets_of_V:
            non_X_var_values = [True if i in subset else False for i in V_indices_no_X]
            V_scenario = V_given_set(V,non_X_var_values,V_indices_no_X)
            V_pos_X_scenario = V_given_set(V_scenario,X,X_indices)
            V_neg_X_scenario = V_given_set(V_scenario,[not i for i in X], X_indices)
            if(outcome_func(U,V_pos_X_scenario) and X == [False]):
                outcome_change_ct += 1
            elif(outcome_func(U,V_neg_X_scenario) and X == [True]):
                outcome_change_ct += 1
        return outcome_change_ct / len(subsets_of_V)


    def adj_responsibility(self,X,X_indices,num_outcomes,outcome_val,outcome_func,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcomes: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This function tests a potential new definition of responsibilty.

        Returns: Boolean_influence(X) * Chockler_Halpern_responsibility(X)

        """
        inf = self.influence(X,X_indices,num_outcomes,outcome_func)
        res = self.responsibility(X,X_indices,outcome_val,outcome_func,num_outcomes,fast)
        return inf*res

    def adj_responsibility_2(self,X,X_indices,num_outcome_var,outcome_val,outcome_func,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcome_var: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This function tests a potential new definition of responsibilty.

        Returns: sum of influence of each var in W that change (if X is causal, o.w. = 0 )

        """
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

    def adj_responsibility_2m(self,X,X_indices,num_outcome_var,outcome_val,outcome_func,fast):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcome_var: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Fast: T/F whether function should return after finding first satifying assignment

        Note: This function returns adj_2 and zultan responsibilty at once (to save time/compute)

        Returns: (a,z) where a = adj_2 definition and z = zultan_responsibility

        """
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


    def zultan_responsibility(self,X,X_indices,outcome_val,outcome_func,num_outcome_var):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Outcome_val: Desired outcome val/setting to test
                Outcome_function: Function that determines outcome given self.U and self.V
                Num_outcome_var: Number of outcomes (sometimes there are mutiple outcomes)

        Note: This function returns the definition of responsibility from Zultan's paper
        Link: https://www.ncbi.nlm.nih.gov/pubmed/23855451

        Returns: { 1, X is pivotal
                    1 / (sum of 1 / c_i where c_i is the number of changes to W in paths to pivotality), X is causality
                    0, o.w. }

        """
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

    def mc_inf_sample(self,X,X_indices,num_outcome_var,outcome_func,n):
        """
        Input:  X: List of values of desired variables to test
                X_indices: List containing indices of X variables in self.V
                Num_outcome_var: Number of outcomes (sometimes there are mutiple outcomes)
                Outcome_function: Function that determines outcome given self.U and self.V
                n: number of samples to take

        Note: This function estimates Boolean influence through Monte Carlo sampling

        Returns: An estimate of Boolean influence over n samples

        """
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
