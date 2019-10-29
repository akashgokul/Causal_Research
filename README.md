# Causal Library
### This repository is a collection of functions defining causation of a variable and the responsibility of a variable within a given model.
#### The python file: causality_and_responsibility.py contains all these functions and has comments to understand each function. The iPython Notebook was simply used for testing purposes.

## Syntax
A model consists of the following elements:
 - A set U of exogenous variables (variables outside of the system)
   - In our case, this U will be a list where U[i] = the boolean value of the i-th exogenous variable
 - A set V of endogenous variables (variables within the system)
   - Similar to U, V will be a list s.t. V[i] = boolean value of the i-th exogenous variable
 - A set R, which tells the ranges of each variables
  - This is included for syntax purposes but our analysis is restricted to variables with range {0,1} i.e. False-True.
  - Likewise, R is a list in our code.
 - A function F, which takes in U,V,R:
  - F models the outcome/scenario and will return True if a given outcome occured (e.g. Canidate A won the election) and False otherwise.

A CausalModel is an instance of said class that has the attributes mentioned above.

## Causal Analysis:
There are 3 iterations on the definition of causality: Halpern-Pearl (2001), Halpern-Pearl (2005), Halpern-Pearl (2014)
To test if variable(s) X (a list) is a cause of some outcome you simply call causality_check, updated_causality_check, or modified_causality_check respectively. Note #1: These are instance methods.
Note #2: -X refers to a list containing the values for the variables of interest (e.g. [True])
-X_indices refer to a list of their indices in V. (e.g. [0])
-num_outcome_var: refers to the index of the outcome variable in V
-outcome_val = expected outcome value (True or False)
-fast: a boolean flag which (if True) tells the analysis to halt early (if the definition is satisfied) or (if False) search the whole space of paths to causality

##Responsibility Attribution
Most of our research concerned the responsibility attribution function and we played around with different definitions.
Within this section you will find the following functions:
- responsibility: This function computes the responsibility of variable X on outcome Y using Chockler & Halpern's definition of responsibility (2004)
-influence: This computes the boolean influence of a variable (doesn't require causality)
- zultan_responsibility: This function uses the definition of responsibility defined in the paper Multiple Counterfactual Pivotality model by Zultan, Gerstenberg, Lagnado (2012). Warning: the description of this model is somewhat ambiguous and I believe that's why many of further testing did not match the paper's results. Unfortunately, they did not code the model.

Note 1: the same variables -> meanings, from the Causal Analysis section, hold here.
