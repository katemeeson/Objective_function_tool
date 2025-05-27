# Objective_function_tool (repository developed April 2025)
#### Author of repository: Kate Meeson
#### Aim: To develop a (Python-based) tool for the prediction of the objective function for lexicographic optimisation that best satisfies a qualitative experimental criteria
# File descriptions
#### Status: In development. If you would like to test this function or be involved in its benchmarking, please get in touch to collaborate: kate.meeson@manchester.ac.uk
#### 1. LEXOFF_function.py is the file containing the LEXOFF function; import into Jupyter notebook as 'import LEXOFF_function as lx' then use as 'lx.lexoff()'
#### 2. LEXicographic_Objective_Function_Finder_(and_example).ipynb gives the source code for the LEXOFF function itself and a usage example for Chinese Hamster Ovary (CHO) cells; function can be copied onto User's notebook and applied to a different model and experimental criteria
#### 3. Lexicographic_constraints_early_model.ipynb, Lexicographic_constraints_late_model.ipynb and Lexicographic_constraints_stationary_model.ipynb are early LEXOFF developments, where an alternative control was explored: instead of comparison to experimental qualitative criteria to determine accuracy of objective function, a comparison to a flux-constrained GEM was used
#### 4. Wrapping_simulated_annealing_into_a_function.ipynb is early developments of a simulated annealing tool that fits an objective function to experimental data, in this notebook it has been applied to yeast using published data from PMID: 35145105
