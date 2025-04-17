def lexoff(model,first_opt,coefficient_1,coefficient_2,criteria):
    # Including the versions I used when developing -->
    import pandas as pd #Version 2.2.3
    import cobra #Version 0.29.0
    from cobra.io import read_sbml_model
    import matplotlib.pyplot as plt #Version 3.10.0
    from cobra import Model as CobraModel
    import warnings
    warnings.filterwarnings('ignore', 
                           message='DataFrame is highly fragmented', 
                           category=pd.errors.PerformanceWarning)
    '''
    GOAL: Performs a lexicographic optimisation on a COBRA model, to predict the objective function(s) that best satisfy experimental criteria
    PROCESS: 
        1. First optimisation is assumed by the User and maximisation or minimisation is specified
        2. The highest priority is placed on the first optimisation, then this objective value is used as additional model constraints for a second optimisation
        3. Second optimisation can be either maximisation or minimisation of any exchange reaction in COBRA model
        4. Comparison of predicted fluxes to an experimental criteria, which has been specified by the User
    
    PARAMETERS:
        - 'model': A COBRA model, ideally omics constrained (e.g. transcriptomics or proteomics)
        - 'first_opt': A reaction ID that is a User assumption of the first objective function
        - 'coefficient_1': either 1 or -1 to specify the maximisation or minimisation of 'first_opt'
        - 'coefficient_2': either 1 or -1 to specify the maximisation or minimisation of every exchange reaction in turn, to be tested as the second optimisation
        - 'criteria': a qualitative experimental criteria, in the form of a dictionary of reaction IDs (keys) and a qualitative flux direction (values) of 1, 0 or -1 for flux in forward direction (i.e. production), no flux or flux in reverse direction (i.e. consumption)
    
    Returns:
        - Dataframe of predicted absolute fluxes, qualitative fluxes and proportion accuracy to experimental criteria for each lexicographic optimisation. Column headers are reaction ID of second optimisation
        - Dictionary: Accuracy scores to criteria (keys) and number of reactions that achieved this score (values)
        - Message reporting maximum accuracy score achieved, the number of reactions that achieved this and an example of the first few reactions IDs with maximum accuracy
        - Bar plot with the distribution of accuracies achieved by optimisations
    '''
    
    # Input validation
    if not isinstance(model, CobraModel):
        raise TypeError("Check inputs: Model must be a cobra.Model.Model instance")
    if not isinstance(criteria, dict):
        raise TypeError("Check inputs: Criteria should be in dictionary format with keys as reaction IDs and values as 1, 0 or -1")
    if coefficient_1 not in (-1, 1):
        raise ValueError("Check inputs: Coefficient 1 must be either 1 or -1")
    if coefficient_2 not in (-1, 1):
        raise ValueError("Check inputs: Coefficient 2 must be either 1 or -1")

    # Ensure exchange reaction boundaries reflect reversibility if the experimental criteria specifies a metabolite uptake
    for k,v in criteria.items():
        if criteria[k] == -1: # If you have specified that there is uptake of this metabolite
            if model.reactions.get_by_id(k).reversibility == False: # But if the inbuilt bounds don't allow uptake of this metabolite...
                model.reactions.get_by_id(k).lower_bound = -1000 # Update the model bounds to allow for metabolite uptake
                print('reaction bounds for:',k,'have been updated from 0,+1000 to -1000,+1000 to allow metabolite uptake')
            else:
                continue # If reversibility is already allowed
        else:
            continue # If criteria specifies 0 flux (0) or production (+1)

    # Print optimisation instructions to confirm to User what they have coded
    if coefficient_1 == 1:
        first_goal = 'maximise'
    if coefficient_1 == -1:
        first_goal = 'minimise'
    if coefficient_2 == 1:
        second_goal = 'maximisation'
    if coefficient_2 == -1:
        second_goal = 'minimisation'
    print(f"Optimisation: User would like to {first_goal} {first_opt}, then explore the {second_goal} of exchange reactions as the second optimisation")
    
    # Set up empty dataframe to append with optimisation results later on
    objectives_df = pd.DataFrame()
    first_col = []
    for k in criteria.keys(): # Add index of reaction IDs for criteria (to insert absolute fluxes)
        first_col.append(k)
    for k in criteria.keys(): # Add index of reaction IDs again (to insert qualitative fluxes)
        first_col.append(k)
    first_col.append('accuracy') # Add index for accuracy to experimental data
    objectives_df['criteria_reactions'] = first_col
    objectives_df = objectives_df.set_index('criteria_reactions')
    criteria_keys = []
    for c in criteria.keys():
        criteria_keys.append(c)
    
    # Make test model and run optimisations
    model_test = model.copy()
    model = model_test
    columns_to_add = {}
    for e in model.exchanges: # Only exchange reactions will be explored within the objective
        with model:
            # First optimisation
            model.objective = {model.reactions.get_by_id(first_opt):coefficient_1}
            objective_sol = model.slim_optimize()
            model.reactions.get_by_id(first_opt).bounds = (objective_sol,objective_sol) # Set objective value as additional constraints
            # Second optimisation
            model.objective = {model.reactions.get_by_id(e.id):coefficient_2}
            solution = model.optimize()
            
            # Fill in dataframe column with absolute solutions, qualitative solutions and accuracy
            # Absolute solutions...
            criteria_solutions = []
            for c in criteria.keys():
                criteria_solutions.append(solution.fluxes[c])
            criteria_col = []
            for c in criteria_solutions:
                criteria_col.append(c)
            qual_values = []
            # Qualitative solution
            for c in criteria_solutions:
                if c>0:
                    criteria_col.append(1)
                    qual_values.append(1)
                if c<0:
                    criteria_col.append(-1)
                    qual_values.append(-1)
                if c==0:
                    criteria_col.append(0)
                    qual_values.append(0)
            true_values = []
            for q in criteria.values():
                true_values.append(q)
            test_values = qual_values
            # Compute and insert accuracy scores
            accuracy = sum(int(a) == int(b) for a, b in zip(true_values, test_values))
            criteria_col.append(accuracy/len(criteria.keys()))
            objectives_df[e.id] = criteria_col
            
    # Sort columns according to accuracy score
    sorted_columns = objectives_df.iloc[-1].sort_values(ascending=True).index
    sorted_df = objectives_df[sorted_columns]
    
    # Create dictionary from dataframe of optimsiation solutions
    accuracies = (sorted_df).iloc[-1]
    unique_accuracies = [float(val) for val in accuracies.unique()]
    accuracy_dict = {val: list(accuracies[accuracies == val].index) for val in unique_accuracies}
    
    # Print message summarising the optimisation results
    max_accuracy = max(unique_accuracies)
    count = len(accuracy_dict[max_accuracy])
    first_few = accuracy_dict[max_accuracy][:3]
    message = f"The maximum accuracy is {max_accuracy} and this has been achieved by {count} unique reactions for second optimisation, including the {second_goal} of: {first_few}"
    print(message)
    
    # Bar plot for distribution of accuracies achieved by optimisations
    x_axis = unique_accuracies
    y_axis = [len(rxns) for rxns in accuracy_dict.values()]
    
    for acc, rxns in accuracy_dict.items():
        print(f"{acc:.2f} accuracy: {len(rxns)} reactions")
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(x_axis, y_axis, color='black', edgecolor='black', width=0.05)
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
    ax.set_xlabel('Accuracy')
    ax.set_ylabel('Frequency of Reactions')
    ax.set_title('Distribution of Optimisation Accuracies')
    ax.set_xlim(-0.05, 1.05)
    ax.set_xticks([0] + x_axis + [1])
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Return variables
    return(sorted_df,accuracy_dict,message,fig)