def optimize_objective_from_qualitative_criteria(
    model,
    objective_reactions,
    qualitative_constraints,
    bounds=None,
    max_iter=1000,
    maxfun=1000
):
    #import relevant packages
    import numpy as np
    import pandas as pd
    import cobra
    from cobra.io import read_sbml_model
    from scipy.optimize import dual_annealing
    from cobra import Model as CobraModel
    import warnings
    warnings.filterwarnings('ignore', 
                           message='DataFrame is highly fragmented', 
                           category=pd.errors.PerformanceWarning)
    """
    Optimizes objective coefficients for a COBRA model to best match qualitative flux constraints.
    Logs each function evaluation's normalized coefficients, accuracy, and flux values.

    Parameters:
    - model: COBRApy model
    - objective_reactions: list of reaction IDs to include in the objective
    - qualitative_constraints: dict of {reaction_id: expected_direction}, where direction ∈ {-1, 1}
    - bounds: list of (min, max) tuples for each reaction coefficient (default: (-1, 1) for all)
    - max_iter: int, number of maximum iterations for dual annealing

    Returns:
    - optimal_objective: dict of {reaction_id: normalized objective coefficient}
    - accuracy: float (0–1), final qualitative match accuracy
    - log_df: pandas DataFrame with columns ['coefficients', 'accuracy', 'fluxes']
    - message summarising the global minimum (maximum accuracy to experimental) and the fitted objective (reaction IDs and coefficients)
    - rounded_fluxes_df which stores the rounded fluxes for the criteria reactions each iterations, so we can analyse convergence
    - fig of the accuracy over time, to give an idea of convergence
    """

    #Check the input data
    if not isinstance(model, CobraModel):
        raise TypeError("Check inputs: Model must be a cobra.Model.Model instance")
    if not isinstance(objective_reactions, list):
        raise TypeError("Check inputs: objective_reactions should be a list of reaction IDs in model")
    if not isinstance(qualitative_constraints, dict):
        raise TypeError("Check inputs: qualitative_constraints should be a dictionary of reactions IDs and expected direction of reaction, i.e. either 1 or -1, depending on whether this reaction is a production or consumption")
    invalid_values = {k: v for k, v in qualitative_constraints.items() if v not in {-1, 0, 1}}
    if invalid_values:
        raise ValueError(
            f"Invalid reaction direction values found: {invalid_values}. "
            "Each value in qualitative_constraints must be -1, 0, or 1."
        )

    for k,v in qualitative_constraints.items():
        if qualitative_constraints[k] == -1: # If you have specified that there is uptake of this metabolite
            if model.reactions.get_by_id(k).reversibility == False: # But if the inbuilt bounds don't allow uptake of this metabolite...
                model.reactions.get_by_id(k).lower_bound = -1000 # Update the model bounds to allow for metabolite uptake
                print('reaction bounds for:',k,'have been updated from 0,+1000 to -1000,+1000 to allow metabolite uptake')
            else:
                continue # If reversibility is already allowed
        else:
            continue # If criteria specifies 0 flux (0) or production (+1)

    print(f"Optimisation: User would like to fit an objective function including reactions: {objective_reactions} to predict a flux distribution best matching experimental data, measuring {len(qualitative_constraints)} reactions")
    
    if bounds is None:
        bounds = [(-1, 1)] * len(objective_reactions)

    selected_qualitative_reactions = list(qualitative_constraints.keys())
    results_log = []

    def evaluate_solution(c_raw):
        # Normalize coefficients so absolute values sum to 1
        if np.sum(np.abs(c_raw)) == 0:
            c = np.zeros_like(c_raw)
        else:
            c = c_raw / np.sum(np.abs(c_raw))

        # Apply to model
        for i, rxn_id in enumerate(objective_reactions):
            model.reactions.get_by_id(rxn_id).objective_coefficient = c[i]

        solution = model.optimize()

        if solution is None:
            accuracy = 0.0
            mismatch_count = len(selected_qualitative_reactions)
            flux_dict = {rxn_id: None for rxn_id in selected_qualitative_reactions}
        else:
            fluxes = np.array([model.solver.variables[rxn_id].primal for rxn_id in selected_qualitative_reactions])
            expected = np.array([qualitative_constraints[rxn_id] for rxn_id in selected_qualitative_reactions])
            rounded_fluxes = np.where(np.abs(fluxes) < 1e-6, 0, np.sign(fluxes).astype(int))
            mismatch_count = np.count_nonzero(rounded_fluxes != expected)
            accuracy = 1 - (mismatch_count / len(expected))
            flux_dict = {rxn_id: flux for rxn_id, flux in zip(selected_qualitative_reactions, fluxes)}

        # Log evaluation
        results_log.append({
            'coefficients': c.tolist(),
            'accuracy': accuracy,
            'fluxes': flux_dict
        })

        return mismatch_count

    # Run the optimization
    result = dual_annealing(evaluate_solution, bounds, maxiter=max_iter, maxfun=maxfun)

    # Final scaling of optimal result
    if np.sum(np.abs(result.x)) == 0:
        scaled_coeffs = np.zeros_like(result.x)
    else:
        scaled_coeffs = result.x / np.sum(np.abs(result.x))

    for i, rxn_id in enumerate(objective_reactions):
        model.reactions.get_by_id(rxn_id).objective_coefficient = scaled_coeffs[i]

    model.optimize()
    final_fluxes = np.array([model.solver.variables[rxn_id].primal for rxn_id in selected_qualitative_reactions])
    expected = np.array([qualitative_constraints[rxn_id] for rxn_id in selected_qualitative_reactions])
    rounded_fluxes = np.where(np.abs(final_fluxes) < 1e-6, 0, np.sign(final_fluxes).astype(int))
    accuracy = 1 - (np.count_nonzero(rounded_fluxes != expected) / len(expected))

    log_df = pd.DataFrame(results_log)

    sampled_df = log_df.iloc[::10]
    fig,ax = plt.subplots(figsize=(11,5))
    ax.plot(sampled_df.index, sampled_df['accuracy'], marker='o', linestyle='-', color='black')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Accuracy')
    ax.set_title('Accuracy Over Time (Every 50th Evaluation)')
    ax.grid(True)

    print("Optimal Coefficients:", dict(zip(objective_reactions,scaled_coeffs)))
    print(f"Final Accuracy: {accuracy:.2%}")

    return(dict(zip(objective_reactions,scaled_coeffs)),accuracy,log_df,fig)