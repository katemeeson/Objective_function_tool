def plot_agreement_heatmaps(function_result, num_columns=100): #function_result: the result of optimize_objective_from_qualitative_criteria_agreement function
    """
    Plots heatmaps of agreement (binary 0/1) for the first and last `num_columns` iterations.

    Parameters:
        agreement_df (pd.DataFrame): DataFrame with reaction IDs as rows and iteration steps as columns.
        num_columns (int): Number of columns (iterations) to plot per subplot (default: 100).
    """
    agreement_df = function_result[4]
    num_cols = agreement_df.shape[1]
    cols_to_plot = min(num_columns, num_cols)

    first_chunk = agreement_df.iloc[:, :cols_to_plot]
    last_chunk = agreement_df.iloc[:, -cols_to_plot:]

    # Define custom colormap
    cmap = mcolors.ListedColormap(["#d73027", "#1a9850"])  # red for 0, green for 1
    bounds = [-0.5, 0.5, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Create subplots
    fig_height = max(6, len(agreement_df) * 0.4)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(cols_to_plot / 4, fig_height))

    # Plot first chunk
    sns.heatmap(
        first_chunk,
        ax=axes[0],
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.2,
        linecolor='white',
        square=False,
        xticklabels=max(1, cols_to_plot // 10),
        yticklabels=True
    )
    axes[0].set_title(f"First {cols_to_plot} Iterations")
    axes[0].set_xlabel("Iteration")
    axes[0].set_ylabel("Reaction ID")

    # Plot last chunk
    sns.heatmap(
        last_chunk,
        ax=axes[1],
        cmap=cmap,
        norm=norm,
        cbar=False,
        linewidths=0.2,
        linecolor='white',
        square=False,
        xticklabels=max(1, cols_to_plot // 10),
        yticklabels=True
    )
    axes[1].set_title(f"Last {cols_to_plot} Iterations")
    axes[1].set_xlabel("Iteration")
    axes[1].set_ylabel("Reaction ID")

    plt.tight_layout()
    plt.show()