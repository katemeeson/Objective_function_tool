def visualising_objective_proportions(function_result): #function_result: the result of optimize_objective_from_qualitative_criteria_agreement function
    result = function_result[0]

    abs_values = {k: abs(v) for k, v in result.items()}
    total = sum(abs_values.values())
    proportions = {k: v / total for k, v in abs_values.items()}
    colors = {
        True: 'tab:blue',   # positive (maximization)
        False: 'tab:orange' # negative (minimization)
    }
    labels = list(proportions.keys())
    sizes = list(proportions.values())
    signs = [result[k] > 0 for k in labels]
    title_parts = [f"{k} ({result[k]:.3f})" for k in labels]
    title_text = "Objective function:\n" + ", ".join(title_parts)
    fig, ax = plt.subplots(figsize=(12, 3))
    left = 0
    for i, size in enumerate(sizes):
        ax.barh(0, size, left=left, color=colors[signs[i]], edgecolor='black')
        if size > 0.05:
            xpos = left + size / 2
            label_text = f"{labels[i]}\n({result[labels[i]]:.3f})"
            ax.text(xpos, 0.15, label_text, ha='center', va='bottom', fontsize=10)
        left += size
    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel('Proportion of objective function')
    ax.set_title(title_text, fontsize=12)
    legend_elements = [
        Patch(facecolor='tab:blue', edgecolor='black', label='Maximization (positive)'),
        Patch(facecolor='tab:orange', edgecolor='black', label='Minimization (negative)')
    ]
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()