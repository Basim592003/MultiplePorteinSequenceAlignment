import matplotlib.pyplot as plt
import seaborn as sns
import os
from Bio import Phylo
import plotly.express as px
import numpy as np

def calculate_conservation_score(alignment):
    """
    Calculate conservation score at each position of the alignment with improved scoring.
    """
    alignment_array = np.array([list(rec) for rec in alignment])
    n_positions = alignment_array.shape[1]
    conservation_scores = []

    # Define amino acid groups for similarity scoring
    aa_groups = {
        'hydrophobic': set('AILMFWYV'),
        'polar': set('STNQ'),
        'positive': set('RHK'),
        'negative': set('DE'),
        'special': set('CGP')
    }

    for i in range(n_positions):
        column = alignment_array[:, i]
        unique_aas, counts = np.unique(column, return_counts=True)
        
        # Calculate basic conservation
        most_common = np.max(counts)
        basic_score = most_common / len(column)
        
        # Add similarity bonus
        similarity_bonus = 0
        for group in aa_groups.values():
            group_count = sum(1 for aa in column if aa in group)
            if group_count > 1:
                similarity_bonus += 0.1 * (group_count / len(column))
        
        # Penalize gaps
        gap_penalty = 0.2 * (column == '-').sum() / len(column)
        
        # Combine scores
        final_score = min(1.0, basic_score + similarity_bonus - gap_penalty)
        conservation_scores.append(final_score)

    return np.array(conservation_scores)



def plot_guide_tree(guide_tree_file):
    if not os.path.exists(guide_tree_file):
        return None

    try:
        guide_tree_data = Phylo.read(guide_tree_file, "newick")
        
        # Create figure with custom background
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Set the background colors first
        ax.set_facecolor('#0F1117')  # Dark navy/black background
        fig.patch.set_facecolor('#0F1117')  # Dark navy/black background
        
        # Draw the tree
        Phylo.draw(guide_tree_data, 
                  axes=ax, 
                  do_show=False,
                  label_func=lambda x: x.name,
                  show_confidence=False)
        
        # Force all elements to be white
        for c in ax.get_children():
            try:
                c.set_color('#ffffff')  # Set everything to white
                if isinstance(c, plt.Line2D):
                    c.set_linewidth(2)  # Make lines thicker
                elif isinstance(c, plt.Text):
                    c.set_fontsize(17)  # Set text size
            except:
                continue
        
        # Remove axes
        ax.axis('off')
        
        # Add title with light color
        
        # Adjust layout
        plt.tight_layout()
        
        return fig
    except Exception as e:
        print(f"Error plotting guide tree: {str(e)}")
        return None
    
def plot_plotly_heatmap(conservation_scores):
    # Create a more detailed heatmap
    fig = px.imshow(
        [conservation_scores],
        color_continuous_scale=[
            [0, 'rgb(255,255,255)'],      # White for low conservation
            [0.3, 'rgb(166,206,227)'],    # Light blue
            [0.6, 'rgb(31,120,180)'],     # Medium blue
            [0.8, 'rgb(178,223,138)'],    # Light green
            [1.0, 'rgb(51,160,44)']       # Dark green for high conservation
        ],
        aspect='auto'
    )

    # Customize layout
    fig.update_layout(

        xaxis_title={
            'text': "Position in Alignment",
            'font': dict(size=16)
        },
        yaxis_showticklabels=False,
        # plot_bgcolor='white',
        width=1000,
        height=300,  # Increased height from 250 to 400
        # margin=dict(l=0, r=0, t=60, b=60),  # Reduced margins
        coloraxis_colorbar=dict(
            title="Conservation<br>Score",
            titleside="right",
            ticks="outside",
            tickmode="array",
            ticktext=["Low", "Medium", "High"],
            tickvals=[0.2, 0.5, 0.8],  # Added middle value to match three labels
            len=0.8,
            x=1.02
        )
    )

    # Add gridlines
    fig.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgrey',
        zeroline=False,
        dtick=10,
        showline=False,  # Remove x-axis line
        mirror=False     # Remove mirror effect
    )

    # Remove y-axis lines
    fig.update_yaxes(
        showline=False,
        mirror=False,
        showgrid=False
    )



    return fig

def format_newick_string(newick_str):
    """Format a Newick tree string to be more readable with indentation."""
    formatted = ""
    indent = 0
    
    for char in newick_str:
        if char == '(':
            formatted += '(\n' + ' ' * (indent + 2)
            indent += 2
        elif char == ',':
            formatted += ',\n' + ' ' * indent
        elif char == ')':
            indent -= 2
            formatted += '\n' + ' ' * indent + ')'
        else:
            formatted += char
            
    return formatted

