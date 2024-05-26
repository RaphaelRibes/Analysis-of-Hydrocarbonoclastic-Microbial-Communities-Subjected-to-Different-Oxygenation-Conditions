import numpy as np
from data.explorer import dt_loader, true_loader
from tqdm.auto import tqdm
import json
import os
from scipy.special import comb
import matplotlib.pyplot as plt
import seaborn as sbn


def calculate_pj(N, N1, N2, j):
    """
    Calculate the probability pj that two species co-occur at exactly j sites.

    Parameters:
    N (int): Total number of sites where either species could occur.
    N1 (int): Number of sites occupied by species 1.
    N2 (int): Number of sites occupied by species 2.
    j (int): Number of sites where both species co-occur.

    Returns:
    float: Probability that species 1 and 2 co-occur at exactly j sites.
    """
    if j < max(0, N1 + N2 - N) or j > min(N1, N2):
        return 0

    numerator = (comb(N, j) * comb(N - j, N2 - j) * comb(N - N2, N1 - j))
    denominator = (comb(N, N1) * comb(N, N2))

    return numerator / denominator


def calculate_probabilities(N, N1, N2, Qobs):
    """
    Calculate Plt, Pgt, and Pet.

    Parameters:
    N (int): Total number of sites where either species could occur.
    N1 (int): Number of sites occupied by species 1.
    N2 (int): Number of sites occupied by species 2.
    Qobs (int): Observed number of sites having both species.

    Returns:
    tuple: Probabilities Plt, Pgt, and Pet.
    """
    p_lt = sum(calculate_pj(N, N1, N2, j) for j in range(0, Qobs))
    p_et = calculate_pj(N, N1, N2, Qobs)
    p_gt = sum(calculate_pj(N, N1, N2, j) for j in range(Qobs + 1, N + 1))

    return p_lt, p_et, p_gt


def coocurence_matrix(N: list, q_obs: np.ndarray, n_site: int):
    """based on the probabilistic model for analysing species co-occurrence
    Joseph A. Veech 10.1111/j.1466-8238.2012.00789.x

    N: number of occurrences for each group
    q_obs: the matrix of the observed number of sites having both species"""
    positions = np.zeros((len(N), len(N)))
    for n, N1 in tqdm(enumerate(N), total=len(N), desc='Calculating positions'):
        for m, N2 in enumerate(N):
            p_lt, p_et, p_gt = calculate_probabilities(n_site, N1, N2, int(q_obs[n, m]))

            if p_lt + p_et < 0.05:
                positions[n, m] = -1
            elif p_gt + p_et < 0.05:
                positions[n, m] = 1
            else:
                positions[n, m] = 0

    return positions


def presence_matrix(data, names, threshold=0.01):
    """
    This function generates a presence matrix for a given dataset, based on a specified threshold.
    The presence matrix is a binary matrix that indicates the presence (1) or absence (0) of a species in a sample.

    Parameters:
    data (DataFrame): The dataset that contains the information about the samples and species.
    names (dict): A dictionary that maps the species names to their corresponding taxonomy.
    threshold (float): A value that determines the cutoff for considering a species as present in a sample.

    Returns:
    tuple: The final presence matrix and the corresponding species names as numpy arrays.
    """

    # Create a reversed dictionary from the names dictionary
    r_names = {v: k for k, v in names.items()}

    # Sort the keys of r_names and store them in raw_names
    raw_names = np.sort([n for n in r_names.keys()])

    # Get the unique samples from the data
    samples = data['Sample'].unique()

    # Create the presence matrix
    matrix = [[int(data[(data['Taxonomy'] == r_names[name]) & (data['Sample'] == sample) & (data['Condition'] == 'OSCIL')]['Abundance'].mean() > threshold) for name in raw_names] for sample in tqdm(samples, desc='Creating the presence matrix')]

    # Transpose the matrix to align with the conventional representation of presence-absence matrices where rows represent species and columns represent samples
    matrix = np.transpose(matrix)

    # Filter out the species that are not present in any sample
    matrice_f, names_f = zip(*[(row, raw_names[i]) for i, row in enumerate(matrix) if np.sum(row) > 0])

    return np.array(matrice_f), np.array(names_f)


def observed_coocurence(matrix: np.ndarray) -> np.ndarray:
    """
    This function calculates the observed co-occurrence of species in a given presence-absence matrix.

    Parameters:
    matrix (np.ndarray): A binary matrix where each row represents a species and each column represents a sample.
                         A value of 1 indicates the presence of the species in the sample, and 0 indicates its absence.

    Returns:
    np.ndarray: A matrix where each element [i, j] represents the number of samples where both species i and j are present.
    """
    q_obs = np.zeros_like(matrix)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            q_obs[i, j] = np.sum((matrix[i] + matrix[j]) == 2)
    return q_obs


def load_data():
    """
    This function checks if the presence-absence matrix and species names are already saved in 'fig3.json' and 'names.json' files.
    If the files exist, it loads the matrix and names from the files.
    If the files do not exist, it calls the `true_loader` function to load the dataset and species names,
    generates the presence-absence matrix using the `presence_matrix` function, and saves the matrix and names to 'fig3.json' and 'names.json' files.

    Returns:
    tuple: The presence-absence matrix and the corresponding species names as numpy arrays.
    """
    if os.path.exists('fig3.json') and os.path.exists('names.json'):
        with open('fig3.json', 'r') as f:
            matrix = np.array(json.load(f))
        with open('names.json', 'r') as f:
            names = np.array(json.load(f))
    else:
        dt, names = true_loader()
        matrix, names = presence_matrix(dt, names)
        with open('fig3.json', 'w') as f:
            json.dump(matrix.tolist(), f)
        with open('names.json', 'w') as f:
            json.dump(names.tolist(), f)
    return matrix, names


def calculate_positions(matrix, dt):
    """
    This function calculates the positions of species co-occurrence based on a given presence-absence matrix and a dataset.
    It first calculates the total number of occurrences for each species and the unique number of samples under the 'OSCIL' condition.
    If the positions have been previously calculated and saved in 'fig3_positions.json', it loads the positions from the file.
    Otherwise, it calculates the positions using the `coocurence_matrix` function and saves them to 'fig3_positions.json'.

    Parameters:
    matrix (np.ndarray): The presence-absence matrix where each row represents a species and each column represents a sample.
    Dt (DataFrame): The dataset that contains the information about the samples and species.

    Returns:
    np.ndarray: The positions of species co-occurrence.
    """
    N = np.sum(matrix, axis=1)
    n_site = dt[dt['Condition'] == 'OSCIL']['Sample'].nunique()
    if os.path.exists('fig3_positions.json'):
        with open('fig3_positions.json', 'r') as f:
            positions = json.load(f)
    else:
        positions = coocurence_matrix(N, observed_coocurence(matrix), n_site)
        with open('fig3_positions.json', 'w') as f:
            json.dump(positions.tolist(), f)
    return np.array(positions)


def calculate_bilan(positions):
    """
    This function calculates the count of positive, negative, and random co-occurrences in the given positions matrix.
    A positive co-occurrence is represented by 1, a negative co-occurrence is represented by -1, and a random co-occurrence is represented by 0 in the positions matrix.

    Parameters:
    positions (np.ndarray): The positions matrix where each element [i, j] represents the co-occurrence of species i and j.

    Returns:
    dict: A dictionary with keys 'positive', 'negative', and 'random', and values representing the count of positive, negative, and random co-occurrences, respectively.
    """
    bilan = {'positive': 0, 'negative': 0, 'random': 0}
    for p in range(len(positions)):
        for q in range(len(positions)):
            if positions[p, q] == 1:
                bilan['positive'] += 1
            elif positions[p, q] == -1:
                bilan['negative'] += 1
            else:
                bilan['random'] += 1
    return bilan


matrix, names = load_data()
dt, _ = true_loader()
positions = calculate_positions(matrix, dt)

# Create a mask for the upper triangle of the positions matrix.
# np.triu creates an array with the same shape as positions, filled with ones in the upper triangle and zeros elsewhere.
# The second argument, k=0, specifies that the diagonal should be included in the upper triangle.
mask = np.triu(np.ones_like(positions), k=0)

# Apply the mask to the positions matrix.
# np.ma.array creates a masked array, which is an array that has a separate Boolean mask that indicates where data is "missing" or "invalid".
# In this case, the mask indicates the upper triangle of the positions matrix, so the resulting masked array, fixed, has the same values as positions in the lower triangle and masked values in the upper triangle.
fixed = np.ma.array(positions, mask=mask)

bilan = calculate_bilan(fixed)

print(bilan)
print(f'Percentage of non random: {100 * (bilan["positive"] + bilan["negative"]) / (bilan["positive"] + bilan["negative"] + bilan["random"]):.2f}%')


def plot_figure(fixed, names):
    """
    This function plots a figure representing the co-occurrence of species based on a given positions matrix and species names.
    The figure is saved in three formats: PNG, SVG, and TIF.

    Parameters:
    fixed (np.ndarray): The positions matrix where each element [i, j] represents the co-occurrence of species i and j.
    names (np.ndarray): The names of the species.

    Returns:
    None
    """
    # Create a new figure and axes with a specified size
    fig, ax = plt.subplots(figsize=(20, 20))

    # Display the positions matrix as an image
    ax.imshow(fixed, vmin=-1, vmax=1, cmap=sbn.diverging_palette(240, 10, as_cmap=True))

    # Remove borders
    for edge, spine in ax.spines.items():  # edging the spine
        spine.set_visible(False)

    # Remove the y-axis and x-axis ticks
    ax.set_yticks([])
    ax.set_xticks([])

    # Load the classifications from the 'metabolisms.json' file
    with open('metabolisms.json', 'r') as f:
        classifications = json.load(f)

    # Define the metabolisms and their corresponding colors
    metabolisms = ['Dégradation d\'hydrocarbure', 'Sulfo-oxidation', 'Fixation de l\'azote', 'Sulfato-réduction', 'Dénitrification']
    col = ['black', 'pink', 'green', 'red', 'blue']

    # Initialize a list to keep track of the first occurrence of each metabolism
    first = [True, True, True, True, True]

    # Loop over each species
    for cla in names:
        # Format the species name
        f_cla = cla.replace('_', ' ').replace('[Desulfobacterium] catecholicum group', 'Desulfobacterium catecholicum')

        # Add the species name to the plot
        ax.text(names.tolist().index(cla), names.tolist().index(cla), f_cla,
                ha='left', va='center', color='black', fontsize=13,
                fontstyle='italic' if 'Unknown ' not in cla else 'normal')

        # Loop over each metabolism
        for decalage, met in enumerate(metabolisms):
            # If the species has the current metabolism, add a dot to the plot
            if classifications[f_cla][decalage] > 0:
                ax.plot(-1 - 1 * decalage, names.tolist().index(cla), 'o',
                        color=col[decalage], label=met.capitalize() if first[decalage] else '', markersize=10,
                        markeredgecolor='black' if col[decalage] == 'pink' else col[decalage])

                # If this is the first occurrence of the current metabolism, set the corresponding element in the 'first' list to False
                if first[decalage]:
                    first[decalage] = False

    # Add a legend to the plot
    plt.legend(loc='upper center', ncol=1, fontsize=13, title='Metabolismes', title_fontsize=20, shadow=True, fancybox=True)

    # Adjust the layout of the plot
    fig.tight_layout()

    # Save the figure in PNG, SVG, and TIF formats
    fig.savefig('fig3.png', dpi=300, transparent=True)
    fig.savefig('fig3.svg', dpi=300)
    fig.savefig('fig3.tif', transparent=True)

    # Display the figure
    plt.show()

    # Close the figure
    plt.close(fig)


plot_figure(fixed, names)

with open('metabolisms.json', 'r') as f:
    classifications = json.load(f)

metabolisms = ['Dégradation d\'hydrocarbure', 'Sulfo-oxidation', 'Fixation de l\'azote', 'Sulfato-réduction',
               'Dénitrification']
colors = ['black', 'pink', 'green', 'red', 'blue']

# Iterate over the range of metabolisms
for itera1 in tqdm(range(5), desc='Creating subplots'):
    for itera2 in range(5):
        # Initialize lists to store the positions, names, and indices of the species
        denit_vs_desulf = []
        new_names = []
        indices = []

        # Iterate over the species names
        for name in names:
            # Format the species name
            f_cla = name.replace('_', ' ').replace('[Desulfobacterium] catecholicum group', 'Desulfobacterium catecholicum')

            # If the species has either of the current metabolisms, add its position to denit_vs_desulf, its name to new_names, and 1 to indices
            if classifications[f_cla][itera1] > 0 or classifications[f_cla][itera2] > 0:
                denit_vs_desulf.append(positions[names.tolist().index(name)])
                new_names.append(name)
                indices.append(1)
            else:
                # If the species does not have either of the current metabolisms, add 0 to indices
                indices.append(0)

        # Initialize a list to store the final positions
        final = []

        # Iterate over the positions in denit_vs_desulf
        for i, row in enumerate(denit_vs_desulf):
            # Add a new list to final
            final.append([])

            # Iterate over the positions in the current row
            for j, col in enumerate(row):
                # If the species at the current position has either of the current metabolisms, add its position to the last list in final
                if indices[j] == 1:
                    final[-1].append(col)

        # Convert final to a numpy array and mask the upper triangle
        denit_vs_desulf = np.array(final)
        mask = np.triu(np.ones_like(final), k=0)
        denit_vs_desulf = np.ma.array(denit_vs_desulf, mask=mask)

        # Create a new figure and axes with a specified size
        fig, ax = plt.subplots(figsize=(10, 10))

        # Display the positions matrix as an image
        im = ax.imshow(denit_vs_desulf, vmin=-1, vmax=1, cmap=sbn.diverging_palette(240, 10, as_cmap=True))

        # Remove the y-axis and x-axis ticks
        ax.set_yticks(np.arange(len(new_names)))
        ax.set_xticks(np.arange(len(new_names)))

        # Remove the plot borders
        for edge, spine in ax.spines.items():
            spine.set_visible(False)

        # Initialize a list to keep track of the first occurrence of each metabolism
        verif = [True, True, True, True, True]

        # Iterate over the new species names
        for name in new_names:
            # Format the species name
            f_cla = name.replace('_', ' ').replace('[Desulfobacterium] catecholicum group', 'Desulfobacterium catecholicum')

            # Add the species name to the plot
            ax.text(new_names.index(name), new_names.index(name), f_cla, ha='left', va='center', color='black')

            # Iterate over the metabolisms
            for i in range(5):
                # If the species has the current metabolism, add a dot to the plot
                if classifications[f_cla][i] > 0:
                    ax.plot(-1 - i, new_names.index(name), 'o', color=colors[i],
                            label=metabolisms[i] if verif[i] else '',
                            markeredgecolor='black' if colors[i] == 'pink' else colors[i])
                    # If this is the first occurrence of the current metabolism, set the corresponding element in the 'verif' list to False
                    if verif[i]:
                        verif[i] = False

        # Remove the y-axis and x-axis ticks
        ax.set_xticks([])
        ax.set_yticks([])

        # Add a legend to the plot
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

        # Save the figure in PNG, SVG, and TIF formats
        fig.savefig(f'fig3 {metabolisms[itera1]} vs {metabolisms[itera2]}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'fig3 {metabolisms[itera1]} vs {metabolisms[itera2]}.svg', dpi=300, bbox_inches='tight')
        fig.savefig(f'fig3 {metabolisms[itera1]} vs {metabolisms[itera2]}.tif', bbox_inches='tight')

        # Close the figure
        plt.close(fig)
