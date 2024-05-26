# Import necessary libraries and modules
from data.explorer import true_loader
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import braycurtis
import matplotlib.pyplot as plt
import numpy as np
import json
import os
from skbio.stats.distance import permanova
from skbio import DistanceMatrix
from tqdm.auto import tqdm

# Load the data
dt, names = true_loader()
arr = []

# Check if the preprocessed data exists
if os.path.exists('fig 2.a wt pretreated.json'):
    # If it does, load the data
    with open('fig 2.a wt pretreated.json', 'r') as f:
        arr = json.load(f)
else:
    # If it doesn't, preprocess the data
    for comm in tqdm(['A', 'B', 'C']):
        temp = [[], [], []]
        for taxo in dt['Taxonomy'].unique():
            for n, sample in enumerate(dt[(dt['Taxonomy'] == taxo)
                                          & (dt['Time'] == 5.)
                                          & (dt['Community'] == comm)
                                          & (dt['Condition'] == "OSCIL")]['Abundance'].tolist()):
                temp[n].append(sample)
        for n in range(3):
            arr.append(temp[n])

    # Save the preprocessed data
    with open('fig 2.a wt pretreated.json', 'w') as f:
        json.dump(arr, f)

# Convert the data to a numpy array
arr = np.array(arr)

# Calculate the Bray-Curtis distance
size = len(arr)
diss_matrix = np.zeros((size, size))  # Initialize the dissimilarity matrix
for i in range(size):
    for j in range(size):
        # Calculate the Bray-Curtis distance between each pair of points
        diss_matrix[i, j] = braycurtis(arr[i], arr[j])

# Perform the PCoA analysis
pcoa_results = pcoa(diss_matrix)

# Get the coordinates of the points in the plane
plot_matrix = []
third_dim = []
for (x, y, z) in zip(pcoa_results.samples['PC1'], pcoa_results.samples['PC2'], pcoa_results.samples['PC3']):
    plot_matrix.append([x, y, z])
    third_dim.append(z)

third_dim = np.array(third_dim)
pos_lim = (min(third_dim), max(third_dim))

# Create a new figure and axes with 3D projection
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Define the positions and colors for each community
pos = [
    plot_matrix[:3],
    plot_matrix[3:6],
    plot_matrix[6:],
]
colors = ['r', 'b', 'g']

# Plot the data
for com, matrix in enumerate(pos):
    for j, p in enumerate(matrix):
        # Plot the points
        ax.scatter(p[0], p[1], p[2], color=colors[com], marker='o',
                   label=f'Community {['A', 'B', 'C'][com]}' if j == 0 else None)
        # Draw a line from the point to the bottom of the plot
        ax.plot([p[0], p[0]], [p[1], p[1]], [p[2], min(third_dim)], color=colors[com],
                linestyle='--', alpha=0.5, lw=0.75)
        # Put a cross at the bottom of the plot
        ax.scatter(p[0], p[1], min(third_dim), color=colors[com], marker='x', alpha=0.5, lw=0.75)

# Set the labels for the axes
percentages = pcoa_results.proportion_explained
ax.set_xlabel(f'PC1: {percentages[0] * 100:.2f}%')
ax.set_ylabel(f'PC2: {percentages[1] * 100:.2f}%')
ax.set_zlabel(f'PC3: {percentages[2] * 100:.2f}%')

# Remove the top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Set the view angle
ax.view_init(45, 45, 0)

# Show the plot
plt.show()

# Save the figure in different formats
fig.savefig('fig 2.a.png', dpi=300, transparent=True)
fig.savefig('fig 2.a.svg', transparent=True)
fig.savefig('fig 2.a.tif', dpi=300, transparent=True)

# Close the figure
plt.close(fig)

# Create a new figure and axes
fig, ax = plt.subplots()

# Plot the explained variance
ax.plot(range(len(pcoa_results.proportion_explained[:5])), [p*100 for p in pcoa_results.proportion_explained[:5]])
ax.set_ylabel('Explained Variance (%)')
ax.set_xticks(range(5))
ax.set_xticklabels(['PC1', 'PC2', 'PC3', 'PC4', 'PC5'], rotation=45)

# Color the area under the first three points in the y-axis
ax.fill_between(range(3), 0, [p*100 for p in pcoa_results.proportion_explained[:3]], color='gray', alpha=0.25)

# Annotate in the middle
ax.annotate(f'Explain {round(sum([p*100 for p in pcoa_results.proportion_explained[:3]]), 2)}%',
            (1, sum([p*100 for p in pcoa_results.proportion_explained[:3]])/6), ha='center')

# Adjust the layout
fig.tight_layout()

# Show the plot
plt.show()

# Save the figure in different formats
fig.savefig('fig 2.a explained variance.png', dpi=300, transparent=True)
fig.savefig('fig 2.a explained variance.svg', transparent=True)
fig.savefig('fig 2.a explained variance.tif', dpi=300, transparent=True)

# Close the figure
plt.close(fig)

# Calculate the global permanova
distances = DistanceMatrix(diss_matrix)
results_same = permanova(distances, grouping=['A']*3 + ['B']*3 + ['C']*3)

# Print the results
print(results_same)
