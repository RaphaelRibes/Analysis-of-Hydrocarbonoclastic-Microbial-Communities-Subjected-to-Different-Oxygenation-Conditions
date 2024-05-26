from data.explorer import dt_loader
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import braycurtis
import matplotlib.pyplot as plt
import numpy as np
import json
import os
from tqdm.auto import tqdm
from PIL import Image

data = dt_loader()

times = {'ANOXIC': [0, 8.0, 10., 11., 15.],
         'OXIC': [5, 5.4, 7, 8.0, 10., 11., 15.],
         'OSCIL': [5, 5.4, 7, 8.0, 10., 11., 15.]}

arr = []
verif = []
if os.path.exists('pretreated_2.b.json'):
    with open('pretreated_2.b.json', 'r') as f:
        arr = json.load(f)
else:
    for cond in ['OXIC', 'ANOXIC', 'OSCIL']:
        for t in times[cond]:
            abundance = [[], [], []]
            samples = data[(data['Condition'] == cond) & (data['Time'] == t)]['Sample'].unique()
            samples.sort()
            for n, sample in enumerate(samples):
                abundance[n // 3].append(
                    data[(data['Condition'] == cond) & (data['Time'] == t) & (data['Sample'] == sample)][
                        'Abundance'].tolist())
                verif.append((sample, len(abundance[n // 3][-1])))
            if t == 0:
                for ab in abundance[0]:
                    arr.append(ab)

            else:
                for ab in abundance:
                    arr.append(np.mean(ab, axis=0).tolist())
    with open('pretreated_2.b.json', 'w') as f:
        json.dump(arr, f)

arr = np.array(arr)

# on calcule la distance de Bray-Curtis
size = len(arr)
diss_matrix = np.zeros((size, size))  # on initialise la matrice de dissimilarité

for i in range(size):
    for j in range(size):
        diss_matrix[i, j] = braycurtis(arr[i], arr[j])

pcoa_results = pcoa(diss_matrix)

plot_matrix = []
third_dim = []
for (x, y, z) in zip(pcoa_results.samples['PC1'], pcoa_results.samples['PC2'], pcoa_results.samples['PC3']):
    plot_matrix.append([x, y, z])
    third_dim.append(z)

third_dim = np.array(third_dim)
pos_lim = (min(third_dim), max(third_dim))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

pos = [[], [], []]
for n, point in enumerate(plot_matrix):
    pos[n % 3].append(point)

form = ['o', 's', '^']
colors = ['b', 'g', 'r']
for n, comm in enumerate(['B', 'C', 'A']):
    for i, p in enumerate(pos[n]):
        if i < 7:
            marker = form[0]
        elif 7 <= i < 14:
            marker = form[1]
        else:
            marker = form[2]
        ax.scatter(p[0], p[1], p[2], c=colors[n], marker=marker, label=comm if i == 0 else None)
        ax.plot([p[0], p[0]], [p[1], p[1]], [p[2], min(third_dim)], color=colors[n],
                linestyle='--', alpha=0.5, lw=0.75)
        # put a cross at the bottom of the plot
        ax.scatter(p[0], p[1], min(third_dim), color=colors[n], marker='x', alpha=0.5, lw=0.75)

# get the percentages of the explained variance
percentages = pcoa_results.proportion_explained

ax.set_xlabel(f'PC1: {percentages[0] * 100:.2f}%')
ax.set_ylabel(f'PC2: {percentages[1] * 100:.2f}%')
ax.set_zlabel(f'PC3: {percentages[2] * 100:.2f}%')

ax.view_init(45, 45, 0)
plt.show()
fig.savefig('fig 2.b.png', dpi=300, transparent=True)
# save a svg
fig.savefig('fig 2.b.svg')
# save a tif
fig.savefig('fig 2.b.tif', dpi=300, transparent=True)

for angle in tqdm(range(0, 360*3 + 1)):
    # Normalize the angle to the range [-180, 180] for display
    angle_norm = (angle + 180) % 360 - 180

    # Cycle through a full rotation of elevation, then azimuth, roll, and all
    elev = azim = roll = 0
    if angle <= 360:
        elev = angle_norm
    elif angle <= 360*2:
        azim = angle_norm
    else:
        elev = azim = roll = angle_norm

    # Update the axis view and title
    ax.view_init(elev, azim, roll)
    plt.draw()
    plt.pause(0.01)
    # fig.tight_layout()
    if angle < 10:
        fig.savefig(f'test\\000{angle}.png', dpi=300)
    elif angle < 100:
        fig.savefig(f'test\\00{angle}.png', dpi=300)
    elif angle < 1000:
        fig.savefig(f'test\\0{angle}.png', dpi=300)
    else:
        fig.savefig(f'test\\{angle}.png', dpi=300)

plt.close(fig)

fig, ax = plt.subplots()
percentages = [p*100 for p in percentages[:5]]
ax.plot(percentages, marker='x')
ax.set_xticks(range(5))
ax.set_xticklabels(['PC1', 'PC2', 'PC3', 'PC4', 'PC5'], rotation=45)
ax.set_ylabel('Explained Variance (%)')
ax.fill_between(range(3), 0, [p * 100 for p in pcoa_results.proportion_explained[:3]], color='gray', alpha=0.25)
ax.annotate(f'Explain {round(sum([p * 100 for p in pcoa_results.proportion_explained[:3]]), 2)}%',
            (2, sum([p * 100 for p in pcoa_results.proportion_explained[:3]]) / 6), ha='center')
plt.tight_layout()
fig.savefig('fig 2.b explained variance.png', dpi=300, transparent=True)
# save a svg
fig.savefig('fig 2.b explained variance.svg')
# save a tif
fig.savefig('fig 2.b explained variance.tif', dpi=300, transparent=True)


# make the gif
images = []
for filename in os.listdir('test'):
    images.append(Image.open(f'test\\{filename}'))

images[0].save('fig 2.b.gif', save_all=True, append_images=images[1:], duration=100, loop=0)
