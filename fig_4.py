import numpy as np
import matplotlib.pyplot as plt
from data.explorer import dt_loader, true_loader
from scipy.interpolate import splev, splrep
from scipy.stats import ttest_ind


# La figure 4.a est définie tel que:
# Abondances relatives des BSR (bleu) et BD (violet) selon le temps sous condition d'oscillation anoxique/oxique.
# La courbe verte montre l'évolution du rapport d'abondance DB/SRB tout au long de l'expérience.

metabolisms = {
    'BSR': ['Desulfatiglans', 'Desulfatirhabdium', 'Desulfobacter', '[Desulfobacterium]_catecholicum_group',
            'Desulfoconvexum', 'Desulfofrigus', 'Desulfopila', 'Desulfosarcina', 'Desulfosarcinaceae', 'Desulfuromonas',
            'Fusibacter', 'Desulfatiferula', 'Desulfobacula', 'Desulfobulbus', 'Desulfocapsa', 'Desulfococcus',
            'Desulfofaba', 'Desulfoluna', 'Desulfomonile', 'Desulfonema', 'Desulforhopalus', 'Desulfospira',
            'Desulfosporosinus', 'Desulfotalea', 'Desulfovibrio', 'Desulfuromusa', 'Dethiosulfatibacter'],
    'BD': ['Ardenticatenales', 'Flavirhabdus', 'Fluviicola', 'Gammaproteobacteria Incertae Sedis', 'Hyphomonadaceae',
           'Marinibacterium', 'Marinobacter', 'Pseudoalteromonas', 'Pseudomonas', 'Saccharospirillaceae',
           'Sulfurimonas', 'Arcobacter', 'Bacillus', 'Thiobacillus', 'Paracoccus', 'Vibrio']
}

bsr = []
bd = []

times = [5, 5.4, 7, 8, 10, 11, 15]

abundancies = {
    'BSR': [],
    'BD': []
}
sd = {
    'BSR': [],
    'BD': []
}
data, names = true_loader()
# reverse the keys and values of names
r_names = {v: k for k, v in names.items()}

for t in times:
    temp_bsr = []
    temp_bd = []
    temp_bd_sd = []
    temp_bsr_sd = []
    for name in metabolisms['BSR']:
        try:
            taxo = r_names[name]
            abundance_data = data[(data['Taxonomy'] == taxo) &
                                  (data['Time'] == t) &
                                  (data['Condition'] == 'OSCIL')]['Abundance']
            temp_bsr.append(abundance_data.mean())
            temp_bsr_sd.append(abundance_data.std())
            bsr.append(name)
        except:
            # print(name)
            pass
    for name in metabolisms['BD']:
        try:
            taxo = r_names[name]
            abundance_data = data[(data['Taxonomy'] == taxo) &
                                  (data['Time'] == t) &
                                  (data['Condition'] == 'OSCIL')]['Abundance']
            temp_bd.append(abundance_data.mean())
            temp_bd_sd.append(abundance_data.std())
            bd.append(name)
        except:
            pass
    abundancies['BSR'].append(np.sum(temp_bsr))
    abundancies['BD'].append(np.sum(temp_bd))
    sd['BSR'].append(np.mean(temp_bsr_sd))
    sd['BD'].append(np.mean(temp_bd_sd))
print(abundancies)
print(sd)
fig, ax = plt.subplots()

# Bar plot for BSR
ax.bar([str(t) for t in times], [ab * 100 for ab in abundancies['BSR']], color='r', label='BSR')

# Error bars for BSR
for n, (srb, s) in enumerate(zip(abundancies['BSR'], sd['BSR'])):
    ax.errorbar(n, srb * 100, yerr=s * 100, fmt='none', color='darkred', capsize=4)

# Bar plot for BD
ax.bar([str(t) for t in times], [ab * 100 for ab in abundancies['BD']], color='b',
       bottom=[ab * 100 for ab in abundancies['BSR']], label='BD')

# Error bars for BD
for n, (srb, bd, s) in enumerate(zip(abundancies['BSR'], abundancies['BD'], sd['BD'])):
    ax.errorbar(n, (bd + srb) * 100, yerr=s * 100, fmt='none', color='darkblue', capsize=4)

# Highlight periods
ax.axvspan(2, 3, color='gray', alpha=0.25)
ax.axvspan(4, 5, color='gray', alpha=0.25)

# Secondary axis for the ratio of BD/BSR
ax2 = ax.twinx()
x = range(len(times))
y = [abundancies['BD'][i] / abundancies['BSR'][i] for i in range(len(times))]
spline = splrep(x, y)
x_interp = np.linspace(min(x), max(x), 100)
y_interp = splev(x_interp, spline)
ax2.plot(x_interp, y_interp, color='g', label='Ratio BD/BSR')

ax2.errorbar(n, np.array(abundancies['BD']) / np.array(abundancies['BSR']),
             yerr=np.array(sd['BD']) / np.array(sd['BSR']), fmt='none', color='darkgreen', capsize=4)

ax2.set_ylabel('Ratio BD/BSR', color='darkgreen')
ax2.tick_params(axis='y', labelcolor='darkgreen')

# Remove top spines
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Set labels
ax.set_xlabel('Temps (jour)')
ax.set_ylabel('Abondance Relative (%)')
# ax.legend()

# Adjust layout
fig.tight_layout()

# Show plot
plt.show()

# Save plot
fig.savefig('figure 4.a.png', dpi=300, transparent=True)
fig.savefig('figure 4.a.svg', transparent=True)
fig.savefig('figure 4.a.tif', dpi=300, transparent=True)
del ax, ax2, fig

# student test
t, p = ttest_ind(abundancies['BSR'][:2], abundancies['BSR'][3:])
print(f"BSR: t: {t}, p: {p}")

ratio = np.array(abundancies['BD']) / np.array(abundancies['BSR'])
t, p = ttest_ind(ratio[:2], ratio[3:])
print(f"Ratio: t: {t}, p: {p}")

########################################################################################################################
# La figure 4.b est définie tel que:
# Progress of the DB/SRB abundance ratio along the experiment under anoxic/oxic oscillation condition (A/O, green) and under either anoxic (A, red) or oxic (O, blue) permanent conditions
# Initialize empty lists for BSR and BD
bsr = []
bd = []

# Initialize dictionaries for abundancies and standard deviations
abundancies = {
    'BSR ANOXIC': [],
    'BSR OSCIL': [],
    'BSR OXIC': [],
    'BD ANOXIC': [],
    'BD OSCIL': [],
    'BD OXIC': []
}
sd = {
    'BSR ANOXIC': [],
    'BSR OSCIL': [],
    'BSR OXIC': [],
    'BD ANOXIC': [],
    'BD OSCIL': [],
    'BD OXIC': []
}

# Iterate over the times
for t in times:
    # Iterate over the conditions
    for cond in ['ANOXIC', 'OSCIL', 'OXIC']:
        # Determine the condition
        c = 'OSCIL' if cond == 'ANOXIC' and t in [5, 5.4, 7] else cond
        # Initialize temporary lists for BSR and BD
        temp_bsr = []
        temp_bd = []
        temp_bd_sd = []
        temp_bsr_sd = []
        # Iterate over the names in BSR
        for name in metabolisms['BSR']:
            try:
                # Get the taxonomy
                taxo = r_names[name]
                # Get the abundance data
                abundance_data = data[(data['Taxonomy'] == taxo) &
                                      (data['Time'] == t) &
                                      (data['Condition'] == c)]['Abundance']
                # Append the mean and standard deviation to the temporary lists
                temp_bsr.append(abundance_data.mean())
                temp_bsr_sd.append(abundance_data.std())
                # Append the name to BSR
                bsr.append(name)
            except:
                pass
        # Iterate over the names in BD
        for name in metabolisms['BD']:
            try:
                # Get the taxonomy
                taxo = r_names[name]
                # Get the abundance data
                abundance_data = data[(data['Taxonomy'] == taxo) &
                                      (data['Time'] == t) &
                                      (data['Condition'] == 'OSCIL')]['Abundance']
                # Append the mean and standard deviation to the temporary lists
                temp_bd.append(abundance_data.mean())
                temp_bd_sd.append(abundance_data.std())
                # Append the name to BD
                bd.append(name)
            except:
                pass

        # Append the sum of the temporary lists to the abundancies dictionary
        abundancies[f'BSR {cond}'].append(np.sum(temp_bsr))
        abundancies[f'BD {cond}'].append(np.sum(temp_bd))
        # Append the mean of the temporary lists to the sd dictionary
        sd[f'BSR {cond}'].append(np.mean(temp_bsr_sd))
        sd[f'BD {cond}'].append(np.mean(temp_bd_sd))

# Create a new figure and axes
fig, ax = plt.subplots()
colors = ['b', 'g', 'r']
colors_std = ['darkblue', 'darkgreen', 'darkred']
conditions = ['ANOXIC', 'OSCIL', 'OXIC']
# Plot the DB/SRB ratio for each condition
for cond in conditions:
    x = range(len(times))
    y = [abundancies[f'BD {cond}'][i] / abundancies[f'BSR {cond}'][i] for i in range(len(times))]
    spline = splrep(x, y)
    x_interp = np.linspace(min(x), max(x), 100)
    y_interp = splev(x_interp, spline)
    ax.plot(x_interp, y_interp, label=f'{cond} BD/BSR', color=colors[conditions.index(cond)])
    for n in range(len(times)):
        ax.errorbar([str(t) for t in times][n],
                    abundancies[f'BD {cond}'][n] / abundancies[f'BSR {cond}'][n],
                    yerr=sd[f'BD {cond}'][n] + sd[f'BSR {cond}'][n],
                    fmt='none', capsize=4, color=colors_std[conditions.index(cond)])

# Highlight the periods between 7 and 8 hours, and between 10 and 11 hours
ax.axvspan(2, 3, color='gray', alpha=0.25)
ax.axvspan(4, 5, color='gray', alpha=0.25)

# Set the labels for the axes
ax.set_xlabel('Temps (jour)')
ax.set_ylabel('Ratio BD/BSR')
# Remove the top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Adjust the layout
fig.tight_layout()
# Save the figure in different formats
fig.savefig('figure 4.b.png', dpi=300, transparent=True)
fig.savefig('figure 4.b.svg', transparent=True)
fig.savefig('figure 4.b.tif', dpi=300, transparent=True)
# Show the plot
plt.show()

# ########################################################################################################################
# La figure 4.c est définie tel que:
# Évolution de l'abondance relative des bactéries sulfurisées oxydantes (DB-SO, rouge) et des bactéries hétérotrophes (DB-H, bleu).
# Le taux d'expression des SRB (violet) correspond au rapport du nombre de copies du transcrit/gène dsrAB par rapport au nombre de copies du transcrit/gène 16S rRNA
# Define the list of sulfur oxidant bacteria
sulfo_oxidant = ['Sulfurimonas', 'Thiobacillus', 'Thiomicrospira', 'Paracoccus', 'Arcobacter']

# Define the SRB expression values
srb_expression = [2.603592457, 3.740459458, 4.125750314, 0.870961393, 1.087523143, 0.068564644, 1.189759624]

# Initialize empty lists for sulfur oxidant and BD bacteria
so = []
bd = []

# Initialize dictionaries for abundancies and standard deviations
abundancies = {
    'BD': [],
    'SO': []
}
sd = {
    'BD': [],
    'SO': []
}

# Iterate over the times
for t in times:
    # Initialize temporary lists for sulfur oxidant and BD bacteria
    temp_bsr = []
    temp_bd = []
    temp_bd_sd = []
    temp_bsr_sd = []

    # Iterate over the names in sulfur oxidant bacteria
    for name in sulfo_oxidant:
        try:
            # Get the taxonomy
            taxo = r_names[name]

            # Get the abundance data
            abundance_data = data[(data['Taxonomy'] == taxo) &
                                  (data['Time'] == t) &
                                  (data['Condition'] == 'OSCIL')]['Abundance']

            # Append the mean and standard deviation to the temporary lists
            temp_bsr.append(abundance_data.mean())
            temp_bsr_sd.append(abundance_data.std())

            # Append the name to BD
            bd.append(name)
        except:
            pass

    # Iterate over the names in BD
    for name in metabolisms['BD']:
        try:
            # Get the taxonomy
            taxo = r_names[name]

            # Get the abundance data
            abundance_data = data[(data['Taxonomy'] == taxo) &
                                  (data['Time'] == t) &
                                  (data['Condition'] == 'OSCIL')]['Abundance']

            # Append the mean and standard deviation to the temporary lists
            temp_bd.append(abundance_data.mean())
            temp_bd_sd.append(abundance_data.std())

            # Append the name to BD
            bd.append(name)
        except:
            pass

    # Append the sum of the temporary lists to the abundancies dictionary
    abundancies['SO'].append(np.sum(temp_bsr))
    abundancies['BD'].append(np.sum(temp_bd))

    # Append the mean of the temporary lists to the sd dictionary
    sd['SO'].append(np.mean(temp_bsr_sd))
    sd['BD'].append(np.mean(temp_bd_sd))

# Create a new figure and axes
fig, ax = plt.subplots()

# Make a barplot with the BSR at the bottom and BD at the top and a curve of the SRB expression
ax.bar([str(t) for t in times], [ab * 100 for ab in abundancies['SO']], color='gold', label='BD-OS')

# Add standard deviation
for n, (srb, s) in enumerate(zip(abundancies['SO'], sd['SO'])):
    ax.errorbar(n, srb * 100, yerr=s * 100, fmt='none', color='orange', capsize=4)

# Add BD bar to the plot
ax.bar([str(t) for t in times], [ab * 100 for ab in abundancies['BD']], color='b',
       bottom=[ab * 100 for ab in abundancies['SO']], label='BD-H')

# Add standard deviation for BD
for n, (srb, bd, s) in enumerate(zip(abundancies['SO'], abundancies['BD'], sd['BD'])):
    ax.errorbar(n, (bd + srb) * 100, yerr=s * 100, fmt='none', color='darkblue', capsize=4)

# Plot the SRB expression
ax2 = ax.twinx()
x = range(len(times))
y = srb_expression
spline = splrep(x, y)
x_interp = np.linspace(min(x), max(x), 100)
y_interp = splev(x_interp, spline)
ax2.plot(x_interp, y_interp,
         color='violet',
         label='Taux d\'expression des BSR')

# Add error bars for SRB expression
for n in range(len(times)):
    ax2.errorbar(n, srb_expression[n],
                 fmt='none', color='darkviolet', capsize=4)

# Set the y-label for the secondary axis
ax2.set_ylabel('Taux d\'expression des BSR', color='darkviolet')

# Color the ticks in the same color
ax2.tick_params(axis='y', labelcolor='darkviolet')

# Remove top spines
ax.spines['top'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Highlight the period between 7 and 8 hours
ax.axvspan(2, 3, color='gray', alpha=0.25)

# Highlight the period between 10 and 11 hours
ax.axvspan(4, 5, color='gray', alpha=0.25)

# Set the labels for the axes
ax.set_xlabel('Temps (jour)')
ax.set_ylabel('Abondance Relative (%)')

# Adjust the layout
fig.tight_layout()

# Save the figure in different formats
fig.savefig('figure 4.c.png', dpi=300, transparent=True)
fig.savefig('figure 4.c.svg', transparent=True)
fig.savefig('figure 4.c.tif', dpi=300, transparent=True)

# Show the plot
plt.show()

# student test
t, p = ttest_ind(abundancies['SO'][:2], abundancies['SO'][3:])
print(f"BD-SO: t: {t}, p: {p}")
t, p = ttest_ind(abundancies['BD'][:2], abundancies['BD'][3:])
print(f"BD-H: t: {t}, p: {p}")
########################################################################################################################
# La figure 4.d est définie tel que:
# Évolution des OTU (unités taxonomiques opérationnelles) les plus abondantes de DB-SO (Sulfurimonas, rouge) et de DB-H
# (Marinobacter, bleu). Sur tous les graphiques, les zones ombrées indiquent les périodes d'oxygénation.
# Define the abundance and standard deviation for Marinobacter and Sulfurimonas
marinobacter = [
    data[(data['Taxonomy'] == 'Bacteria|Proteobacteria|Gammaproteobacteria|Alteromonadales|Marinobacteraceae|Marinobacter') &
         (data['Time'] == t) &
         (data['Condition'] == 'OSCIL')]['Abundance'].mean()
    for t in times]  # Calculate the mean abundance of Marinobacter for each time point
sulfurimonas = [
    data[(data['Taxonomy'] == 'Bacteria|Campilobacterota|Campylobacteria|Campylobacterales|Sulfurimonadaceae|Sulfurimonas') &
         (data['Time'] == t) &
         (data['Condition'] == 'OSCIL')]['Abundance'].mean()
    for t in times]  # Calculate the mean abundance of Sulfurimonas for each time point
marinobacter_sd = [
    data[(data['Taxonomy'] == 'Bacteria|Proteobacteria|Gammaproteobacteria|Alteromonadales|Marinobacteraceae|Marinobacter') &
         (data['Time'] == t) &
         (data['Condition'] == 'OSCIL')]['Abundance'].std()
    for t in times]  # Calculate the standard deviation of Marinobacter abundance for each time point
sulfurimonas_sd = [
    data[(data[
              'Taxonomy'] == 'Bacteria|Campilobacterota|Campylobacteria|Campylobacterales|Sulfurimonadaceae|Sulfurimonas') &
         (data['Time'] == t) &
         (data['Condition'] == 'OSCIL')]['Abundance'].std()
    for t in times]  # Calculate the standard deviation of Sulfurimonas abundance for each time point

# Create a new figure and axes
fig, ax = plt.subplots()

# Plot the abundance of Marinobacter
x = range(len(times))  # Define the x-axis as the range of time points
y = [m * 100 for m in marinobacter]  # Define the y-axis as the abundance of Marinobacter
spline = splrep(x, y)  # Create a B-spline representation of the data
x_interp = np.linspace(min(x), max(x), 100)  # Create an array of evenly spaced values over the range of x
y_interp = splev(x_interp, spline)  # Evaluate the B-spline at the specified points
ax.plot(x_interp, y_interp, color='navy', label='Marinobacter')  # Plot the interpolated data
ax.errorbar([str(t) for t in times], [m * 100 for m in marinobacter], yerr=[m for m in marinobacter_sd], fmt='none',
            color='midnightblue', capsize=4)  # Add error bars to the plot

# Plot the abundance of Sulfurimonas
x = range(len(times))  # Define the x-axis as the range of time points
y = [s * 100 for s in sulfurimonas]  # Define the y-axis as the abundance of Sulfurimonas
spline = splrep(x, y)  # Create a B-spline representation of the data
x_interp = np.linspace(min(x), max(x), 100)  # Create an array of evenly spaced values over the range of x
y_interp = splev(x_interp, spline)  # Evaluate the B-spline at the specified points
ax.plot(x_interp, y_interp, color='gold', label='Sulfurimonas')  # Plot the interpolated data
ax.errorbar([str(t) for t in times], [s * 100 for s, m in zip(sulfurimonas, marinobacter)],
            yerr=[s for s in sulfurimonas_sd], fmt='none', color='y', capsize=4)  # Add error bars to the plot

# Highlight the periods between 7 and 8 hours, and between 10 and 11 hours
ax.axvspan(2, 3, color='gray', alpha=0.25)  # Highlight the period between 7 and 8 hours
ax.axvspan(4, 5, color='gray', alpha=0.25)  # Highlight the period between 10 and 11 hours

# Set the labels for the axes
ax.set_xlabel('Temps (jour)')  # Set the label for the x-axis
ax.set_ylabel('Abondance Relative %')  # Set the label for the y-axis

# Remove the top and right spines
ax.spines['top'].set_visible(False)  # Remove the top spine
ax.spines['right'].set_visible(False)  # Remove the right spine

# Adjust the layout
fig.tight_layout()  # Adjust the layout of the figure

# Save the figure in different formats
fig.savefig('figure 4.d.png', dpi=300, transparent=True)  # Save the figure as a PNG file
fig.savefig('figure 4.d.svg', transparent=True)  # Save the figure as a SVG file
fig.savefig('figure 4.d.tif', dpi=300, transparent=True)  # Save the figure as a TIFF file

# Show the plot
plt.show()  # Display the figure

# Close the figure
plt.close(fig)  # Close the figure
