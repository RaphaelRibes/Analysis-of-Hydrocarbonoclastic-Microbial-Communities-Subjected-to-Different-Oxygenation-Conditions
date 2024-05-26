import pandas as pd


def dt_loader(path='data\\abundancebacter.csv'):
    """Charge les données avec un format plus facile à manipuler"""
    dt = pd.read_csv(path, sep=';')

    # Renomme la première colonne en 'ID'
    dt = dt.rename(columns={dt.columns[0]: 'id'})

    # Renomme 'barcode' en 'Barcode'
    dt = dt.rename(columns={'barcode': 'Barcode'})

    # Retire tous les 'Com' de la colonne 'Community'
    dt['Community'] = dt['Community'].str.replace('Com', '')

    # Néttoie la colonne 'Time'
    times = dt['Time'].unique().tolist()
    for t in times:
        if '+10h' in t:
            dt.loc[dt['Time'] == t, 'Time'] = '5.4'
        else:
            dt.loc[dt['Time'] == t, 'Time'] = t.replace('T', '')

    # Convertit la colonne 'Time' en float
    dt['Time'] = dt['Time'].astype(float)

    return dt


def true_loader():
    """
    This function loads the data using the dt_loader function, modifies the taxonomy, and creates a dictionary of names.

    Returns:
    dt (DataFrame): The modified data.
    names (dict): A dictionary where the keys are the original taxonomy and the values are the modified names.
    """

    # Load the data using the dt_loader function
    dt = dt_loader()

    # Combine the columns 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus' into a single column 'Taxonomy'
    # The values in 'Taxonomy' are the values in the other columns joined by a '|'
    # If 'uncultured' or 'Unknown_Family' is in the data, it is removed
    dt['Taxonomy'] = dt[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']].apply(lambda x: '|'.join(x),
                                                                                          axis=1)

    # Drop the columns 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'OTU', 'Barcode', 'Primer'
    dt = dt.drop(columns=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'OTU', 'Barcode', 'Primer'])

    # Initialize an empty list for the filtered data
    f_dt = []

    # For each unique value in 'Taxonomy', append the rows with that value to f_dt
    for d in dt['Taxonomy'].unique():
        f_dt.append(dt[dt['Taxonomy'] == d])

    # Initialize an empty dictionary for the names
    names = {}

    # For each unique value in 'Taxonomy', modify the name and add it to the names dictionary
    # The key is the original taxonomy and the value is the modified name
    for name in dt['Taxonomy'].unique():
        f_name = name.split('|')
        while f_name[-1] == 'Unknown_Family' or f_name[-1] == 'uncultured':
            f_name.pop(-1)
        if len(f_name) == 6:
            f_name = f_name[-1]
        else:
            f_name = 'Unknown ' + f_name[-1]
        names[name] = f_name

    # Return the modified data and the names dictionary
    return dt, names


if __name__ == '__main__':
    data = dt_loader('abundancebacter.csv')
    for col in data.columns:
        print(f"{col}: {len(data[col].unique())}")
