def get_parcellation_labels_from_LUT(parcellation_lookup_table):

    #================================
    # Parcellation labels
    # This function get the labels corresponding to the values in the look-up table. 
    # This function expects a FreeSurfer-style lookup table, i.e. label numbers in the first column,
    # corresponding labels in the second column
    #================================

    """
    inputs:    
    parcellation_lookup_table: filename of the look-up table

    outputs:
    pandas dataframe containing the label number and the corresponding label
    """
    # Getting the labels and order
    import pandas as pd
    import re

    numbers = list()
    labels = list()

    file = open(parcellation_lookup_table,'r')

    for line in file:
        if len(line) > 4 and not re.search('#',line):
            content = list()

            entries = line.split(' ')
            for entry in entries:
                if entry:
                    content.append(entry)

            if len(content[0]) > 3:
                parts = content[0].split('\t')
                content[0] = parts[0]
                content[1] = parts[1]

            numbers.append(content[0])
            labels.append(content[1])

    return pd.Series(labels,index=numbers,name='label')

def get_parcellation_labels_from_gpickle(label_dict):

    #================================
    # Parcellation labels
    # This function reads the labels of ROIs from a gpickle dictionary
    #================================

    """ 
    inputs:
    label_dict: gpickle dictionary containing the ROI information

    outputs:
    pandas series containing the labels and index numbers of all ROIs
    """

    import pandas as pd
    import networkx as nx
    
    numbers = list()
    labels = list()

    Labels = nx.read_gpickle(label_dict)
    for label in Labels:
        numbers.append(label)
        labels.append(Labels[label]['labels'])
    
    return pd.Series(labels,index=numbers,name='label')

def remove_non_cortical_ROIs(labels,adjmat):

    #================================
    # Remove non-cortical label from a list of parcellation labels
    # This function removes all non-cortical ROIs from a list of ROIs. 
    # The cortical ROIs have to be labelled with 'ctx-' to be recognised.
    # The corresponding entries are also removed from the adjacency matrix.
    #================================

    """
    inputs:
    labels: list of labels
    adjmat: adjacency_matrix in the same order as the labels

    outputs:
    list of labels without non-cortical ROIs
    adjacency matrix without entries for non-cortical ROIs
    """

    import numpy as np
    import re 
    
    counter = 0
    non_cortical_labels = list()

    for label in labels:
        if not re.search('ctx',label):
            non_cortical_labels.append(counter)

        counter += 1

    labels = np.delete(labels,non_cortical_labels)
    adjmat = np.delete(adjmat,non_cortical_labels,axis=0)
    adjmat = np.delete(adjmat,non_cortical_labels,axis=1)
    
    return (labels,adjmat)

def get_label_order_from_file(label_order_file):
    #================================
    # Get labels from a saved file
    # This function reads a label file and return the labels in the same order as a list
    # This may be useful for re-order labels according to an external file
    #================================

    """ 
    inputs:
    label_order_file: text file containing the labels

    outputs:
    list of labels in order of appearance in the text file
    """

    import pandas as pd
    new_order = pd.read_csv(label_order_file,names=['Label'])
    new_order = new_order['Label'].values.tolist()
    return new_order