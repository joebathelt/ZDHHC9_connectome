#===============================================
# Analysis of connectome data using graph theory 

# Utility functions
def get_parcellation_labels(label_dict):
    import networkx as nx
    import pandas as pd
    numbers = list()
    labels = list()

    Labels = nx.read_gpickle(label_dict)
    for label in Labels:
        numbers.append(label)
        labels.append(Labels[label]['labels'])
    
    return pd.Series(labels,index=numbers,name='label')

def remove_non_cortical_ROIs(labels,adjmat):
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
    
    return (labels, adjmat)

def generate_ROI_file(FreeSurfer_ROI_file):
	"""
	This script generates a dictionary of ROIs found in the FreeSurfer parcellation filename
	"""
	from nipype.interfaces.freesurfer import MRIConvert
	mc = MRIConvert()
	mc.inputs.in_file = FreeSurfer_ROI_file
	mc.inputs.out_type = 'niigz'
	mc.run()

	import nipype.interfaces.cmtk as cmtk
	rg = cmtk.ROIGen()
	rg.inputs.aparc_aseg_file = FreeSurfer_ROI_file.split('.')[0] + '_out.nii.gz'
	rg.inputs.use_freesurfer_LUT = True
	out_file = rg.run()

	return out_file

#===============================================
# Calculating graph measures in each group
## Loading the connectome data
import scipy.io as sio
import pandas as pd
import numpy as np
import bct

groups = ['patient','control']

for group in groups:
    if group == 'patient':
        participants = ['z1','z2','z3','z4','z5','z6','z8']
    elif group == 'control':
        participants = ['c1','c2','c3','c5','c6','c7','c8']
        
    all_measures = np.empty(shape=[68,len(participants),5])
    adjmats =  np.empty(shape=[68,68,len(participants)])
    counter = 0

    for participant in participants:
        adjmat = sio.loadmat(participant + '_FA.mat')
        adjmat = adjmat['adjacency_matrix']
        labels = get_parcellation_labels(generate_ROI_file(FreeSurfer_ROI_file)).values
        labels,adjmat = remove_non_cortical_ROIs(labels,adjmat)
        all_measures[:,counter,0] = bct.degrees_und(adjmat)
        all_measures[:,counter,1] = bct.strengths_und(adjmat)
        all_measures[:,counter,2] = bct.clustering_coef_wu(adjmat)
        all_measures[:,counter,3] = bct.betweenness_wei(adjmat)
        all_measures[:,counter,4] = bct.efficiency_wei(adjmat,local=True)
        adjmats[:,:,counter] = adjmat
        counter += 1
        
        
    mean_measures = np.mean(all_measures,axis=1)
    if group == 'patient':
        patient = pd.DataFrame(mean_measures, index=labels,columns=['patient.NodeDegree','patient.Strength','patient.ClustCoeff','patient.BetweenCent','patient.LocEff'])
        patient_measures = all_measures
        patient_adjmats = adjmats
    elif group == 'control':
        control = pd.DataFrame(mean_measures, index=labels,columns=['control.NodeDegree','control.Strength','control.ClustCoeff','control.BetweenCent','control.LocEff'])
        control_measures = all_measures
        control_adjmats = adjmats
       
# Getting the expression values
import pandas as pd
lh_dataframe = pd.read_csv('ZDHHC9_expression_lh.csv',sep=',')
rh_dataframe = pd.read_csv('ZDHHC9_expression_rh.csv',sep=',')

expression = pd.concat([lh_dataframe,rh_dataframe])
expression = expression.sort(columns='ZDHHC9')

# Combining expression values and graph metrics by the label of the cortical region
combined = expression.join(patient)
combined = combined.join(control)
combined = combined.drop(['annotationNumber','r','g','b','a','hexColor','Label.name'],axis=1)
combined.to_csv('Graph_measures_and_expression.csv')

#===============================================
# Comparison of edge weight by region
control_left = list()
control_right = list()
control_interhem = list()
control_total = list()

patient_left = list()
patient_right = list()
patient_interhem = list()
patient_total = list()

for i in range(0,len(participants)):
    control_right.append(np.mean(control_adjmats[0:34,0:34,i]))
    patient_right.append(np.mean(patient_adjmats[0:34,0:34,i]))
    
    control_left.append(np.mean(control_adjmats[34:,34:,i]))
    patient_left.append(np.mean(patient_adjmats[34:,34:,i]))
    
    control_interhem.append(np.mean(np.vstack([control_adjmats[0:34,34:,i],control_adjmats[34:,0:34,i]])))
    patient_interhem.append(np.mean(np.vstack([patient_adjmats[0:34,34:,i],patient_adjmats[34:,0:34,i]])))
    
    control_total.append(np.mean(control_adjmats[:,:,i]))
    patient_total.append(np.mean(patient_adjmats[:,:,i]))

left = np.hstack([control_left,patient_left])
right = np.hstack([control_right,patient_right])
interhem = np.hstack([control_interhem,patient_interhem])
total = np.hstack([control_total,patient_total])
participants = ['c1','c2','c3','c5','c6','c7','c8','z1','z2','z3','z4','z5','z6','z8']

import pandas as pd
df = pd.DataFrame()
df['Participant.ID'] = participants 
df['Participant.Group'] = np.hstack([np.repeat('control',7),np.repeat('ZDHHC9',7)])
df['left_weight'] = left
df['right_weight'] = right
df['interhemispheric_weight'] = interhem
df['total'] = total
df.to_csv('weight_by_region.csv')




