# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 01:26:31 2024

@author: Odysseas Makariou
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import classification_report
from sklearn.metrics import cohen_kappa_score
import seaborn as sns

# General Plot settings
font = 'calibri'
fontSize = 15

# concensus plot settings
line_color = 'k'
line_w = 1
line_s = '-'
fill_color = 'r'
fill_opacity = 0.3

subjects = ["426A","448A","444A","428A","455A","434A","435A","447A","454A"]  # Set subject names
subject_num = len(subjects)

# Organise data for all subjects
data = [] 
data_wake =[]
data_nrem = []
data_rem =[]

# List to store (subject, Cohen's Kappa) tuples
kappa_values_general = []
kappa_values_wake = []
kappa_values_nrem = []
kappa_values_rem = []

for subject in subjects:
    # Obtain data from manual sleepscoring spreadsheet
    results = pd.read_excel("Manual_SleepScoring_Christia_6pm_to_8pm_ODY_2.xlsx", sheet_name=subject) # Extract respective spreadsheet
    autoscore = results["Autoscore"].values # Extract autoscore data for particular subject
    manuscore = results["Manuscore"].values # Extract manuscore data for particular subject
    epoch = results["Epoch"].values # Extract epoch data for particular subject

    # Identify inconsistencies in results data (NaN)
    nan_idx_auto = np.argwhere(np.isnan(autoscore)).flatten()
    nan_idx_manu = np.argwhere(np.isnan(manuscore)).flatten()
    nan_idx_combine = np.unique(np.concatenate((nan_idx_auto, nan_idx_manu)))

    # Remove inconsistencies from data
    autoscore = np.delete(autoscore, nan_idx_combine)
    manuscore = np.delete(manuscore, nan_idx_combine)
    epoch = np.delete(epoch, nan_idx_combine)
    
    # Generate general Cohen's Kappa for individual subject
    kappa_general = cohen_kappa_score(manuscore,autoscore)
    print(f"General Cohen's kappa for {subject}: {kappa_general:.5f}")
    
    # Generate Cohen's Kappa for each brain state for individual subject
    # Wake
    manuscore_wake = manuscore.copy()
    autoscore_wake = autoscore.copy()
    # Reassign variables for Wake (3) and not Wake (4) entries
    manuscore_wake[manuscore_wake != 2] = 4
    manuscore_wake[manuscore_wake == 2] = 3
    autoscore_wake[autoscore_wake != 2] = 4
    autoscore_wake[autoscore_wake == 2] = 3
    # Calculate Cohen's Kappa value for Wake state for individual subject
    kappa_individual_wake = cohen_kappa_score(manuscore_wake,autoscore_wake)
    print(f"Cohen's kappa for {subject} - Wake: {kappa_individual_wake:.5f}")
    
    # NREM
    manuscore_nrem = manuscore.copy()
    autoscore_nrem = autoscore.copy()
    # Reassign variables for NREM (3) and not NREM (4) entries
    manuscore_nrem[manuscore_nrem != 1] = 4
    manuscore_nrem[manuscore_nrem == 1] = 3
    autoscore_nrem[autoscore_nrem != 1] = 4
    autoscore_nrem[autoscore_nrem == 1] = 3
    # Calculate Cohen's Kappa value for Wake state for individual subject
    kappa_individual_nrem = cohen_kappa_score(manuscore_nrem,autoscore_nrem)
    print(f"Cohen's kappa for {subject} - NREM: {kappa_individual_nrem:.5f}")
    
    # REM
    manuscore_rem = manuscore.copy()
    autoscore_rem = autoscore.copy()
    # Reassign variables for REM (3) and not REM (4) entries
    manuscore_rem[manuscore_rem != 0] = 4
    manuscore_rem[manuscore_rem == 0] = 3
    autoscore_rem[autoscore_rem != 0] = 4
    autoscore_rem[autoscore_rem == 0] = 3
    # Calculate Cohen's Kappa value for Wake state for individual subject
    kappa_individual_rem = cohen_kappa_score(manuscore_rem,autoscore_rem)
    print(f"Cohen's kappa for {subject} - REM: {kappa_individual_rem:.5f}")
    
    # Append individual subject's Cohen's Kappa value to list
    kappa_values_general.append((subject,kappa_general))
    kappa_values_wake.append((subject,kappa_individual_wake))
    kappa_values_nrem.append((subject,kappa_individual_nrem))
    kappa_values_rem.append((subject,kappa_individual_rem))

    # Append sleepscoring data
    data.append({"autoscore": autoscore, "manuscore": manuscore, "epoch": epoch})

# Calculate averaged Cohen's Kappa value for each brain state and general with standard error
kappa_values_array = np.array([value[1] for value in kappa_values_general])
kappa_averaged_general = np.mean(kappa_values_array)
sem_general = np.std(kappa_values_array) / np.sqrt(len(kappa_values_array))

kappa_values_array = np.array([value[1] for value in kappa_values_wake])
kappa_averaged_wake = np.mean(kappa_values_array)
sem_wake = np.std(kappa_values_array) / np.sqrt(len(kappa_values_array))

kappa_values_array = np.array([value[1] for value in kappa_values_nrem])
kappa_averaged_nrem = np.mean(kappa_values_array)
sem_nrem = np.std(kappa_values_array) / np.sqrt(len(kappa_values_array))

kappa_values_array = np.array([value[1] for value in kappa_values_rem])
kappa_averaged_rem = np.mean(kappa_values_array)
sem_rem = np.std(kappa_values_array) / np.sqrt(len(kappa_values_array))

print(f"Average Cohen's kappa - General: {kappa_averaged_general:.5f}, SEM: {sem_general:.5f}")
print(f"Average Cohen's kappa - Wake: {kappa_averaged_wake:.5f}, SEM: {sem_wake:.5f}")
print(f"Average Cohen's kappa - NREM: {kappa_averaged_nrem:.5f}, SEM: {sem_nrem:.5f}")
print(f"Average Cohen's kappa - REM: {kappa_averaged_rem:.5f}, SEM: {sem_rem:.5f}")
    
# Combine data from all subjects
autoscore_all = data[0]["autoscore"]
manuscore_all = data[0]["manuscore"]

for i in range(0, subject_num-1):
    autoscore_all = np.concatenate((autoscore_all, data[i+1]["autoscore"]))
    manuscore_all = np.concatenate((manuscore_all, data[i+1]["manuscore"]))
    
# Generate confusion matrix for manuscore and autoscore data
confusion_matrix = metrics.confusion_matrix(manuscore_all,autoscore_all)
sns.heatmap(confusion_matrix, cmap='Blues', annot=True, fmt='g')
plt.xticks(np.arange(3)+0.5,["WAKE", "NREM", "REM"])
plt.yticks(np.arange(3)+0.5,["WAKE", "NREM", "REM"])
plt.xlabel("Autoscore")
plt.ylabel("Manuscore")
plt.title("Confusion Matrix")   

# Generate summary report
report = classification_report(manuscore_all,autoscore_all, output_dict=True)
final_report = pd.DataFrame(report).transpose()
print(final_report)

# Generate Cohen's Kappa for the collective of all subjects
kappa_all_general = cohen_kappa_score(manuscore_all,autoscore_all)
print(f"Overall Cohen's kappa: {kappa_all_general:.5f}")

# Generate Cohen's Kappa for each brain state for the collective of all subjects
# Wake
manuscore_all_wake = manuscore_all.copy()
autoscore_all_wake = autoscore_all.copy()
# Reassign variables for REM (3) and not REM (4) entries
manuscore_all_wake[manuscore_all_wake != 2] = 4
manuscore_all_wake[manuscore_all_wake == 2] = 3
autoscore_all_wake[autoscore_all_wake != 2] = 4
autoscore_all_wake[autoscore_all_wake == 2] = 3

# NREM
manuscore_all_nrem = manuscore_all.copy()
autoscore_all_nrem = autoscore_all.copy()
# Reassign variables for REM (3) and not REM (4) entries
manuscore_all_nrem[manuscore_all_nrem != 1] = 4
manuscore_all_nrem[manuscore_all_nrem == 1] = 3
autoscore_all_nrem[autoscore_all_nrem != 1] = 4
autoscore_all_nrem[autoscore_all_nrem == 1] = 3

# REM
manuscore_all_rem = manuscore_all.copy()
autoscore_all_rem = autoscore_all.copy()
# Reassign variables for REM (3) and not REM (4) entries
manuscore_all_rem[manuscore_all_rem != 0] = 4
manuscore_all_rem[manuscore_all_rem == 0] = 3
autoscore_all_rem[autoscore_all_rem != 0] = 4
autoscore_all_rem[autoscore_all_rem == 0] = 3

kappa_all_wake = cohen_kappa_score(manuscore_all_wake,autoscore_all_wake)
kappa_all_nrem = cohen_kappa_score(manuscore_all_nrem,autoscore_all_nrem)
kappa_all_rem = cohen_kappa_score(manuscore_all_rem,autoscore_all_rem)
print(f"Overall Cohen's kappa - Wake: {kappa_all_wake:.5f}")
print(f"Overall Cohen's kappa - NREM: {kappa_all_nrem:.5f}")
print(f"Overall Cohen's kappa - REM: {kappa_all_rem:.5f}")




plt.show()