#Changes for each rat: start_sample, end_sample, filename, state_data_path, selected_times

import mne
from pylab import *
from mne import io
import pandas as pd
import itertools
import os
import os.path
import numpy as np
from numpy import *
#import parameters
import os
import matplotlib.pyplot as plt
#prm = parameters.Parameters()
#from OpenEphys import *
import mne
import xlrd

import sys
import itertools
mne.viz.set_browser_backend('matplotlib', verbose=None)


number_of_channels = 16
sample_rate = 250.4
sample_datatype = 'int16'
display_decimation = 1

"To set start and end times, put sample start and end below"
start_sample=23737922
end_sample=25540802

tmin = start_sample/sample_rate
tmax = end_sample/sample_rate

"To load the data, put file location and name below using double back to front slash"
filename="C:\\TAINI_1048_UBE3A_426_EM1_TransA-2024_01_09-0001.dat"
subject_name = '426'
#loading .dat files for all mice
def load_dat(file_path):

    '''Load a .dat file by interpreting it as int16 and then de-interlacing the 16 channels'''

    print("Loading_" + file_path)

    # Load the raw (1-D) data
    dat_raw = np.fromfile(file_path, dtype=sample_datatype)

    # Reshape the (2-D) per channel data
    step = number_of_channels * display_decimation
    dat_chans = [dat_raw[c::step] for c in range(number_of_channels)]

    # Build the time array
    data=np.array(dat_chans)
    del(dat_chans)
    #channel_names=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16']
    global channel_names
    global channel_types
    channel_names=['EMG', '2', '3', '4', 'MC-L',
                           '6', '7', 'VC-L', 'MC-R', 'VC-R',
                           '11', '12', '13', '14', '15', '16']
    channel_types=['emg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg','eeg']

    # This creates the info that goes with the channels, which is names, sampling rate, and channel types
    info = mne.create_info(channel_names, sample_rate, channel_types)

    # This makes the object that contains all the data and info about the channels. Computations like plotting, averaging, power spectrums can be performed on this object
    custom_raw = mne.io.RawArray(data, info)
   
    # For testing code, take only a small segment of the data
    #custom_raw = custom_raw.crop(tmin=0, tmax=60) # this is for testing
    return custom_raw
custom_raw = load_dat(filename)
custom_raw = custom_raw.crop(tmin, tmax)


# Load sleep state score
state_data_path = "C:\\426.csv"
def load_state_data(state_data_path):

    # Load sleep state score from .csv
    global state_data
    state_data = pd.read_csv(state_data_path, delimiter=",") # read .csv file with sleep score
    state_data['onset'] = state_data.index*5 # make column of time (index of the dataframe)
    state_data['duration'] = 5 # set duration of event - set as 5 since the sleep score is calculated in epochs of 5 seconds

    state_data = state_data[['onset', 'duration', 'Manuscore']]

sleep_states = load_state_data(state_data_path)

#create annotations
annotations = mne.Annotations(onset = np.asarray(state_data["onset"]), duration =  np.asarray(state_data["duration"]),description = np.asarray(state_data["Manuscore"]))

#Plotting recording with annotations to select epochs
annotated_custom_raw = custom_raw.copy().set_annotations(annotations)

def plot_annotated_plot():
    annotated_custom_raw.plot(None, 60, 0, 16, scalings = "auto", order=[0,4,7,8,9], show_options = "true")

#convert annotations to events
events = mne.events_from_annotations(annotated_custom_raw)
array_events = list(events[0])

#create epochs for events
event_id = {"Wake": 1, "NREM":2, "REM":3}
epoch = mne.Epochs(raw = annotated_custom_raw, events = array_events, event_id = event_id, tmin=-0.1, tmax = 4.9, picks='eeg')

# calculate overall power spectra and plot for average
def plot_wake_psd():
    #calculate power
    wake_psd = epoch['Wake'].compute_psd(method="welch", fmax = 125)
    #plotting average
    psds, freqs = wake_psd.get_data(), wake_psd.freqs
    mean_psd = np.mean(psds, axis=0)
    psds = mean_psd.T
    plt.figure(figsize=(10, 6))
    plt.plot(freqs, psds)
    plt.yscale("log")  
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (dB)')
    plt.grid(True)
    plt.ylim(bottom=1)
    plt.xlim(left=0, right=50)
    plt.show()

def plot_REM_psd():
    #calculate power
    REM_psd = epoch['REM'].compute_psd(method="welch", fmax = 125)
    #plotting average
    psds, freqs = REM_psd.get_data(), REM_psd.freqs
    mean_psd = np.mean(psds, axis=0)
    psds = mean_psd.T
    plt.figure(figsize=(10, 6))
    plt.plot(freqs, psds)
    plt.yscale("log")  
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (dB)')
    plt.grid(True)
    plt.ylim(bottom=1)
    plt.xlim(left=0, right=50)
    plt.show()

def plot_NREM_psd():
    #calculate power
    NREM_psd = epoch['NREM'].compute_psd(method="welch", fmax = 125)
    #plotting average
    psds, freqs = NREM_psd.get_data(), NREM_psd.freqs
    mean_psd = np.mean(psds, axis=0)
    psds = mean_psd.T
    plt.figure(figsize=(10, 6))
    plt.plot(freqs, psds)
    plt.yscale("log")  
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density (dB)')
    plt.grid(True)
    plt.ylim(bottom=1)
    plt.xlim(left=0, right=50)
    plt.show()

 
#!!!Once you've run the first part of this script and change the selected_times - remove this comment and the three quotations above and at the end - to be able to run this second half.!!!


#used for noting times for each state  
#selected_wake_time = [1030, 1035, 5550, 7025, 7030]
#selected_NREM_time = [4195, 4200, 4205, 4210, 4215]
#selected_REM_time = [2755, 2760, 3135, 3140, 3145]

#update with visually selected time points! best to keep in order
#order specifc - therefore must
selected_times = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 3120, 3125, 3130, 3135, 3140, 3145, 3150, 3155, 3160, 3165]

epochs_list = []
for time in selected_times:
    tmax = 5.0
   
    #converts onset times to index
    start_index = int(time  * sample_rate)
    end_index = int((time + tmax) * sample_rate)
   
    # Sub-select data around the index
    epoch_data = annotated_custom_raw.get_data()[:, start_index:end_index]
   
    # Append the epoch to the list
    epochs_list.append(epoch_data)

# Concatenate the epochs in the list to create a new Raw object
select_data = np.concatenate(epochs_list, axis=1)
# Create a new Raw object from the concatenated data
info = mne.create_info(channel_names, sample_rate, channel_types)
select_annotated_raw = mne.io.RawArray(select_data, annotated_custom_raw.info)
#plot new rat object to check
select_annotated_raw.plot(None, 60, 0, 16, scalings = "auto", order=[0,4,7,8,9], show_options = "true")

#re-create_annotations
# Create an empty list to store annotations
annotations_list = []

# Define parameters for the first set of annotations (5x5 seconds with label 'label_1')
#change label to match activity with the earliest time points
start_time = 0
duration = 5
label_1 = 'NREM'

# Create annotations for the first set
#change the expected range(n)
for i in range(10):
    onset = start_time + i * duration
    annotations_list.append((onset, duration, label_1))

# Define parameters for the second set of annotations (5x5 seconds with label 'label_2')
#change label to match activity with the middle time points
start_time = 50  # Assuming the first set ended at time 25
label_2 = 'REM'

# Create annotations for the second set
for i in range(10):
    onset = start_time + i * duration
    annotations_list.append((onset, duration, label_2))
   
# Define parameters for the second set of annotations (5x5 seconds with label 'label_2')
#change label to match activity with the latest time points
start_time = 100  # Assuming the first set ended at time 25
label_3 = 'Wake'

# Create annotations for the second set
for i in range(10):
    onset = start_time + i * duration
    annotations_list.append((onset, duration, label_3))


# Convert the list of annotations to mne.Annotations object
annotations = mne.Annotations(onset=[item[0] for item in annotations_list],
                              duration=[item[1] for item in annotations_list],
                              description=[item[2] for item in annotations_list])

global select_annotated_custom_raw

# Print the annotations
print(annotations)
select_annotated_custom_raw = select_annotated_raw.copy().set_annotations(annotations)

#check by plotting with annotations
select_annotated_custom_raw.plot(None, 60, 0, 16, scalings = "auto", order=[0,4,7,8,9], show_options = "true")


#create events
select_events = mne.events_from_annotations(select_annotated_custom_raw)
select_array_events = list(select_events[0])

#create epochs for events - SELECT CHANNEL
epoch = mne.Epochs(raw = select_annotated_custom_raw, events = select_array_events, event_id = event_id, tmin=-0.1, tmax = 4.9, picks='MC-R')


#calculating power spectra using MNE epoch.compute_psd() function for each state
#Wake
wake_psd = epoch['Wake'].compute_psd(method="multitaper", normalization='full', fmax = 125)
psds, freqs = wake_psd.get_data(), wake_psd.freqs
#plot of average power
mean_psd = np.mean(psds, axis=0)
psds = mean_psd.T
plt.figure(figsize=(10, 6))
plt.plot(freqs, psds)
plt.yscale("log")  
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density (dB)')
plt.grid(True)
plt.ylim(bottom=1)
plt.xlim(left=0, right=50)
plt.show()

# export data fro each brain state in a csv file for analysis
# Flatten freqs and psds arrays
wake_freqs_flat = freqs.flatten()
wake_psds_flat = psds.flatten()
# Create a DataFrame
df = pd.DataFrame({'freqs': wake_freqs_flat, 'psd': wake_psds_flat})
# Define the CSV file path
csv_filename = f"{subject_name}_PSD_Wake_data.csv"
# Save the DataFrame to a CSV file
df.to_csv(csv_filename, index=False)
print(f"Data saved to {csv_filename}")

#NREM
nrem_psd = epoch['NREM'].compute_psd(method="welch", fmax = 125)
psds, freqs = nrem_psd.get_data(), nrem_psd.freqs
#plot of average power
mean_psd = np.mean(psds, axis=0)
psds = mean_psd.T
plt.figure(figsize=(10, 6))
plt.plot(freqs, psds)
plt.yscale("log")  
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density (dB)')
plt.grid(True)
plt.ylim(bottom=1)
plt.xlim(left=0, right=50)
plt.show()

# export data fro each brain state in a csv file for analysis
# Flatten freqs and psds arrays
nrem_freqs_flat = freqs.flatten()
nrem_psds_flat = psds.flatten()
# Create a DataFrame
df = pd.DataFrame({'freqs': nrem_freqs_flat, 'psd': nrem_psds_flat})
# Define the CSV file path
csv_filename = f"{subject_name}_PSD_NREM_data.csv"
# Save the DataFrame to a CSV file
df.to_csv(csv_filename, index=False)
print(f"Data saved to {csv_filename}")

#REM
rem_psd = epoch['REM'].compute_psd(method="multitaper", normalization='full', fmax = 125)
psds, freqs = rem_psd.get_data(), rem_psd.freqs
#plot of average power
mean_psd = np.mean(psds, axis=0)
psds = mean_psd.T
plt.figure(figsize=(10, 6))
plt.plot(freqs, psds)
plt.yscale("log")  
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density (dB)')
plt.grid(True)
plt.xlim(left=0, right=50)
plt.show()

# export data fro each brain state in a csv file for analysis
# Flatten freqs and psds arrays
rem_freqs_flat = freqs.flatten()
rem_psds_flat = psds.flatten()
# Create a DataFrame
df = pd.DataFrame({'freqs': rem_freqs_flat, 'psd': rem_psds_flat})
# Define the CSV file path
csv_filename = f"{subject_name}_PSD_REM_data.csv"
# Save the DataFrame to a CSV file
df.to_csv(csv_filename, index=False)
print(f"Data saved to {csv_filename}")


#must repeat for the power spectra for each state (REM, NREM and Wake)
#power within specific frequency bands
delta_band = (1, 5)
theta_band = (5, 10)
sigma_band = (10, 16)
gamma_band = (30, 48)

#nrem
#obtain psd for frequency band
psd_data = {
    'delta': [],
    'theta': [],
    'sigma': [],
    'gamma': []
    }

#obtain psds (power data) and frequency from overall power spectra on selected REM epochs
psds, freqs = nrem_psd.get_data(), nrem_psd.freqs
psds = psds.T

# Calculate average power in each frequency band
delta_power = np.mean(psds[(freqs >= delta_band[0]) & (freqs < delta_band[1])])
theta_power = np.mean(psds[(freqs >= theta_band[0]) & (freqs < theta_band[1])])
sigma_power = np.mean(psds[(freqs >= sigma_band[0]) & (freqs < sigma_band[1])])
gamma_power = np.mean(psds[(freqs >= gamma_band[0]) & (freqs < gamma_band[1])])

# Store average power for each frequency band
psd_data['delta'].append(delta_power)
psd_data['theta'].append(theta_power)
psd_data['sigma'].append(sigma_power)
psd_data['gamma'].append(gamma_power)

# Convert lists to numpy arrays for easier manipulation
for band in psd_data:
    psd_data[band] = np.array(psd_data[band])

# Print average power in frequency bands
for band, power_data in psd_data.items():
    print(f'Average power in {band} band across epochs:', (math.log10(np.mean(power_data))))
         
#Plot average power in frequency bands
import seaborn as sns

#create dataframe with power
df = pd.DataFrame.from_dict(psd_data)

#seaborn figure customisation (Melissa Fasol Medium)
fig, axs = plt.subplots(1, 1, figsize=(20,15), sharex = True, sharey = True)
#barplot
sns.barplot(data=df)
sns.despine #remove figure border
axs.set_yscale('log')
#axs.set_ylim((10**0, 10**6)) #sets y axis scale
axs.set_xlabel("Frequency Band")
axs.set_ylabel("log$_{10}$ Power ($\mu V^2$)")


