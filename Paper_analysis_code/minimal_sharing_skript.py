"""
@author Sophia Gimple
"""
"""
This is the example skript for aligning Data recorded from the Medtronic Percept with EEG data recorded seperately based on TENS artefacts (we propose using stimulation settings 80Hz, 0.5ms, 1.5-2mA).
This skript includes functions for loading and cleaning the Percept data and automatically detecting TENS artefacts in a specific frequency band from the data.
Additionally the skript provides a function for automatically detect TENS artefacts in EEG data. It is assumed that the EEG data has been previously loaded and is saved in a numpy array with dimension NxM with N equal to the number of recorded channels.
Furthermore, a crosscorrelation function tests for crosscorrelation between the filtered EEG and Percept signal (extracting only a power trace that corresponds to the TENS stimulation frequency) and returns the computed lag that can be used for aligning the Percept and EEG data.


"""

import  os
import mne
import pyxdf
import sys
import numpy as np
import copy
from configparser import ConfigParser
from datetime import date
from os.path import join, exists
from os import listdir
import os
from mne.io import read_raw_fieldtrip
import mne
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import specgram
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, hilbert
from scipy import signal
import scipy.stats
from sklearn.preprocessing import normalize
from scipy import signal
import scipy.stats as stats
from sklearn.preprocessing import normalize
import math

def load_Percept(datapath,file_name,recording_number_percept=[0],recording_mode_percept="BrainSenseTimeDomain", plotting=1,fill_in_window=300):
    """
   #loads the Percept recordings given the relevant information.
   input:
   datapath: folder path where the percept file is saved
   file_name: name of the relevant file
   recording_number_percept: List of all recordings that should be exported. The recordings need to be of the same length e.g. recorded from different contacts at the same time.
   recording_mode_percept: Recording mode that was used for recording the relevant data. Options: BrainSenseTimeDomain or IndefiniteStreaming. BrainSenseTimeDomain is recommended.
   plotting: Boolean value, if set to 1 the raw data is plotted. Otherwise no data is plotted.
   fill_in_window: In case of missing values each value is imputed as the mean of the nan "fill_in_window" values before and after the missing value.
    """
    json_object=set_data_path(os.path.join(datapath, file_name))
    #load all recordings listed in recording_number_percept
    for ind,numb in enumerate(recording_number_percept):

        
        lfp_arr_elec1,Fs=new_lfp_checking_and_cleaning_all(json_object,recording_mode_percept,recording_number_percept[ind])
        nan_percentage=sum(np.isnan(lfp_arr_elec1))/len(np.isnan(lfp_arr_elec1))
        print("Electrode " + str(numb)+ "NAN percentage"+str(sum(np.isnan(lfp_arr_elec1))/len(np.isnan(lfp_arr_elec1))))
        print("Minutes recording:"+str((len(lfp_arr_elec1)/Fs)/60))
        if plotting:
            plt.figure()
            plt.plot(lfp_arr_elec1,color="orange")
           
        lfp_arr_elec1_corr=copy.deepcopy(lfp_arr_elec1)
    
        nan_val=0


    #if there are nan values in the percept recording they are filled in using the mean of the fill_in_window recording points before and after the missing value
        for i in range(0,len(lfp_arr_elec1)):
            if np.isnan(lfp_arr_elec1[i]):
                nan_val=nan_val+1

                comp_part=lfp_arr_elec1[max(0,i-fill_in_window):min(len(lfp_arr_elec1),i+fill_in_window)]
                mean_comp_part=np.nanmean(comp_part)
                lfp_arr_elec1_corr[i]=mean_comp_part
        if ind==0:
                lfp_elec_arr=np.zeros((len(recording_number_percept),len(lfp_arr_elec1_corr)))
        lfp_elec_arr[ind,0:len(lfp_arr_elec1_corr)]=lfp_arr_elec1_corr       

    return lfp_elec_arr,json_object,Fs,nan_percentage


def new_lfp_checking_and_cleaning_all(json_object,string,recording_num):
    """
 The function is adapted from the Percept Toolbox from Jeroen Habets and others
 Input:
 json_object: loaded datastructure as saved in .json format by the Medtronic toolbox
 string: recording modality of the relevant recording BrainSenseStreaming or IndefiniteStreaming
 recording_num: relevant recording number (integer)
    """
   
    #extract the streaming data, sampling rate and ticks from the structure
    streaming_data=json_object[string][recording_num]
    Fs = streaming_data['SampleRateInHz']
    ticksMsec = convert_list_string_floats(streaming_data['TicksInMses'])
    # compute the difference between time ticks
    ticksDiffs = np.diff(np.array(ticksMsec))

    #extract the package sizes
    packetSizes = convert_list_string_floats(streaming_data['GlobalPacketSizes'])

#In indefinite streaming commonly the data is seperated into two packages sent at the same time (time difference 0 between two packages). If this is the case we combine them into one package to be able to process this similarly to the BrainSenseStreaming data
    if string=="IndefiniteStreaming":
        

        packetSizes = convert_list_string_floats(streaming_data['GlobalPacketSizes'])

        for i_msec, msec in enumerate(ticksMsec):
            if i_msec==0:
                prev_msec=msec
                ticksDiffs_it=[]
                packageSize_it=[]
                ticksMsec_it=[]
                ticksMsec_it.append(ticksMsec[i_msec])
                #ticksDiffs_it.append(0)
                packageSize_it.append(packetSizes[i_msec])
            elif msec< prev_msec:
                print("error, there seems to be a problem with the packages. The new package is dated earlier than the previous package.")
            else:
                if msec!=prev_msec:
                    ticksDiffs_it.append(msec-prev_msec)
                    packageSize_it.append(packetSizes[i_msec])
                    ticksMsec_it.append(ticksMsec[i_msec])
                if msec==prev_msec:
                    packageSize_it[-1]=packageSize_it[-1]+packetSizes[i_msec]
                prev_msec=msec

        ticksDiffs=ticksDiffs_it
        packetSizes=packageSize_it
        ticksMsec=ticksMsec_it 


    values, counts = np.unique(ticksDiffs, return_counts=True)
    #We compare if the most common time difference between ticks matches with the expected tick difference
    det_ticksDiffs= values[counts.argmax()]
    print(f"The real time difference between packages is {det_ticksDiffs}")

    #We compute the mean package size of the loaded data
    calc_package_Sizes=Fs/((1/det_ticksDiffs)*1000)
    print(f"The real package size is {calc_package_Sizes}")
 

    #We check if there is data missing (irregularities in the calculated time differences)
    data_is_missing = (ticksDiffs != det_ticksDiffs).any()

    #We extract the lfp data
    lfp_data = streaming_data['TimeDomainData']


#we calculate the package sizes for uneven and even packages as they can differ
    
    if calc_package_Sizes%1==0:
        p_size_differs=0
        even_pack_size=calc_package_Sizes
        uneven_pack_size=calc_package_Sizes
    else:
        
        p_size_differs=1
        even_pack_sizes=packetSizes[::2]

        values, counts = np.unique(even_pack_sizes, return_counts=True)
        even_pack_size= values[counts.argmax()]

        uneven_pack_sizes=packetSizes[1::2]
        values, counts = np.unique(uneven_pack_sizes,return_counts=True)
        uneven_pack_size= values[counts.argmax()]
        print(f"package sizes differ between even {even_pack_size} and uneven {uneven_pack_size} packages")
    #We inform the user if data is missing or not        
    if data_is_missing:
        print('LFP Data is missing!! NaNs are filled in')
    else:
        print('No LFP data missing based on timestamp '
            'differences between data-packages')   

   
    #We check if there are too big packages and whether the package size variates between even and uneven packages and print the results
   
    too_large_package=np.nonzero(np.array(packetSizes)>max(even_pack_size,uneven_pack_size))[0]
    print("The following packages are too large:"+ str(too_large_package)) 


# We combine all data in one array (new_lfp_arr)   
    data_length_ms = ticksMsec[-1] + det_ticksDiffs - ticksMsec[0]  # length of a pakcet in milliseconds is always 250
    data_length_samples = int(data_length_ms / 1000 * Fs) + 1  # add one to calculate for 63 packet at end

    new_lfp_arr = np.array([np.nan] * data_length_samples)
    #if the first package is too large we want to handle this by cutting from the beginning, so we check for too big packages in the first and the consecutive packages
    number_wrong_first_packages=0
    if len(too_large_package)>0 and too_large_package[0]==0:
        consecutive=1
        it=0
        
        while consecutive and len(too_large_package)-1>it:
            if it+1==too_large_package[it+1]:
                it=it+1
            else:
                consecutive=0
    

        number_wrong_first_packages=it 
        first_right_package=it+1
        
        current_packetSizes=0
        for k in range(0,number_wrong_first_packages+1): current_packetSizes=current_packetSizes+int(packetSizes[k])

        print(f'The first {number_wrong_first_packages+1} packages are cutdown by {current_packetSizes-(it+1)*calc_package_Sizes} samples in total')
        start_pos_old_arr=int(current_packetSizes-(it+1)*calc_package_Sizes)
        temp_data=lfp_data[start_pos_old_arr:current_packetSizes]

        new_lfp_arr[:len(temp_data)]=temp_data

        i_old_arr=current_packetSizes
        i_new_arr = int((it+1)*calc_package_Sizes)
    #if a package that is not one of the first consecutive packages is too large we cut from the end. We have so for not observed this case
    else:
        current_packetSizes=int(packetSizes[0])
        temp_data=lfp_data[0:current_packetSizes]
        new_lfp_arr[:len(temp_data)]=temp_data

        i_old_arr=current_packetSizes
        i_new_arr = int(current_packetSizes)

    i_packet=number_wrong_first_packages+1

#For all other packages we also check the package size but cut from the end in case they are too big
    for i_diff, diff in enumerate(ticksDiffs[number_wrong_first_packages:]):

        if diff == det_ticksDiffs:
            # only lfp values, no nans if distance was 250 ms
            supposed_packetSize = int(packetSizes[i_packet])
            real_packetSize=supposed_packetSize
            

            # in case of very rare TOO LARGE packetsize (there is MORE DATA than expected based on the first and last timestamps), we cut it to the size we would expect in case of even or uneven packages
            if real_packetSize > math.ceil(calc_package_Sizes):
                print(f'UNKNOWN TOO LARGE DATAPACKET (# {i_diff}) IS CUTDOWN BY {real_packetSize - 63} samples')

                if i_packet%2==0:
                    supposed_packetSize = even_pack_size
                elif i_packet%2==1:
                    supposed_packetSize=uneven_pack_size
            
            elif real_packetSize < math.floor(calc_package_Sizes):
                print(f"potential data missing in packet # {i_diff}")
                if i_packet%2==0:
                    supposed_packetSize = even_pack_size
                elif i_packet%2==1:
                    supposed_packetSize=uneven_pack_size
            
            temp_data=lfp_data[i_old_arr:i_old_arr+real_packetSize]
            if i_new_arr+len(temp_data)>len(new_lfp_arr):
                 print("the lfp array is to short by"+str(i_new_arr+len(temp_data)-len(new_lfp_arr)))
                 new_lfp_arr[
                i_new_arr:int(i_new_arr + len(temp_data))
                ] = temp_data[0:-1*(i_new_arr+len(temp_data)-len(new_lfp_arr))]
            else:
                 new_lfp_arr[
                i_new_arr:int(i_new_arr + len(temp_data))
                ] = temp_data
            i_old_arr += int(real_packetSize)
            i_new_arr += int(supposed_packetSize)
            i_packet += 1
        #If data is missing we keep the initialized nan values
        elif diff>det_ticksDiffs:
            print(f'diff {i_diff} too big add NaNs by skipping {diff-det_ticksDiffs}')
            supposed_packetSize = int(packetSizes[i_packet])
            real_packetSize=supposed_packetSize

            
            msecs_missing = (diff - det_ticksDiffs)  # difference if one packet is missing is 500 ms
            
            secs_missing = msecs_missing / 1000
           
            samples_missing = int(secs_missing * Fs)
            print("num_samples_missing"+ str(samples_missing))
            # no filling with NaNs, bcs array is created full with NaNs
            i_new_arr += int(samples_missing)  # shift array index up by number of NaNs left in the array
            temp_data=lfp_data[i_old_arr:i_old_arr+real_packetSize]
            new_lfp_arr[
                int(i_new_arr):int(i_new_arr + len(temp_data))
            ] = temp_data
            i_old_arr += int(real_packetSize)
            i_new_arr += int(supposed_packetSize)
            i_packet += 1
    
    # cut off nan values at the end
    nonnan=0
    nan_count=0

    while nonnan==0:
            if np.isnan(new_lfp_arr[-1]):
                new_lfp_arr = new_lfp_arr[:-1]
                nan_count=nan_count+1
            
            else:
                nonnan=1
    return new_lfp_arr, Fs

def automatic_artefact_detection(data,sr,start_freq,end_freq,sig_val=1,rel_chan=0,len_art=10):
    """
    This function automatically detects artefacts in a specific frequency range in the signal.
    Input:
    data: data array from which the artefact will be detected. Each row is interpreted as a different channel.
    start_freq: Lower boundary of the frequency of the artefact e.g. 75, 79...
    end_freq: Upper boundary of the frequency of the relevant artefact e.g. 85, 81, ...
    sig_val: The artefacts are detected by comparing the mean of the values in a sliding window (m1) to the global mean of all data points (m_ges). sig_val is the minimum required distance in standard variations (std_ges) from the global mean
    for classification as an artefact. That means that for every classified point with index i the sliding window of interval [i:i+len_art] the mean value m1=mean(data[0,i:i+len(art)]) fulfills the condition  m1>mges+ sig_val* sd_ges.
    rel_chan: The channelnumber (saved as a row) in the data array that should be used for artefact detection. E.g. 0 corresponds to row 0 in the matrix data
    len_art: Size of the sliding window (in seconds) used for artefact detection
    """
    filtered_data=filter_frequency(data,sr,start_freq,end_freq)

    ges_mean=np.mean(filtered_data[rel_chan,:])
    ges_std=np.std(filtered_data[rel_chan,:])
    arr=np.zeros(len(filtered_data[rel_chan,:]))
    for i in range(0,len(filtered_data[rel_chan,:])-int(len_art*int(sr))):
        current_part=filtered_data[rel_chan,i:i+int(len_art*int(sr))]
        current_mean=np.mean(current_part)
        if current_mean>(ges_mean+sig_val*ges_std):
            arr[i]=1
    arr_start_art=[]
    for j in range(1,len(arr)):
        if arr[j]==1 and arr[j-1]==0:
            arr_start_art.append(j)
    plt.figure()
    plt.plot(stats.zscore(filtered_data[rel_chan,:]))

    plt.plot(arr,color="red")
    plt.show()
    return filtered_data,arr

#This function filters only the relevant frequency
def filter_frequency(array,sr,low_freq,high_freq):
    """
    This function filters the data keeping only the specified frequency range between low_freq and high_freq
    Input:
    array: input array or matrix
    sr: sampling rate of the input data
    low_freq: lower boundary of the relevant frequency band
    high_freq: upper boundary of the relevant frequency band
    """

    e=len(array)
    x=np.arange(0,sr*e)/sr
    b, a = butter(3, [low_freq/(sr/2),high_freq/(sr/2)], btype='bp')
    #filtered=filtfilt(b,a,array)
    filtered=filtfilt(b,a,array,padlen=0)
    array_filtered=np.abs(hilbert(filtered))
    return array_filtered

def set_data_path(datapath):
    """
    helperfunction for loading the .json file from datapath datapath. It returns the loaded data structure
    """

    if exists(join(datapath)):
            with open(datapath, 'r') as f:
                json_object = json.loads(f.read())
        
    else:
            raise ValueError(f'JSON file ({datapath}) not found ')
    return json_object


def convert_list_string_floats(string_list):
    """
    helperfunction that converts a string into a list of floats by splitting when delimiter ',' is in string
    """
    try:
        floats = [float(v) for v in string_list.split(',')]
    except:
        floats = [float(v) for v in string_list[:-1].split(',')]

    return floats

def crosscorrelate_and_compute_lags_and_scores(lfp_sync_part,eeg_sync_part,sr_lfp,sr_eeg,plotting):
    """
    #function that computes a lag between the percept and the eeg signal by crosscorrelating the two signals. It selects the best Percept channel and the best EEG channel specified in the channel_num_arr based on crosscorrelation scores.
    A positive lag corresponds to the Percept artefact being later than the corresponding EEG artefact.
    input:
    lfp_sync_part: Percept signal with the right frequency band extracted.
    eeg_sync_part: EEG signal with the right frequency band extracted
    sr_lfp: Sampling rate of the LFP signal.
    sr_eeg: Sampling rate of the EEG signal.
    plotting: 1- The function plots the original signals, the aligned signals, the crosscorrelation scores per lag and a histogram showing the distribution of scores.
    """

    if sr_eeg !=sr_lfp:
         print("The sampling rates between the eeg and lfp don't seem to match. The lfp data is automatically resampled.")
         lfp_sync_part=rowwise_resample(lfp_sync_part,sr_lfp,sr_eeg)
    eeg_sync_part1=copy.deepcopy(eeg_sync_part)

    score=np.zeros((lfp_sync_part.shape[0],eeg_sync_part.shape[0]))
    scores=np.zeros((lfp_sync_part.shape[0],eeg_sync_part.shape[0],len(lfp_sync_part[0,:])+len(eeg_sync_part[0,:])-1))
    lag_vals=np.zeros((lfp_sync_part.shape[0],eeg_sync_part.shape[0]))

    for it_lfp in range(0,lfp_sync_part.shape[0]):
            #normalize the Percept signal
            a=lfp_sync_part[it_lfp,:]

            norm_a = np.linalg.norm(a)
            a = a / norm_a

            for it_eeg in range(0,eeg_sync_part.shape[0]):

                #normalize the EEG signal
                b=eeg_sync_part1[it_eeg,:]
                norm_b = np.linalg.norm(b)

                b = b / norm_b

                #correlate and save the results in the matrices lag_vals, scores and score

                correlation = np.correlate(a, b, mode="full")

                lags = signal.correlation_lags(a.size, b.size, mode="full")
                lag = lags[np.argmax(correlation)]
                score[it_lfp,it_eeg]=np.max(correlation)
                scores[it_lfp,it_eeg,:]=correlation
                lag_vals[it_lfp,it_eeg]=lag

    #ind shows the best EEG and Percept channels
    ind=np.unravel_index(np.argmax(score), np.array(score).shape)
    print("The best LFP and EEG channels for crosscorrelation are"+str(ind[0])+" and "+str(ind[1]))
    lag=-1*(int(lag_vals[ind]))
    print("lag",lag," maximum correlation",np.max(score))

    #print("lfp electrode"+ str(ind))

        #align the Percept and EEG data based on the lag

    if lag>0:
                    eeg_sync_part_new=eeg_sync_part1[:,lag:]
                    lfp_sync_part_new=lfp_sync_part
                    
                    
    elif lag==0:
                    
                    eeg_sync_part_new=eeg_sync_part1
                    lfp_sync_part_new=lfp_sync_part

    else:

                    lfp_sync_part_new=lfp_sync_part[:,-lag:]

                    eeg_sync_part_new=eeg_sync_part1

        #plot the non-aligned EEG and Percept signal
    if plotting:

            plt.figure()

            plt.plot(stats.zscore(eeg_sync_part1[ind[1],:]),color="blue",alpha=0.7)
            plt.plot(stats.zscore(lfp_sync_part[ind[0]]),color="orange",alpha=0.7)
            #plot the non-aligned EEG and Percept signal
            plt.title("EEG and Percept before aligning")

            #plot the aligned EEG and Percept signal
            plt.figure()

            plt.plot(stats.zscore(eeg_sync_part_new[ind[1],:]),color="blue",alpha=0.7)


            plt.plot(stats.zscore(lfp_sync_part_new[ind[0],:]),color="orange",alpha=0.7)
            plt.title("EEG and Percept after aligning")
            
            #plot all scores of the best EEG and LFP signal channels
            plt.figure()
            plt.plot(scores[ind[0],ind[1],:])
            plt.title("Scores of the best LFP and EEG channel correlation")
            #plot the distributions of the scores of the best EEG and LFP channels
            plt.figure()
            x,y,z= plt.hist(scores[ind[0],ind[1],:], bins=30, color='skyblue', edgecolor='black')
            plt.vlines(np.mean(scores[ind[0],ind[1],:]),ymin=0,ymax=x.max())
            plt.vlines(np.mean(scores[ind[0],ind[1],:])+3*np.std(scores[ind[0],ind[1],:]),ymin=0,ymax=x.max())
            plt.title("Distribution of the scores of the best LFP and EEG channel correlation")
    return lag, ind, np.max(correlation)

def rowwise_resample(data,sr_old,sr_new):
    """
    This function resamples input data (array or 2D matrix) rowwise to a given sampling rate
    data: Input data (array or 2D matrix)
    sr_old: Sampling rate of the input data
    sr_new: Sampling rate to which the input data should be resampled
    """
    data=np.array(data)
    if np.array(data).ndim==1:
        upsample_data=signal.resample(data,int((len(data)/sr_old)*sr_new))
    elif np.array(data).ndim==2:
        for i in range(data.shape[0]):
            if i==0:
                upsample_data_temp=signal.resample(data[0,:],int((len(data[0,:])/sr_old)*sr_new))  
                upsample_data=np.zeros((data.shape[0],len(upsample_data_temp)))
                upsample_data[0,:]= upsample_data_temp
            else:
                upsample_data[i,:]=signal.resample(data[i,:],int((len(data[i,:])/sr_old)*sr_new))

    else:
         print("error the function rowwise_resample can only resample 1 or 2 dimensional data structures")
    return upsample_data