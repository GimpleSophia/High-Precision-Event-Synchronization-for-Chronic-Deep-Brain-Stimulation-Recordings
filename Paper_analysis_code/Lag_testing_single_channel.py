##This script computes offsets (lags) for predefined artifact groups 
# iterating over all EEG channels individually
#it extracts the mean/max/min offsets for using each EEG channel alone
import os
import pandas as pd
import sys


from config import DATA_PATH
from minimal_sharing_skript import *

path = DATA_PATH
#List of participants included in this analysis
participant_list=["P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11"]
# EEG sampling rates
sr_eeg=256
#EEG channels (excluding non-electrophysiological e.g. accelorometer channels)
eeg_channel_num=[0,21]
#create a list of all eeg channels
eeg_channel_list=range(eeg_channel_num[0],eeg_channel_num[1])
#change the task category you want to analyse (sync_test for synchronization testing or GoNOGo for the cognitvie task)
task_category="sync_test"
#task_category="GoNOGo"

#Initialize lists for collecting relevant information
#list of current eeg channel
list_channel=[]
#mean value of all offset differences in samples when synchronizing only with this eeg channel
all_lag_list_mean=[]
#max value of all offset differences in samples when synchronizing only with this eeg channel
all_lag_list_max=[]
#min value of all offset differences in samples when synchronizing only with this eeg channel
all_lag_list_min=[]
#mean value of all offset differences in milliseconds when synchronizing only with this eeg channel
all_lag_list_mean_ms=[]
#max value of all offset differences in milliseconds when synchronizing only with this eeg channel
all_lag_list_max_ms=[]
#min value of all offset differences in milliseconds when synchronizing only with this eeg channel
all_lag_list_min_ms=[]


#iterate over all eeg_channels
for eeg_chan in eeg_channel_list:
    #create a temporary list with the offset differences in samples and milliseconds of all recordings
    print(eeg_chan)
    lag_list_diff=[]
    lag_list_diff_ms=[]
    # iterate over all participants
    #Data is saved in the following structure:
    # DATA_PATH
    #          "P0x" folder for each participant
    #               GoNOGo (task folder where applicable)
    #                   -P0XGoNOGoY... task run for participant X and task number X
    #               sync_test (synchronization test folder where applicable)
    #                   -P0Xsync_test_ZZ... run for participant X and stimulation condition ZZ (indef/off/on/treat)

    for p_i,parti in enumerate(participant_list):
        print(parti)
        #find the subfolder corresponding to the single recordings and iterate over them
        parti_path=path+parti+"/"+task_category+"/"
        if os.path.isdir(parti_path):
            


            subfolders = [ f.path for f in os.scandir(parti_path) if f.is_dir() ]
            #for each recording
            for f_i,subfolder in enumerate(subfolders):
                #load the eeg data
                filtered_eeg=np.load(subfolder+"/eeg_filtered_sliced.npy")
                #load the lfp data
                filtered_lfp=np.load(subfolder+"/lfp_filtered_sliced.npy")
                #load the group limits of the first, second and third artefact group
                artefact_limits=np.load(subfolder+"/group_limits.npy")
                
                #lfp_channels=np.load(subfolder+"/lfp_channels.npy")
                #we determine the correct start and end index for selecting the right eeg channel
                chan_start_eeg=eeg_chan
                chan_end_eeg=eeg_chan+1
                #For artefact group 1 compute the crosscorrelation, lag1_t= offset of group 1, ind_1= lfp and eeg channel indices selected for crosscorrelation, corr_1= corresponding correlaiton value of group 1
                lag1_t,ind_1,corr_1=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[0]:artefact_limits[1]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[0]:artefact_limits[1]],sr_eeg,sr_eeg,plotting=0)
                 # if we look at sync_test recordings we have 3 artefact groups
                if task_category=="sync_test":
                    #For artefact group 2 compute the crosscorrelation, lag2_t= offset of group 2, ind_2= lfp and eeg channel indices selected for crosscorrelation, corr_2= corresponding correlaiton value of group 2
                    lag2_t,ind_2,x=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[1]:artefact_limits[2]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[1]:artefact_limits[2]],sr_eeg,sr_eeg,plotting=0)
                    #For artefact group 3 compute the crosscorrelation, lag3_t= offset of group 3, ind_3= lfp and eeg channel indices selected for crosscorrelation, corr_3= corresponding correlaiton value of group 3
                    lag3_t,ind_3,x=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[2]:artefact_limits[3]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[2]:artefact_limits[3]],sr_eeg,sr_eeg,plotting=0)
                    #save the offsets in a list
                    lag_list=[lag1_t,lag2_t,lag3_t]
                # if we look at GoNOGo task recordings we have only 2 artefact groups
                else:
                    #For artefact group 2 compute the crosscorrelation, lag2_t= offset of group 2, ind_2= lfp and eeg channel indices selected for crosscorrelation, corr_2= corresponding correlaiton value of group 2
                    lag2_t,ind_2,x=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[2]:artefact_limits[3]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[2]:artefact_limits[3]],sr_eeg,sr_eeg,plotting=1)
                    #save the offsets in a list
                    lag_list=[lag1_t,lag2_t]
                #extract and append the difference in offsets both in samples and milliseconds 
                lag_list_diff.append(np.max(lag_list)-min(lag_list))
                lag_list_diff_ms.append(((np.max(lag_list)-min(lag_list))/256)*1000)
    #append all the information per channel to the previously initialized lists
    list_channel.append(eeg_chan)
    all_lag_list_mean_ms.append(np.mean(lag_list_diff_ms))
    all_lag_list_min_ms.append(np.min(lag_list_diff_ms))
    all_lag_list_max_ms.append(np.max(lag_list_diff_ms))
    all_lag_list_mean.append(np.mean(lag_list_diff))
    all_lag_list_min.append(np.min(lag_list_diff))
    all_lag_list_max.append(np.max(lag_list_diff))
#combine all the information to dataframes 
df_lags_single_channel=pd.DataFrame({"channel":list_channel,
                    "mean_lag":all_lag_list_mean,
                    "min_lag":all_lag_list_min,
                    "max_lag":all_lag_list_max,
                    "mean_lag_ms":all_lag_list_mean_ms,
                    "min_lag_ms":all_lag_list_min_ms,
                    "max_lag_ms":all_lag_list_max_ms,
                    
                    })  
 

#save the dataframe with all the information 
df_lags_single_channel.to_csv(path+"dataframe_lag_shift_single_channel"+task_category+"from channel"+str(eeg_channel_num[0])+" to "+str(eeg_channel_num[1])+".csv")   


