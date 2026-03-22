##This script computes offsets (lags) for predefined artifact groups 
# it extracts a dataframe with information about offsets, scores, impedance, recording mode
#usable for sync-test (3 artefact groups) or the GoNOGo task (2 artefact groups)


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
eeg_channels_nums=[0,21]

#change the task category you want to analyse (sync_test for synchronization testing or GoNOGo for the cognitvie task)
#task_category="sync_test"
task_category="GoNOGo"
#We include lfp impedance values that are specific to the participant, impedances for electrode contact pairs of directional leads are either
#aggregated by taking the minimal (1) the parallel combination (2) or maximum impedance (3)
impedance_ratio=1

#Initialize lists for collecting relevant information

#list of participant abbreviations per task recording
list_participant=[]

#list of task category name per recording
list_task_category=[]

# list of recording names
list_file=[]

# list of recording modalities (BrainSenseTimeDomain recording with DBS stimulation on (0mA)/off (0mA)/treat (at treatment stimulation settings) or indefinite streaming )
list_recording_modality_string=[]

#list of recording modalities indices (indef=0, treat=1, on=2, off=3)
list_recording_modality=[]

#list of computed offset values for artefact group 1
list_lag1=[]
#list of computed offset values for artefact group 2
list_lag2=[]
#list of computed offset values for artefact group 3
list_lag3=[]

#list of impedance value for the selected lfp channel pair for artefact group 1
list_imp1=[]
#list of impedance value for the selected lfp channel pair for artefact group 2
list_imp2=[]
#list of impedance value for the selected lfp channel pair for artefact group 3
list_imp3=[]

#list of values with maximal sample difference between offset values
list_lag_difference=[]
#list of values with maximal difference in ms between offset values
list_lag_difference_ms=[]
#list of selected lfp channel for synchronization of artefact group 1
list_lfp_channel_1=[]
#list of selected lfp channel for synchronization of artefact group 2
list_lfp_channel_2=[]
#list of selected lfp channel for synchronization of artefact group 3
list_lfp_channel_3=[]
#list of selected eeg channel for synchronization of artefact group 1
list_eeg_channel_1=[]
#list of selected eeg channel for synchronization of artefact group 2
list_eeg_channel_2=[]
#list of selected eeg channel for synchronization of artefact group 3
list_eeg_channel_3=[]
#list of maximal correlation value of crosscorrelation during synchronization for artefact group 1
list_correlation_1=[]
#list of maximal correlation value of crosscorrelation during synchronization for artefact group 2
list_correlation_2=[]
#list of maximal correlation value of crosscorrelation during synchronization for artefact group 3
list_correlation_3=[]
#list of the number of recording contacts (channels) in the lfp recording
list_number_electrodes_lfp=[]

#We iterate over participants
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
            #load the lfp channel names
            lfp_channels=np.load(subfolder+"/lfp_channels.npy")
            #load the impedances of all recording channel pairs
            if impedance_ratio==1:
                lfp_impedances=np.load(subfolder+"/lfp_impedances_monopolar_ratio.npy")
            else:
                lfp_impedances=np.load(subfolder+"/lfp_impedances_monopolar_diff.npy")
            # define the relevant eeg channels as defined above
            chan_start_eeg=eeg_channels_nums[0]
            chan_end_eeg=eeg_channels_nums[1]
            #For artefact group 1 compute the crosscorrelation, lag1_t= offset of group 1, ind_1= lfp and eeg channel indices selected for crosscorrelation, corr_1= corresponding correlaiton value of group 1
            lag1_t,ind_1,corr_1=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[0]:artefact_limits[1]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[0]:artefact_limits[1]],sr_eeg,sr_eeg,plotting=0)
            # select the impedance of the optimal lfp channel pair
            imp_1=lfp_impedances[ind_1[0]]
            #append the number of lfp channels
            list_number_electrodes_lfp.append(filtered_lfp.shape[0])
            # if we look at sync_test recordings we have 3 artefact groups
            if task_category=="sync_test":
                #For artefact group 2 compute the crosscorrelation, lag2_t= offset of group 2, ind_2= lfp and eeg channel indices selected for crosscorrelation, corr_2= corresponding correlaiton value of group 2
                lag2_t,ind_2,corr_2=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[1]:artefact_limits[2]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[1]:artefact_limits[2]],sr_eeg,sr_eeg,plotting=0)
                 # select the impedance of the optimal lfp channel pair
                imp_2=lfp_impedances[ind_2[0]]
                #For artefact group 3 compute the crosscorrelation, lag3_t= offset of group 3, ind_3= lfp and eeg channel indices selected for crosscorrelation, corr_3= corresponding correlaiton value of group 3
                lag3_t,ind_3,corr_3=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[2]:artefact_limits[3]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[2]:artefact_limits[3]],sr_eeg,sr_eeg,plotting=0)
                 # select the impedance of the optimal lfp channel pair
                imp_3=lfp_impedances[ind_3[0]]
                print("imp check",imp_1, imp_2, imp_3, max(lfp_impedances))
                #append the offset, impedance values to the corresponding initialized lists
                list_lag3.append(lag3_t)
                list_imp3.append(imp_3)
                list_correlation_3.append(corr_3)
                #save the offsets and indices in a list
                lag_list=[lag1_t,lag2_t,lag3_t]
                ind_list=[ind_1,ind_2,ind_3]
                list_lfp_channel_3.append(ind_3[0])
                list_eeg_channel_3.append(ind_3[1])
                
                #extract the recording modality from the subfolder name
                if "indef" in subfolder:
                    list_recording_modality.append(0)
                    list_recording_modality_string.append("indef")
                elif "treat" in subfolder:
                    list_recording_modality.append(1)
                    list_recording_modality_string.append("treat")
                elif "off" in subfolder:
                    list_recording_modality.append(2)
                    list_recording_modality_string.append("off")
                else:
                    list_recording_modality.append(3)
                    list_recording_modality_string.append("on")
                
                    

                
            # if we look at GoNOGo task recordings we have only 2 artefact groups
            else:
                #For artefact group 2 compute the crosscorrelation, lag2_t= offset of group 2, ind_2= lfp and eeg channel indices selected for crosscorrelation, corr_2= corresponding correlaiton value of group 2
                lag2_t,ind_2,corr_2=crosscorrelate_and_compute_lags_and_scores(filtered_lfp[:,artefact_limits[2]:artefact_limits[3]],filtered_eeg[chan_start_eeg:chan_end_eeg,artefact_limits[2]:artefact_limits[3]],sr_eeg,sr_eeg,plotting=0)
                # select the impedance of the optimal lfp channel pair
                imp_2=lfp_impedances[ind_2[0]]
                #save the offsets and indices in a list
                lag_list=[lag1_t,lag2_t]
                ind_list=[ind_1,ind_2]
                #extract the recording modality from the subfolder name and the number of lfp channels (3 per electrode for indefinite streaming)
                if len(lfp_channels)>2:
                    list_recording_modality_string.append("indef")
                elif "Go1" in subfolder:
                    list_recording_modality_string.append("treat")

                elif "Go4" in subfolder:
                    list_recording_modality_string.append("treat")
                elif "On_0" in subfolder:
                    list_recording_modality_string.append("On")
                else:
                    print ("labeling missing",subfolder)
                    list_recording_modality_string.append("Off")
                
            #append the offset and impedance values and correlation to the corresponding initialized lists    
            list_lag1.append(lag1_t)
            list_lag2.append(lag2_t)
            list_imp1.append(imp_1)
            list_imp2.append(imp_2)
            list_correlation_1.append(corr_1)
            list_correlation_2.append(corr_2)
            
            list_lfp_channel_1.append(ind_1[0])
            list_lfp_channel_2.append(ind_2[0])
            list_eeg_channel_1.append(ind_1[1])
            list_eeg_channel_2.append(ind_2[1])
            
            #extract and append the difference in offsets both in samples and milliseconds
            list_lag_difference.append(max(lag_list)-min(lag_list))
            list_lag_difference_ms.append(((max(lag_list)-min(lag_list))/256)*1000)
            #also save the corresponding participant, category and recording names
            list_participant.append("P0"+str(p_i+1))
            list_task_category.append(task_category)
            subfolder_name=subfolder.split("/")[-1]
            list_file.append(subfolder_name)
#combine all the information to dataframes 
if task_category=="sync_test":
    df_lags=pd.DataFrame({"participant":list_participant,
                    "recording":list_file,
                    "task":list_task_category,
                    "num lfp electrodes": list_number_electrodes_lfp,
                    "recording_modality":list_recording_modality_string,
                    "lag 1":list_lag1,
                    "lag 2":list_lag2,
                    "lag 3":list_lag3,
                    "lag difference":list_lag_difference,
                    "lag difference ms":list_lag_difference_ms,
                    "channel lfp 1":list_lfp_channel_1,
                    "channel lfp 2":list_lfp_channel_2,
                    "channel lfp 3":list_lfp_channel_3,
                    "impedance lfp 1":list_imp1,
                    "impedance lfp 2":list_imp2,
                    "impedance lfp 3":list_imp3,
                    "channel eeg 1":list_eeg_channel_1,
                    "channel eeg 2":list_eeg_channel_2,
                    "channel eeg 3":list_lfp_channel_3,
                    "correlation 1": list_correlation_1,
                    "correlation 2": list_correlation_2,
                    "correlation 3": list_correlation_3,
                    })  
else:
    df_lags=pd.DataFrame({"participant":list_participant,
                    "recording":list_file,
                    "task":list_task_category,
                    "num lfp electrodes": list_number_electrodes_lfp,
                    "recording_modality":list_recording_modality_string,
                    "lag 1":list_lag1,
                    "lag 2":list_lag2,
                    "impedance lfp 1":list_imp1,
                    "impedance lfp 2":list_imp2,
                    "lag difference":list_lag_difference,
                    "lag difference ms":list_lag_difference_ms,
                    "channel lfp 1":list_lfp_channel_1,
                    "channel lfp 2":list_lfp_channel_2,

                    "channel eeg 1":list_eeg_channel_1,
                    "channel eeg 2":list_eeg_channel_2,
                    "correlation 1": list_correlation_1,
                    "correlation 2": list_correlation_2,


                    })  

#save the dataframe with all the information 
if impedance_ratio==1:
    df_lags.to_csv(path+"dataframe_lag_shift_impedance_num_ratio_"+str(task_category)+"from channel"+str(eeg_channels_nums[0])+" to "+str(eeg_channels_nums[1])+".csv")   
else:
    df_lags.to_csv(path+"dataframe_lag_shift_impedance_num_diff_"+str(task_category)+"from channel"+str(eeg_channels_nums[0])+" to "+str(eeg_channels_nums[1])+".csv")   
