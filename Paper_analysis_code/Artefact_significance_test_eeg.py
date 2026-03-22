import os
import numpy as np
from scipy import stats
import pandas as pd
from config import DATA_PATH


path = DATA_PATH

#List of participants included in this analysis
participant_list=["P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11"]
# task category here synchronization test =sync_test
task_category="sync_test"

#Initialize lists for collecting relevant information
#list of participant abbreviations per task recording
list_participant=[]
# list of recording names
list_recording=[]
#list of the eeg channel names
list_channel=[]
#list of t-values
list_t_value=[]
#list of p-values
list_p_value=[]
#list of mean values of the artefact timewindow
list_mean_artefact=[]
#list of mean values of the no artefact timewindows
list_mean_no_artefact=[]
#list of standard deviation values of the artefact timewindow
list_std_artefact=[]
#list of standard deviation values of the no artefact timewindows
list_std_no_artefact=[]
#EEG channels (excluding non-electrophysiological e.g. accelorometer channels)
eeg_channels_nums=[0,21]

#We iterate over participants
#Data is saved in the following structure:
# DATA_PATH
#          "P0x" folder for each participant
#               GoNOGo (task folder where applicable)
#                   -P0XGoNOGoY... task run for participant X and task number X
#               sync_test (synchronization test folder where applicable)
#                   -P0Xsync_test_ZZ... run for participant X and stimulation condition ZZ (indef/off/on/treat)

for p_i,parti in enumerate(participant_list):
    #find the subfolder corresponding to the single recordings and iterate over them
    parti_path=path+parti+"/"+task_category+"/"
    if os.path.isdir(parti_path):
        
        subfolders = [ f.path for f in os.scandir(parti_path) if f.is_dir() ]
        #for each recording
        for f_i,subfolder in enumerate(subfolders):
            #load the eeg data
            #load the eeg channels
            eeg_channels=np.load(subfolder+"/eeg_channels.npy")
            #load the start timewindow filtered and hilbert transformed eeg signal without arefacts
            eeg_start_segment=np.load(subfolder+"/no_artefact_start_segment_eeg.npy")
            #load the start end end filtered and hilbert transformed timewindow eeg signal without arefacts
            eeg_end_segment=np.load(subfolder+"/no_artefact_end_segment_eeg.npy")
            #load the timewindow filtered and hilbert transformed signal with arefacts
            eeg_artefact_segment=np.load(subfolder+"/artefact_segment_eeg.npy")
            #iterate over all relevant eeg channels
            for i_c, chan in enumerate(eeg_channels[eeg_channels_nums[0]:eeg_channels_nums[1]]):
                    #extract the start and end segment without artefacts for this channel
                    eeg_start_segment_here=eeg_start_segment[i_c]
                    eeg_end_segment_here=eeg_end_segment[i_c]
                    #append the two segments
                    eeg_no_artefact_here=np.concatenate((eeg_start_segment_here,eeg_end_segment_here))
                    #extract the segment with artefacts for this channel
                    eeg_artefact_here=eeg_artefact_segment[i_c,:]
                    #compute and append the mean of the artefact signal
                    list_mean_artefact.append(np.mean(eeg_artefact_here))
                    #compute and append the mean of the non-artefact signal
                    list_mean_no_artefact.append(np.mean(eeg_no_artefact_here))
                    #compute and append the standard deviation of the artefact signal
                    list_std_artefact.append(np.std(eeg_artefact_here))
                    #compute and append the standard deviation of the non-artefact signal
                    list_std_no_artefact.append(np.std(eeg_no_artefact_here))
                    #compute a ttest comparing artefact and non-artefact power signals and append the t and p-values
                    t_value,p_value=stats.ttest_ind(eeg_no_artefact_here,eeg_artefact_here)
                    list_t_value.append(t_value)
                    list_p_value.append(p_value)
                    #append the participant name
                    list_participant.append(parti)
                    #append the channel name
                    list_channel.append(chan)
                    #append the recording name
                    subfolder_name=subfolder.split("/")[-1]
                    list_recording.append(subfolder_name)
#combine all information in a dataframe               
df_eeg_artefact_test=pd.DataFrame({"participant":list_participant,
                                    "recording":list_recording,
                                    "channel":list_channel,
                                    "mean artefact":list_mean_artefact,
                                    "mean no artefact": list_mean_no_artefact,
                                    "std artefact":list_std_artefact,
                                    "std no artefact": list_std_no_artefact,
                                    "t-value": list_t_value,
                                    "p-value":list_p_value  
                                    })  
#save the dataframe
df_eeg_artefact_test.to_csv(path+"dataframe_eeg_artefact_testing.csv")   


