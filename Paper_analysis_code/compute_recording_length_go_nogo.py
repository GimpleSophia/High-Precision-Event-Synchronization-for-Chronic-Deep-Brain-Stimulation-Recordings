import os
import pandas as pd
import sys
import numpy as np
from config import DATA_PATH

sr_eeg=256


path = DATA_PATH
#List of participants included in this analysis
participant_list=["P01","P02","P03","P04","P05","P06","P07","P08","P09","P10","P11"]
#task category is GoNOGo here
task_category="GoNOGo"

#Initialize lists for collecting relevant information
#list of the lfp recording length in minutes
list_lfp_recording_length=[]
#list of the minimal recording length (between lfp or eeg recording) in minutes
list_min_recording_length=[]
#list of the participant abbreviations
list_participant=[]
#list of the recording names
list_file=[]
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
            #append participant names
            list_participant.append("P0"+str(p_i+1))
            subfolder_name=subfolder.split("/")[-1]
            #append the recording name
            list_file.append(subfolder_name)
            #append the number of lfp channels
            list_number_electrodes_lfp.append(filtered_lfp.shape[0])
            #append the lfp recording length in minutes
            list_lfp_recording_length.append((len(filtered_lfp[0,:])/sr_eeg)/60)
            #append the minimal lfp or eeg recording length in minutes
            list_min_recording_length.append((min(len(filtered_lfp[0,:]),len(filtered_eeg[0,:]))/sr_eeg)/60)

#combine all information in a dataframe
df_rl=pd.DataFrame({"participant":list_participant,
                    "recording":list_file,
                    "lfp recording length":list_lfp_recording_length,
                    "min recording length":list_min_recording_length,
                    "num lfp electrodes": list_number_electrodes_lfp
            
                    })  
#save the dataframe
df_rl.to_csv(path+"dataframe_recording_length"+task_category+"from channel.csv")   
