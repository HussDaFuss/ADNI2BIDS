import pdb
from os.path import join, exists, isdir, dirname
from os import makedirs
import pandas as pd
import os
import glob
import numpy as np
import csv
import json
import subprocess
from datetime import datetime

import bids


def pick_desired_modalites(path_data, dataset_name, df_collec, mybids, dict_file_names, mri_df=None, pet_df=None, adnimerge_df=None):
    '''
    This function selects all desired modalities of the ADNI-downloades Dataset (Resting State fMRI, MP-Rage and PET). It also returns a dataframe with all remaining images and their corresponding Subjects, Visits and Modalities.
    
    Parameters
    ----------
    path_data : TYPE str 
        DESCRIPTION. Root directory where all other files and directoies are nested.
    dataset_name : TYPE str
        DESCRIPTION. Name of the Dataset batch that are currently beeing added to the Data Base
    df_collec : TYPE Pandas DataFrame
        DESCRIPTION. Dataframe containing all image with its corresponding data. It is downloaded from Ida Loni - Adni webpage along with the actual image data and the idaSearch_3_21_2022.csv file.
        Column Fields: Image Data ID,"Subject","Group","Sex","Age","Visit","Modality","Description","Type","Acq Date","Format","Downloaded", Alternative_ID, Dicom_path, nifti_path, session
    mybids : TYPE BIDSLayout object
        DESCRIPTION. BIDSLayout object containing the already existing BIDS DataBase in the path_data/BIDS directory.
    dict_file_names : TYPE dict
        DESCRIPTION. Dictionary Containing the desired modalities and it's suffix according to the BIDS structure naming criteria.
    mri_df : TYPE Pandas DataFrame, optional
        DESCRIPTION. DataFrame containing all MRI images present in the ADNI repository along with their corresponding metadata. This Dataframe is downloaded from ADNI as a MRIList.csv. The default is None.
    pet_df : TYPE Pandas Dataframe, optional
        DESCRIPTION. DataFrame containing all PET images present in the ADNI repository along with their corresponding metadata. This Dataframe is downloaded from ADNI as a PETList.csv. The default is None.
    adnimerge_df : TYPE Pandas DataFrame, optional
        DESCRIPTION. DataFrame containing Key ADNI tables merged into one table [ADNI1,GO,2,3]. The default is None.

    Returns
    -------
    rf : TYPE DataFrame
        DESCRIPTION. DataFrame where all Remaining Files in the ADNI Dataset file structure are listed and detailed.
        Columns: 'Image' , 'Modality', 'Subject', 'Session, dicom_path, subject_id, session_id, nifti_path, json_path
    excluded_modalities : TYPE DataFrame
        DESCRIPTION. 

    '''

    

    # Aquí es pot llegir tot el que tinc i fer una petita base de dades.
    # valid_modalities = ['PET', 'fMRI', 'MP', 'florbetapir', 'MT1', 'fcMRI', 'AV45', 'T2', 'Hipp']

    rf = pd.DataFrame(columns=df_collec.columns.tolist() + ['dicom_path', 'nifti_path', 'subject_id', 'session_id'])
    # Here all modalities that aren't fMRI, MP Rage or pet are discarded. The folders in the ADNI Database are deleted.
    excluded_modalities = []
    for subject_folder in os.listdir(join(path_data, dataset_name)):
        if not os.path.isdir(join(path_data, dataset_name, subject_folder)): continue
        for modality_folder in os.listdir(join(path_data, dataset_name, subject_folder)):
            if any([vm.lower() in modality_folder.lower() for vm in dict_file_names.keys()]) and os.path.isdir(join(path_data, dataset_name, subject_folder, modality_folder)):
                for session_folder in os.listdir(join(path_data, dataset_name, subject_folder, modality_folder)):
                    for imageid in os.listdir(join(path_data,dataset_name, subject_folder, modality_folder, session_folder)):
                        final_img_path = join(path_data, dataset_name, subject_folder, modality_folder, session_folder, imageid)

                        if not os.path.isdir(final_img_path): continue

                        # Convert only images that are in the samples file and ignore those that aren't even if
                        # they are downloaded. Ths way we could have all data on the same "Download" folder and
                        # just use different sample files to process it.
                        image_df = df_collec.loc[df_collec['Image Data ID'] == imageid]
                        if len(image_df) > 0:
                            image_df = image_df.iloc[0]
                        else:
                            continue
                        subject_id = 'ADNI' + replace_underscore(image_df['Subject'])


                        # Calcular l'identificador de sessio
                        if pet_df is not None and (len(pet_df[pet_df['Image ID'] == int(imageid[1:])]['Visit_Month_Number']) == 1):
                            adni_subject = adnimerge_df.loc[adnimerge_df['PTID'] == subject_folder]
                            date_format = "%Y-%m-%d"
                            scan_date = datetime.strptime(pet_df.loc[pet_df['Image ID'] == int(imageid[1:])]['Scan Date'].iloc[0], date_format)
                            diff_days = [(scan_date - datetime.strptime(d, date_format)).days for d in adni_subject['EXAMDATE']]
                            idx_adni_subject = np.argmin(np.abs(np.array(diff_days)))
                            session_id = 'm' + str(adni_subject.iloc[idx_adni_subject]['Month']).zfill(3)

                        elif mri_df is not None and (len(mri_df[mri_df['IMAGEUID'] == int(imageid[1:])]['Visit_Month_Number']) == 1):
                            adni_subject = adnimerge_df.loc[adnimerge_df['PTID'] == subject_folder]
                            date_format = "%Y-%m-%d"
                            scan_date = datetime.strptime(mri_df.loc[mri_df['IMAGEUID'] == int(imageid[1:])]['SCANDATE'].iloc[0], date_format)
                            diff_days = [(scan_date - datetime.strptime(d, date_format)).days for d in adni_subject['EXAMDATE']]
                            idx_adni_subject = np.argmin(np.abs(np.array(diff_days)))
                            session_id = 'm' + str(adni_subject.iloc[idx_adni_subject]['Month']).zfill(3)
                            # session_id = 'm' + str(mri_df.loc[mri_df['IMAGEUID'] == int(imageid[1:])].iloc[0]['Visit_Month_Number']).zfill(3)

                        elif adnimerge_df is not None and (len(adnimerge_df[adnimerge_df['IMAGEUID'] == int(imageid[1:])]['Month']) == 1):
                            session_id = 'm' + str(adnimerge_df.loc[adnimerge_df['IMAGEUID'] == int(imageid[1:])].iloc[0]['Month']).zfill(3)

                        else:
                            sid = image_df['Visit']
                            if sid in ['bl', 'sc']:
                                sid = '000'
                            elif 'y' in sid:
                                year = int(sid.split('y')[1])
                                sid = str(year*12).zfill(3)
                            elif isinstance(sid, int):
                                sid = str(int(sid)).zfill(3)

                            session_id = 'm' + str(sid)

                        file_dict = [d for k, d in dict_file_names.items() if k in modality_folder][0]

                        image_df['dicom_path'] = final_img_path
                        image_df['subject_id'] = subject_id
                        image_df['session_id'] = session_id
                        image_df['nifti_path'] = mybids.build_path({**{'subject': subject_id, 'session': session_id, 'extension': 'nii.gz'}, **file_dict})
                        image_df['json_path'] = mybids.build_path({**{'subject': subject_id, 'session': session_id, 'extension': 'json'}, **file_dict})

                        rf = pd.concat([rf,image_df.to_frame().T])
            else:
                # if 'PET_Brain' in modality_folder:
                #     pdb.set_trace()
                if modality_folder not in excluded_modalities:
                    excluded_modalities.append(modality_folder)
            #     shutil.rmtree(join(path_data, dataset_name, subject_folder, modality_folder),
            #                   ignore_errors=False)

    # rf = rf.reset_index()
    # rf = rf.rename(columns={"index": "image_ID"})
    return rf, excluded_modalities


def remove_repeated_files(sample_filter_df, mybids):
    '''
    

    Parameters
    ----------
    sample_filter_df : TYPE
        DESCRIPTION.
    mybids : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''

    if 'ImageID' not in mybids.get_entities():
        return sample_filter_df
        # To remove
        # print('\n Hello! It seems that it is your first time creating the dataset. Press "y" to confirm or "n" to exit; otherwise'
        #       ' just make sure that all json files related to images contain the field "ImageID".')
        # continue_flag = input('Are you sure to continue? (Continue: y. Quit: n)')
        # if continue_flag == 'n':
        #     exit()
        # else:
        #     return sample_filter_df

    existing_image_id = [iid for iid in sample_filter_df['Image Data ID'] if len(mybids.get(ImageID=iid))>0]
    desired_downloaded_files = sample_filter_df[~sample_filter_df['Image Data ID'].isin(existing_image_id)]
        
    return desired_downloaded_files


# Create the dataset description. json
def create_dataset_description_json(dataset_name, main_folder_path):
    '''
    This function creates a .json file with the following fields:
        Name, BIDSVersion, License, Authors, Acknowledgements, HowToAcknowledge,
        Funding, ReferencesAndLinks and DatasetDOI

    Parameters
    ----------
    dataset_name : TYPE str
        DESCRIPTION. Name of the dataset that has been download from ADNI.
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.

    Returns
    -------
    None.

    '''
    data_descr = {}
    data_descr['Name'] = dataset_name
    data_descr['BIDSVersion'] = ''
    data_descr['License'] = ''
    data_descr['Authors'] = ['']
    data_descr['Acknowledgements'] = ''
    data_descr['HowToAcknowledge'] = ''
    data_descr['Funding'] = ['']
    data_descr['ReferencesAndLinks'] = ['']
    data_descr['DatasetDOI'] = ''
    data_descr_path = join(main_folder_path , 'dataset_description.json')

    json_object = json.dumps(data_descr, indent=4)
    with open(data_descr_path, 'w') as outfile:
        outfile.write(json_object)


def create_participants_tsv(main_folder_path, info_csv, mybids, adnimerge_df):
    '''
    

    Parameters
    ----------
    main_folder_path : TYPE
        DESCRIPTION.
    info_csv : TYPE
        DESCRIPTION.
    mybids : TYPE
        DESCRIPTION.
    adnimerge_df : TYPE
        DESCRIPTION.

    Returns
    -------
    participants_dict : TYPE
        DESCRIPTION.
    info_csv : TYPE
        DESCRIPTION.

    '''
    
    '''

    This function creates a .tsv file with the following column headers:
        'participant_id','SubjectID', 'AlternativeID', 'Group', 'Sex' and 'Age'

    Parameters
    ----------
    info_data : TYPE str
        DESCRIPTION. Path and name of the .csv file containing the downloaded ADNI dataset metadata.
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.

    Returns
    -------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    info_csv : TYPE DataFrame
        DESCRIPTION. Matrix containing all metadata of the downloaded ADNI dataset.
        It's columns are 'Data ID'',"Subject","Group","Sex","Age","Visit","Modality",
        "Description","Type","Acq Date","Format","Downloaded"
        This variable is created in 'crear_participants_tsv' function and is used in:
            Folder cration (Find sessions that we have images from to only create those)

    '''

    # CREATE DATASET VARS
    # participants.tsv
    col_names = ['participant_id', 'num_subject', 'adni_id', 'sex', 'ethnic', 'race', 'apoe4', 'education','group_bl', 'age_bl']


    write_header_flag = True
    participants_dict = {}
    existing_subjects = []
    if exists(join(main_folder_path, 'participants.tsv')):
        write_header_flag = False
        participants_dict = {p['participant_id']: p for p in mybids.get_file('participants.tsv').get_df().to_dict(orient='records')}
        existing_subjects = [p['adni_id'] for pid, p in participants_dict.items()]


    to_write = []
    subjs = info_csv['Subject'].unique()
    subjs_to_write = [s for s in subjs if s not in existing_subjects]
    for i, subject_id in enumerate(subjs_to_write):
        subject_bids_id = 'ADNI' + replace_underscore(subject_id)
        # subject_df = info_csv.loc[info_csv['Subject'] == subject_id]

        baseline_visit = np.argmin(adnimerge_df.loc[adnimerge_df['PTID'] == subject_id]['AGE'])
        metadata_bl = adnimerge_df.iloc[baseline_visit]
        d = {
            'participant_id': subject_bids_id,
            'num_subject': str(i + len(existing_subjects)).rjust(5, '0'),
            'adni_id': subject_id,
            'sex':  metadata_bl['PTGENDER'],
            'education': metadata_bl['PTEDUCAT'],
            'apoe4': metadata_bl['APOE4'],
            'ethnic': metadata_bl['PTETHCAT'],
            'race': metadata_bl['PTRACCAT'],
            'age_bl': -metadata_bl['AGE'],
            'group_bl': metadata_bl['DX'],
        }

        participants_dict[subject_bids_id] = d
        to_write.append(d)

    if to_write:
        with open(join(main_folder_path, 'participants.tsv'), 'a') as csv_file:
            tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
            if write_header_flag:
                tsv_writer.writeheader()
            tsv_writer.writerows(to_write)

    return participants_dict, info_csv


def create_participants_json(bids_filepath):
    '''
    This function creates a .json file with the following fields:
        Participant Id, Subject ID, AlternativeID, Group, Sex and Age.

    Parameters
    ----------
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.

    Returns
    -------
    None.

    '''

    # participants.json

    
    if exists( join(bids_filepath, 'participants.json')): 
        return
    
    data_descr={
        "participant_id":{
            "Description": "Participant identifier in BIDS format.",
            },
        "num_subject":{
            "Description": "Number of subject ordered by the inclusion in the BIDS dataset.",
            },
        "sex":{
            "Description": "sex of the participant as reported by the participant",
            "Levels": {
                "M": "male",
                "F": "female"
                    }
            },
        "group_bl":{
            "Description": 'Basleine diagnostic label.',
            "Levels": {
                "CN": '?'##### WHat are the levels?
                }
            },
        "age_bl":{
            "Description": "Age of the participant at the baseline visit",
            "Units": "years"
            },
        "adni_id ":{
            "Description": "Identification as downloaded from ADNI Datast"
            }
        }
    data_descr_path = join(bids_filepath, 'participants.json')
    json_object = json.dumps(data_descr, indent=4)
    with open(data_descr_path, 'w') as outfile:
        outfile.write(json_object)


def create_participants_file(bids_filepath, sample_filter_df, mybids, adnimerge_df):
    '''
    

    Parameters
    ----------
    bids_filepath : TYPE
        DESCRIPTION.
    sample_filter_df : TYPE
        DESCRIPTION.
    mybids : TYPE
        DESCRIPTION.
    adnimerge_df : TYPE
        DESCRIPTION.

    Returns
    -------
    subjs : TYPE
        DESCRIPTION.
    info_csv : TYPE
        DESCRIPTION.

    '''
    '''

    This function calls the two functions 'create_participants_tsv' and
    'create_participants_json' that create the .tsv and .json participants files

    Parameters
    ----------
    main_folder_path : TYPE str
        DESCRIPTION. Path to the main folder where all BIDS structure is contained.
    info_data : TYPE str
        DESCRIPTION. Path and name of the .csv file containing the downloaded ADNI dataset metadata.

    Returns
    -------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    info_csv : TYPE DataFrame
        DESCRIPTION. Matrix containing all metadata of the downloaded ADNI dataset.
        It's columns are 'Data ID'',"Subject","Group","Sex","Age","Visit","Modality",
        "Description","Type","Acq Date","Format","Downloaded"
        This variable is created in 'crear_participants_tsv' function and is used in:
            Folder cration (Find sessions that we have images from to only create those)

    '''
    subjs, info_csv = create_participants_tsv(bids_filepath, sample_filter_df, mybids, adnimerge_df)
    create_participants_json(bids_filepath)

    return subjs, info_csv
    # Subjs is used in subjects Folder, subject var
    # (subject_id_sessions.tsv) is used in creating BIDS structure folders (Find sessions that we have images
    # from to only create those)

    # info_csv is used in creating BIDS structure folders (Find sessions that we have images
    # from to only create those)


# Mirar que no hi hagin participants repetits. Ho veiem a traves de les sessions
# i a través de la imatge concreta.

# CREATE SUBJECT FOLDERS
def create_subject_folders(subject_dict, main_folder_path):
    '''
    

    Parameters
    ----------
    subject_dict : TYPE
        DESCRIPTION.
    main_folder_path : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    '''
    This function creates all subject folders following the BIDS guidelines.

    Parameters
    ----------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    main_folder_path : TYPE str
       DESCRIPTION. Path to the main folder where all BIDS structure is contained.
    subjects_prefix : TYPE str
        DESCRIPTION.Text prefix with which each subject file has been named 'sub-'
    ID_list : TYPE Numpy Object array
        DESCRIPTION. List of all image data ID. Set by the instruction:
            ID_list = df_collec['Image Data ID'].unique()

    Returns
    -------
    None.
    it adds the 'Alternative_ID' column in df_collec variable.

    '''
    for sub_id in subject_dict.keys():
        if not isdir(join(main_folder_path, 'sub-' + sub_id)):
            makedirs(join(main_folder_path, 'sub-' + sub_id))


# CREATE SUBJECT VARS
# subject_id_sessions.tsv
def subject_id_sessions_tsv(remaining_files_df, participants_dict, main_folder_path, mybids, adnimerge_df):
    '''
    

    Parameters
    ----------
    remaining_files_df : TYPE
        DESCRIPTION.
    participants_dict : TYPE
        DESCRIPTION.
    main_folder_path : TYPE
        DESCRIPTION.
    mybids : TYPE
        DESCRIPTION.
    adnimerge_df : TYPE
        DESCRIPTION.

    Returns
    -------
    remaining_files_df : TYPE
        DESCRIPTION.

    '''
    '''
    This function creates .tsv file that contains the metadata of all sessions

    Parameters
    ----------
    subjs : TYPE NumPy Object array
        DESCRIPTION. List containing the participant Id of each subject in the .csv file pointed in
        the variable 'info_data' Example: 003_S_6264
        This variable is created in 'crear_participants_tsv' function and is used in:
            subjects Folder, subject var (subject_id_sessions.tsv)
            and folder cration (Find sessions that we have images from to only create those)
    main_folder_path : TYPE str
       DESCRIPTION. Path to the main folder where all BIDS structure is contained.
    subjects_prefix : TYPE str
        DESCRIPTION.Text prefix with which each subject file has been named 'sub-'
    visits_df : TYPE dataframe
        DESCRIPTION. Dictionary containing the visits in the ADNI protocol and their code.
        Columns: 	Phase	ID	VISCODE	VISNAME	VISORDER
    idas_df : TYPE dataframe
        DESCRIPTION. Dataframe containing the information of the 'idasearch.csv' file, downloaded from ida webpage
        It contains the information on all subjects, subject's sessions and their description
        Columns: Subject_ID, Projet, Phase, Sex, Weight, Research_Group, APOE A1, APOE A2, Visit, Study Date, Archive Date, Age, Modality, Description, Imaging Protocol.
    dataset_name : TYPE str
        DESCRIPTION. Name given by the ADNI database once downloaded. in this case is "\ADNI_50.0-59.0"
    mrilist_df : TYPE dataframe
        DESCRIPTION. List of all MRI images and their description.
        Columns: Unnamed: 0(Index), TYPE, SUBJECT, VISIT, MAGSTRENGTH, SEQUENCE, SCANDATE, STUDYID, IMAGEUID, Visit_Month_Number
    petlist_df : TYPE dataframe
       DESCRIPTION. List of all PET images and their description.
       Columns: Unnamed: 0(Index), TYPE, SUBJECT, VISIT, MAGSTRENGTH, SEQUENCE, SCANDATE, STUDYID, IMAGEUID, Visit_Month_Number

    Returns
    -------
    None.

    '''

    col_names = ['session_id', 'session_description', 'age', 'time_to_bl_days', 'time_to_bl_months', 'time_to_bl_years',
                 'dx', 'phase', 'site', 'viscode', 'icv', 'fdg', 'abeta', 'ptau', 'ttau', 'mmse', 'cdrsb', 'adas11',
                 'adas13', 'faq']

    for i, (subject_id, subject_dict) in enumerate(participants_dict.items()):

        if subject_id == 'ADNI168S6085':
            pass
        sess_desc_path_json = join(main_folder_path, 'sub-' + subject_id,  'sub-' + subject_id + '_sessions.json')
        sess_desc_path_tsv = join(main_folder_path, 'sub-' + subject_id,  'sub-' + subject_id + '_sessions.tsv')
        
        # adni_project = ((idas_df[idas_df['Subject ID'] == subject_dict['adni_id']]['Phase']).unique()[0]).replace(' ', '')
        # viscodes_project = visits_df[visits_df['Phase'] == adni_project]['VISCODE']
        # visits_project = visits_df[visits_df['Phase'] == adni_project]['VISNAME']
        # visits_phase = visits_df[visits_df['Phase'] == adni_project]
        write_header_flag = False
        if not exists(sess_desc_path_tsv): write_header_flag = True
        if not exists(sess_desc_path_json):
            data_descr = {}
            data_descr['session_id'] = 'session identifier'
            data_descr['session_description'] = 'session description as in VISITS.csv file from the ADNI database'
            data_descr['age'] = ''
            data_descr['time_to_bl_days'] = ''
            data_descr['time_to_bl_months'] = ''
            data_descr['time_to_bl_years'] = ''
            data_descr['dx'] = ['']
            data_descr['phase'] = ''
            data_descr['site'] = ''
            data_descr['viscode'] = ''
            data_descr['icv'] = ''
            data_descr['fdg'] = ''
            data_descr['abeta'] = ''
            data_descr['ptau'] = ''
            data_descr['ttau'] = ''
            data_descr['mmse'] = ''
            data_descr['cdrsb'] = ''
            data_descr['adas11'] = ''
            data_descr['adas13'] = ''
            data_descr['faq'] = ''


            json_object = json.dumps(data_descr, indent=4)
            with open(sess_desc_path_json, 'w') as outfile:
                outfile.write(json_object)

        # if exists(sess_desc_path_tsv): continue
        subject_df = remaining_files_df[remaining_files_df['Subject'] == subject_dict['adni_id']]
        sessions_unique = subject_df['session_id'].unique()

        existing_sessions = []
        sessions_df = mybids.get(subject=subject_id, extension='tsv', regex_search='sessions')
        if sessions_df:
            sessions_df = sessions_df[0].get_df()
            existing_sessions = sessions_df['session_id'].to_list()

        sess_to_write = [s for s in sessions_unique if s not in existing_sessions]

        session_list = []
        for sess_viscode in sess_to_write:
            sess_description = subject_df['Visit'].unique()#visits_phase[visits_phase['VISCODE'] == sess_viscode]['VISNAME'].tolist()[0]

            # To remove: now these folders are created just before the DCM2NII if they do not exists.
            # sess_id = 'ses-' + sess_viscode
            # sess_df = subject_df[subject_df['session_id'] == sess_viscode]
            # sess_description = sess_description.replace(' ', '')
            #
            # for _, image_df in sess_df.iterrows():
            #     modality = image_df['Modality']
            #     if not exists(join(main_folder_path, 'sub-' + subject_id, sess_id, dict_modalities[modality])):
            #         makedirs(join(main_folder_path, 'sub-' + subject_id, sess_id, dict_modalities[modality]))
            session_adnimerge = adnimerge_df.loc[(adnimerge_df['PTID'] == subject_dict['adni_id']) & (adnimerge_df['Month'] == int(sess_viscode.split('m')[1]))]
            if len(session_adnimerge) > 1:
                sess_metadata = session_adnimerge.iloc[0]
                d = {
                    'session_id': sess_viscode,
                    'session_description': '_'.join(sess_description),
                    'age': sess_metadata['AGE'],
                    'time_to_bl_days': sess_metadata['Years_bl'] * 365.25,
                    'time_to_bl_months': sess_metadata['Month_bl'],
                    'time_to_bl_years': sess_metadata['Years_bl'],
                    'dx': sess_metadata['DX'],
                    'phase': sess_metadata['COLPROT'],
                    'site': sess_metadata['SITE'],
                    'viscode': sess_metadata['VISCODE'],
                    'icv': sess_metadata['ICV'],
                    'mmse': sess_metadata['MMSE'],
                    'fdg': sess_metadata['FDG'],
                    'abeta': sess_metadata['ABETA'],
                    'ptau': sess_metadata['PTAU'],
                    'ttau': sess_metadata['TAU'],
                    'cdrsb': sess_metadata['CDRSB'],
                    'adas11': sess_metadata['ADAS11'],
                    'adas13': sess_metadata['ADAS13'],
                    'faq': sess_metadata['FAQ'],
                }
            else:
                d = {
                    'session_id': sess_viscode,
                    'session_description': '_'.join(sess_description),
                    'age': '',
                    'time_to_bl_days': '',
                    'time_to_bl_months': '',
                    'time_to_bl_years': '',
                    'dx': '',
                    'phase': '',
                    'site': '',
                    'viscode': '',
                    'icv': '',
                    'mmse': '',
                    'fdg': '',
                    'abeta': '',
                    'ptau': '',
                    'ttau': '',
                    'cdrsb': '',
                    'adas11': '',
                    'adas13': '',
                    'faq': '',
                }

            session_list.append(d)

        if len(session_list) > 0:
            with open(sess_desc_path_tsv, 'a') as csv_file:
                tsv_writer = csv.DictWriter(csv_file, fieldnames=col_names, delimiter='\t')
                if write_header_flag:
                    tsv_writer.writeheader()
                tsv_writer.writerows(session_list)

    return remaining_files_df

# Convertir de DCM a NIFTI
def convert_DCM2NIFTI(sample_filter_df, pet_df, mri_df):
    '''
    

    Parameters
    ----------
    sample_filter_df : TYPE
        DESCRIPTION.
    pet_df : TYPE
        DESCRIPTION.
    mri_df : TYPE
        DESCRIPTION.

    Returns
    -------
    failed_images : TYPE
        DESCRIPTION.

    '''
    '''
    This function converts all the downloaded images from Dicom format to NIFTI format.

    Parameters
    ----------
    df_collec : TYPE dataframe
        DESCRIPTION. Dataframe containing all image with its corresponding data. It is downloaded from Ida Loni - Adni webpage along with the actual image data and the idaSearch_3_21_2022.csv file.
        Column Fields: Image Data ID,"Subject","Group","Sex","Age","Visit","Modality","Description","Type","Acq Date","Format","Downloaded", Alternative_ID, Dicom_path, nifti_path, session

    converter : TYPE interfaces.dcm2nii.Dcm2niix
        DESCRIPTION. Dcm2niix object of nipype.interfaces.dcm2nii.Dcm2niix module
        converter = Dcm2niix()

    Returns
    -------
    None.

    '''

    # TODO: Rename + add stuff to the json (for example, the ADNI image ID).
    # Descargar librarias BIDS =====> conda install pybids not working in windows.
    failed_images = []
    for _, image_df in sample_filter_df.iterrows():
        dicom_path = image_df['dicom_path']
        nifti_path = image_df['nifti_path']
        json_path = image_df['json_path']

        if not exists(dirname(nifti_path)):
            makedirs(dirname(nifti_path))

        if os.path.isdir(dicom_path) and not exists(nifti_path):
            subprocess.run(["dcm2niix", "-z", "y", "-o", dirname(nifti_path), dicom_path])

            # Rename .nii.gz and .json files ------>
            rename_image(nifti_path, json_path, image_df['Image Data ID'])

        if exists(json_path):
            # Add stuff to .json file
            fill_image_json(image_df, json_path, pet_df, mri_df)
        else:
            failed_images.append(dicom_path)


    return failed_images
# C:\Users\Usuari\Documents\Master_Biomedic\TFM\Proves\Prova_2_Roser C:/Users/Usuari/Documents/Master_Biomedic/TFM/Proves/Prova_2_Roser

# Execute functions that create the .json and .tsv files needed in the participants file
# ==============================================================================
# =================            UTILS           =================================
# ==============================================================================

###
#Sujeto, session, tipo adq.

def replace_underscore(word):
    '''
    

    Parameters
    ----------
    word : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return word.replace('_', '')

def rename_image(nifti_filepath, json_filepath, image_id):
    '''
    

    Parameters
    ----------
    nifti_filepath : TYPE str
        DESCRIPTION. 
    json_filepath : TYPE str
        DESCRIPTION.
    image_id : TYPE str
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    '''
    This function gets execute right after nifti files are created. It changes 
    its name to a BIDS accepted name.

    Parameters
    ----------
    image_df : TYPE List
        DESCRIPTION. Row of sample_filter_df DataFrame as a list. 
    nifti_file : TYPE str
        DESCRIPTION. BIDS accepted name for nifti file.
    json_file : TYPE str
        DESCRIPTION. BIDS accepted name for json file.

    Returns
    -------
    None.

    '''
    images_in_folder = os.listdir(dirname(nifti_filepath))
    for image in images_in_folder:
        if image_id not in str(image):
            continue

        elif 'json' in image:
            os.rename(join(dirname(nifti_filepath), image), json_filepath)

        elif 'nii.gz' in image:
            os.rename(join(dirname(nifti_filepath), image), nifti_filepath)

        else:
            print('There is a non .json and non .nii.gz file in the following directory: '+ nifti_filepath)
                
def fill_image_json(image_df, json_filepath, pet_df, mri_df):
    '''
    

    Parameters
    ----------
    image_df : TYPE List
        DESCRIPTION. Row of sample_filter_df DataFrame as a list
    json_filepath : TYPE
        DESCRIPTION.
    pet_df : TYPE
        DESCRIPTION.
    mri_df : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # import json
    # Llegir JSON:
    json_file = open(json_filepath, 'r')
    json_dict = json.load(json_file)
    json_dict['ImageID'] = image_df['Image Data ID']

    # Escriure JSON
    if image_df['Modality'].lower()=='pet':
        json_dict['TracerName'] = 'TracerName'
        json_dict['TracerRadionuclide'] = 'TracerRadionuclide'
        json_dict['ModeOfAdministration'] = 'ModeOfAdministration'
        json_dict['InjectedRadioactivity'] = 666
        json_dict['InjectedRadioactivityUnits'] = 'MBq'
        json_dict['InjectedMass'] = 666
        json_dict['InjectedMassUnits'] = 'ug'
        json_dict['SpecificRadioactivity'] = 666
        json_dict['SpecificRadioactivityUnits'] = 'MBq/ug'
        json_dict['TimeZero'] = '00:00:00'
        json_dict['ScanStart'] = 0
        json_dict['InjectionStart'] = 0
        json_dict['FrameTimesStart'] = [0, 10, 20, 30, 40, 50, 60, 80, 100, 120, 140, 160, 180, 240, 300, 360, 420, 480, 540, 660, 780, 900, 1020]
        json_dict['AcquisitionMode'] = 'list mode'
        json_dict['ImageDecayCorrected'] = True
        json_dict['ImageDecayCorrectionTime'] = 0
        json_dict['ReconMethodName'] = '"3D-OSEM-PSF'
        json_dict['ReconMethodParameterLabels'] = ["subsets","iterations"]
        json_dict['ReconMethodParameterUnits'] = ["none","none"]
        json_dict['ReconMethodParameterValues'] = [666, 666]
        json_dict['ReconFilterType'] = 'none'
        json_dict['ReconFilterSize'] = 666
        json_dict['AttenuationCorrection'] = "[137Cs]transmission scan-based"
        json_dict['InfusionRadioactivity'] = 666
        json_dict['InfusionStart'] = 666
        json_dict['InfusionSpeed'] = 666
        json_dict['InfusionSpeedUnits'] = 'InfusionSpeedUnits'
        json_dict['InjectedVolume'] = 666
        if len(pet_df.loc[pet_df['Image ID'] == int(json_dict['ImageID'][1:])]) > 0:
            json_dict['ScanDate'] = pet_df.loc[pet_df['Image ID'] == int(json_dict['ImageID'][1:])].iloc[0]['Scan Date']
        else:
            json_dict['ScanDate']  = ''
        json_object = json.dumps(json_dict, indent=4)
        json_file = open(json_filepath, "w")
        json_file.write(json_object)
        
    elif image_df['Modality'].lower()=='fmri':
        json_dict['TaskName'] = 'Resting State'
        if len(mri_df.loc[mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])]) > 0:
            json_dict['ScanDate'] = mri_df.loc[mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])].iloc[0]['SCANDATE']
        else:
            json_dict['ScanDate']  = ''

        json_object = json.dumps(json_dict, indent=4)
        # json_absolute_path = str(nifti_path)+'/'+str(json_file)
        json_file = open(json_filepath, "w")
        json_file.write(json_object)
    
    elif image_df['Modality'].lower()=='mri':
        #ADNI identificador (I########)-----> d'aquí hauriem de veure si la imatge nova ja hi és a la base de dades.
        #Data d'aquisició
        if len(mri_df.loc[mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])]) > 0:
            json_dict['ScanDate'] = mri_df.loc[mri_df['IMAGEUID'] == int(json_dict['ImageID'][1:])].iloc[0]['SCANDATE']
        else:
            json_dict['ScanDate'] = ''
        json_object = json.dumps(json_dict, indent=4)
        # json_absolute_path = str(nifti_path)+'/'+str(json_file)
        json_file = open(json_filepath, "w")
        json_file.write(json_object)

    else:
        print('Image '+image_df['Image Data ID']+ ' unwanted image modality ('+image_df['Modality']+') in folder. Subject:'+image_df['subject_id'])
    

#
# # TODO implement this function
# def get_PET_contrast(image_df, dict_file_names):
#     if 'av45' in image_df['Description'].lower() or 'florbetapir' in image_df['Description'].lower():
#         sufix = 'av45'
#
#     elif 'fdg' in image_df['Description'].lower():
#         sufix = 'fdg'
#
#     elif 'tau' in image_df['Description'].lower():
#         sufix = 'tau'
#
#     elif 'dynamic' in image_df['Description'].lower():
#         if '6x5min' in image_df['Description'].lower():
#             sufix = 'dynamic6x5min'
#         elif '4x5min' in image_df['Description'].lower():
#             sufix = 'dynamic4x5min'
#         else:
#             sufix = dict_file_names['PET'] + 'dynamic'
#
#     else:
#         sufix = None
#
#     return sufix
#
# def get_suffix(image_df, dict_file_names):
#     '''
#     This function is used inside the convert_DCM2NIFTI function before the actual format convertion. It stablishes
#     the suffix of an image file given it's modality.
#
#
#     Parameters
#     ----------
#     image_df : TYPE List
#         DESCRIPTION. Row of sample_filter_df DataFrame as a list.
#     dict_file_names : TYPE Dict
#         DESCRIPTION. BIDS accepted modality names for PET, fMRI and MRI
#
#     Returns
#     -------
#     sufix : TYPE str
#         DESCRIPTION. BIDS accepted suffix for each image given it's modality.
#
#     '''
#
#     if 'fmri' in image_df['Modality'].lower():
#
#         if 'fmri' in image_df['Description'].lower():
#
#             sufix = dict_file_names['fMRI']#'task-rest_' + dict_file_names['fMRI']
#
#         else:
#             sufix = '----'
#     #
#     #             rename(name_func, i, df_collec, sufix)
#     #
#     elif 'mri' in image_df['Modality'].lower():
#
#         if 'sense2' in image_df['Description'].lower():
#             sufix = 'T1w'
#         elif 'ts_2' in image_df['Description'].lower():
#             sufix = '_TS_2'
#         elif 'ts_3' in image_df['Description'].lower():
#             sufix = '_TS_3'
#         else:
#             sufix = 'T1w'
#
#     elif 'pet' in image_df['Modality'].lower():
#         sufix = dict_file_names['PET']
#     else:
#         sufix = 'none'
#
#     return sufix
#

