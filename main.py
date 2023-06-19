from fn_utils import *
import bids

# Input variables
subjects_prefix = 'sub-'
path_data = '/home/acasamitjana/data/'
dataset_name = 'ADNI55-60'
sample_filename = 'def_sample-55-60y_1_30_2023.csv'
# idasearch_filename = 'idaSearch_1_19_2023.csv'

bids_filepath = join(path_data, 'BIDS')
sample_filepath = join(path_data, 'util-files', sample_filename)
# idasearch_filepath = join(path_data, dataset_name, idasearch_filename)
mrilist_filepath = join(path_data, 'util-files', 'MRILIST-with-Visit.csv')
petlist_filepath = join(path_data, 'util-files', 'PET_META_LIST-with-Visit.csv')
visits_filepath = join(path_data, 'util-files', 'VISITS.csv')
adnimerge_filepath = join(path_data, 'util-files', 'ADNIMERGE.csv')

sample_df = pd.read_csv(sample_filepath)
# idas_df = pd.read_csv(idasearch_filepath)
mrilist_df = pd.read_csv(mrilist_filepath)
petlist_df = pd.read_csv(petlist_filepath)
visits_df = pd.read_csv(visits_filepath)
adnimerge_df = pd.read_csv(adnimerge_filepath, low_memory=False)

# To write
data_info_to_put_in_readme = 'This dataset consists of a collection of 3 subjects downloaded from ADNI public dataset. \
    \nThese subjects were selected so as to have all 3 image modalities: MRI (T1w), fMRI (resting state), PET (AV45). \
    \nThis dataset was dowload to test some functions for conversion, preprocessing...'
data_info_to_put_in_changes = '2022-02-16 \
    \nFirst Creationg of BIDS Structure'



# dict_modalities = {'PET': 'pet', 'fMRI': 'func', 'MRI': 'anat'}
dict_file_names = {
    # 'PET': 'pet',
    'fMRI': {'suffix': 'bold', 'task': 'rest'},
    'fcMRI': {'suffix': 'bold', 'task': 'rest'},
    # 'Hipp': 't2hipp',
    'MPRAGE': {'suffix': 'T1w'},
    'MPR': {'suffix': 'T1w'},
    'MP-': {'suffix': 'T1w'},
    'IR-SPGR': {'suffix': 'T1w'},
    'IR-FSPGR': {'suffix': 'T1w'},
    'AV45': {'trc': 'fbp', 'suffix': 'pet'},
    'florbetapir': {'trc': 'fbp', 'suffix': 'pet'},
    'FBB': {'trc': 'fbb', 'suffix': 'pet'},
    'Florbetaben': {'trc': 'fbb', 'suffix': 'pet'},
    'PIB': {'trc': 'pib', 'suffix': 'pet'},
    'FLAIR': {'suffix': 'FLAIR'},
    'T2': {'suffix': 'T2w'},
}

id_list = sample_df['Image Data ID'].unique()
subject_list = sample_df['Subject'].unique()



########-------- START CREATING STRUCTURE --------########
# Create BIDS main folder and read the dataset
if not isdir(bids_filepath): makedirs(bids_filepath)

# create README
readme_path = join(bids_filepath, 'README')
if not os.path.isfile(readme_path):
    file_obj = open(readme_path, 'a')
    file_obj.write(data_info_to_put_in_readme)
    file_obj.close()

# CHANGES
changes_path = join(bids_filepath, 'CHANGES')
if not os.path.isfile(changes_path):
    file_ch_obj = open(changes_path, 'a')
    file_ch_obj.write(data_info_to_put_in_changes)
    file_ch_obj.close()

# This function creates the dataset description.json file inside the main folder.
create_dataset_description_json(dataset_name, bids_filepath)
# It works correctly

################
# Let's start! #
################
mybids = bids.layout.BIDSLayout(root=bids_filepath)#, validate=False)

# This function removes all non necessary modalities and returns a DataFrae listing all images in the current working dataset: 'dataset_name'
sample_filter_df, excluded_modalities = pick_desired_modalites(path_data, dataset_name, sample_df, mybids, dict_file_names, mri_df=mrilist_df, pet_df=petlist_df, adnimerge_df=adnimerge_df)
# It works correctly

# This function remove files that altready exist
sample_filter_df = remove_repeated_files(sample_filter_df, mybids)
# It works correctly


# This function creates the participants.json and a participants.tsv files inside the main folder.
participants_dict, info_csv = create_participants_file(bids_filepath, sample_filter_df, mybids, adnimerge_df)
# It works correctly



# This functin creates a folder for each subject inside the main folder.
create_subject_folders(participants_dict, bids_filepath)
# It works correctly
# Mirar en dataframe creado si el sujeto a evaluar está ya existente en la base de datos.


# This function creates a .json and a .tsv where sessiona are to be described inside each participant folder.
sample_filter_df = subject_id_sessions_tsv(sample_filter_df, participants_dict, bids_filepath, mybids, adnimerge_df)


# This functions takes all images from the Downloaded dataset, converts it from DICOM to NIFTI and places it inside each session folder.
failed_images = convert_DCM2NIFTI(sample_filter_df, petlist_df, mrilist_df)

print('\n\n\n####################')
print('##### Summary ######')
print('####################')

print('Total excluded modalities: ' + str(len(excluded_modalities)))
print('\n'.join(excluded_modalities))

print('\n Total Images: ' + str(len(sample_df)))
print('\n Total Filtered Images: ' + str(len(sample_filter_df)))
print('\n Total Converted Images: ' + str(len(sample_filter_df) - len(failed_images)))
print('\n Total Failed Images: ' + str(len(failed_images)))
print('\n -'.join(failed_images))