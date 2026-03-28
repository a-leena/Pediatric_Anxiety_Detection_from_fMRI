from pathlib import Path
import os
from glob import glob
import json
import numpy as np
import pandas as pd
from ast import literal_eval
import itertools
import matplotlib.pyplot as plt
from nilearn.image import load_img, iter_img, math_img, mean_img, new_img_like, resample_to_img, concat_imgs
from nilearn.plotting import (plot_img, plot_stat_map, plot_design_matrix, 
plot_contrast_matrix, plot_roi, plot_connectome, show, find_parcellation_cut_coords)
from nilearn.masking import compute_epi_mask, unmask
from nilearn.glm.first_level import FirstLevelModel
from nilearn.glm.second_level import SecondLevelModel, non_parametric_inference
from nilearn.glm import threshold_stats_img
from nilearn.reporting import get_clusters_table
from nilearn.maskers import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nilearn import datasets
import warnings
warnings.filterwarnings("ignore")

# Paths
BASE = '..'
DATASET_PATH = os.path.join(BASE, 'Pediatric_Anxiety_Disorder')
HOLDOUTS_PATH = os.path.join(BASE, 'dataset_holdouts')
RESAMPLED_ATLAS_PATH = os.path.join(BASE,'resampled_atlases')

# WM & CSF Threshold
WM_CSF_PROB_THRESH = 0.95

# GLM & Contrast Testing HYPERPARAMETERS
HRF_MODEL = 'spm'
DRIFT_MODEL = 'cosine'
HIGH_PASS = 0.01
LOW_PASS = 0.08
STANDARDIZE = 'zscore_sample'
SMOOTHING = 6

# FMRI PROPERTIES
TASK_TR = 2.3
REST_TR = 2
SLICE_TIME_REF = 0.0
PREP_RESOLUTION = 2

# TEs for Cohort 2
TEs = [14.8, 28.4, 42.0]

# CONTRASTS
CONTRASTS = {
    "congruent_effect": "congruent",
    "incongruent_effect": "incongruent",
    "neutral_effect": "neutral",
    "incongruent_vs_congruent": "incongruent - congruent",
    "incongruent_vs_neutral": "incongruent - neutral",
    "neutral_vs_congruent": "neutral - congruent",
    "task_vs_baseline": "(congruent + incongruent + neutral) / 3",
    "error_monitoring": "error - (congruent + incongruent + neutral) / 3"
}

def save_output_file(file, path, filename):
    os.makedirs(path, exist_ok=True)
    filepath = os.path.join(path, filename)
    if filename.endswith('csv'):
        file.to_csv(filepath, index=False)
    elif filename.endswith('nii') or filename.endswith('nii.gz'):
        file.to_filename(filepath)
    elif filename.endswith('npy'):
        np.save(filepath, file)
    elif filename.endswith('json'):
        with open(filepath, 'w') as f:
                json.dump(file, f, indent=4)
    return filepath

def load_file(path, filename):
    filepath = os.path.join(path, filename)
    # print(f"Loading {filepath}")
    if filename.endswith('csv'):
        return pd.read_csv(filepath)
    if filename.endswith('nii') or filename.endswith('nii.gz'):
        return load_img(filepath)
    if filename.endswith('npy'):
        return np.load(filepath)
    if filename.endswith('json'):
        with open(filepath, 'r') as file:
            return json.load(file)

def safe_parse(val):
    if isinstance(val, str):
        return literal_eval(val)
    return val
    
def get_sub_info_list(sub_list, phenotype_table):
    subjects_info = []
    for subject in sub_list:
        sub_id = 'sub-0'+str(subject)
        sub_info = {
            'subject_id': str(subject),
            'cohort': 1 if phenotype_table[phenotype_table['participant_id']==sub_id]['COHORT'].values[0]==1 else 2,
            'diagnosis': 'Anxiety' if phenotype_table[phenotype_table['participant_id']==sub_id]['KSADS_MAIN_DIAGNOSIS'].values[0]== 'ANX' else 'Healthy'
        }
        subjects_info.append(sub_info)
    return subjects_info
    
def get_maps_tables(sub_info, get_task=True, get_task_mean=True, get_rest=True, get_probmaps=True, get_anat=False):
    # Paths
    subject_path = os.path.join(DATASET_PATH, 'sub-'+sub_info['subject_id'], 'ses-1')
    func_path = os.path.join(subject_path, 'func')
    if get_task:
        run1_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-1_bold.nii')
        run2_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-2_bold.nii')
        run1_events_path = os.path.join(func_path, f'sub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-1_events.tsv')
        run2_events_path = os.path.join(func_path, f'sub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-2_events.tsv')
        run1_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-1_bold.txt')
        run2_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-2_bold.txt')
    if get_task_mean:
        mean_path = os.path.join(func_path, f'wmeanasub-{sub_info['subject_id']}_ses-1_task-TAU{sub_info['cohort']}_run-1_bold.nii')
    if get_rest:
        if sub_info['cohort']==1:
            rest_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-rest_bold.nii')
            rest_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-rest_bold.txt')
        else:
            rest_echo1_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-rest_echo-1_bold.nii')
            rest_echo1_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-rest_echo-1_bold.txt')
            rest_echo2_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-rest_echo-2_bold.nii')
            rest_echo2_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-rest_echo-2_bold.txt')
            rest_echo3_path = os.path.join(func_path, f'wrasub-{sub_info['subject_id']}_ses-1_task-rest_echo-3_bold.nii')
            rest_echo3_confounds_path = os.path.join(func_path, f'rp_asub-{sub_info['subject_id']}_ses-1_task-rest_echo-3_bold.txt')
    if get_probmaps:
        probmaps_path = os.path.join(subject_path, 'anat')
        wm_probmap_path = os.path.join(probmaps_path, f'c2sub-{sub_info['subject_id']}_ses-1_T1w.nii')
        csf_probmap_path = os.path.join(probmaps_path, f'c3sub-{sub_info['subject_id']}_ses-1_T1w.nii')
    if get_anat:
        anat_path = os.path.join(subject_path, 'anat', f'sub-{sub_info['subject_id']}_ses-1_T1w.nii')
    
    # Loading images
    if get_task:
        run1 = load_img(run1_path)
        run2 = load_img(run2_path)
    if get_task_mean:
        mean_fmri = load_img(mean_path)
    if get_rest:
        if sub_info['cohort']==1:
            rest = load_img(rest_path)
        else:
            rest_echo1 = load_img(rest_echo1_path)
            rest_echo2 = load_img(rest_echo2_path)
            rest_echo3 = load_img(rest_echo3_path)
            rest = math_img(
                f"({TEs[0]}*img1 + {TEs[1]}*img2 + {TEs[2]}*img3)/{sum(TEs)}",
                img1=rest_echo1, img2=rest_echo2, img3=rest_echo3
            )
    if get_probmaps:
        wm_probmap = load_img(wm_probmap_path)
        csf_probmap = load_img(csf_probmap_path)
    if get_anat:
        anat = load_img(anat_path)
    
    # Loading tables
    if get_task:
        run1_events_table = pd.read_table(run1_events_path)[['onset','duration','trial_type']]
        run2_events_table = pd.read_table(run2_events_path)[['onset','duration','trial_type']]
        run1_confounds_table = pd.read_table(run1_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')
        run2_confounds_table = pd.read_table(run2_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')
    if get_rest:
        if sub_info['cohort']==1:
            rest_confounds_table = pd.read_table(rest_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')
        else:
            rest_confounds_table = pd.read_table(rest_echo1_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')
            # rest_echo2_confounds_table = pd.read_table(rest_echo2_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')
            # rest_echo3_confounds_table = pd.read_table(rest_echo3_confounds_path, header=None, names=['tx', 'ty', 'tz', 'rx', 'ry', 'rz'], sep='\\s+')

    # Returning results
    if get_task and not get_task_mean and not get_rest and not get_anat:
        return run1, run2, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table
    if get_task and get_task_mean and not get_rest and not get_probmaps and not get_anat:
        return run1, run2, mean_fmri, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table
    if get_task and get_task_mean and not get_rest and get_probmaps and not get_anat:
        return run1, run2, mean_fmri, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, wm_probmap, csf_probmab
    if get_task and get_task_mean and get_rest and not get_probmaps and not get_anat:
        return run1, run2, mean_fmri, rest, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, rest_confounds_table
    if get_task and get_task_mean and get_rest and get_probmaps and not get_anat:
        return run1, run2, mean_fmri, rest, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, rest_confounds_table, wm_probmap, csf_probmab
    if get_task and get_task_mean and get_rest and not get_probmaps and get_anat:
        return run1, run2, mean_fmri, rest, anat, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, rest_confounds_table
    if get_task and get_task_mean and get_rest and get_probmaps and get_anat:
        return run1, run2, mean_fmri, rest, anat, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, rest_confounds_table, wm_probmap, csf_probmab
    if get_task_mean and not get_rest and not get_anat:
        return mean_fmri
    if get_task_mean and get_rest and not get_probmaps and not get_anat:
        return mean_fmri, rest, rest_confounds_table
    if get_task_mean and get_rest and get_probmaps and not get_anat:
        return mean_fmri, rest, rest_confounds_table, wm_probmap, csf_probmab
    if get_task_mean and get_rest and not get_probmaps and get_anat:
        return mean_fmri, rest, anat, rest_confounds_table
    if get_task_mean and get_rest and get_probmaps and get_anat:
        return mean_fmri, rest, anat, rest_confounds_table, wm_probmap, csf_probmab
    if get_rest and not get_anat:
        return rest, rest_confounds_table
    if get_rest and get_anat:
        return rest, anat, rest_confounds_table
    if get_anat:
        return anat
    if get_probmaps:
        return wm_probmap, csf_probmab

def get_confounds(wm_probmap, csf_probmap, mean_fmri, for_task=True, for_rest=True, other_task_confounds_dfs=[], other_rest_confounds_df=None, task_fmri_runs=[], rest_fmri=None, threshold=0.9):
    # create masks for WM & CSF
    wm_mask = math_img(f"img >= {threshold}", img=wm_probmap)
    csf_mask = math_img(f"img >= {threshold}", img=csf_probmap)
    wm_mask_resampled = resample_to_img(
        wm_mask, mean_fmri, interpolation='nearest'
    )
    csf_mask_resampled = resample_to_img(
        csf_mask, mean_fmri, interpolation='nearest'
    )

    # obtain voxel-wise WM & CSF signals as 2D matrices
    if for_task and len(task_fmri_runs)>0:
        wm_masker = NiftiMasker(mask_img=wm_mask_resampled, standardize=False)
        csf_masker = NiftiMasker(mask_img=csf_mask_resampled, standardize=False)
        task_wm_signals, task_csf_signals = [], []
        for task_fmri_run in task_fmri_runs:
            task_wm_signals.append(wm_masker.fit_transform(task_fmri_run))
            task_csf_signals.append(csf_masker.fit_transform(task_fmri_run))
    if for_rest and rest_fmri is not None:
        wm_masker = NiftiMasker(mask_img=wm_mask_resampled, standardize=False)
        csf_masker = NiftiMasker(mask_img=csf_mask_resampled, standardize=False)
        rest_wm_signal = wm_masker.fit_transform(rest_fmri)
        rest_csf_signal = csf_masker.fit_transform(rest_fmri)

    # reducing WM & CSF signals to 1D
    if for_task:
        task_wm_means, task_csf_means = [], []
        for task_wm_signal, task_csf_signal in zip(task_wm_signals, task_csf_signals):
            task_wm_means.append(task_wm_signal.mean(axis=1))
            task_csf_means.append(task_csf_signal.mean(axis=1))
    if for_rest:
        rest_wm_mean = rest_wm_signal.mean(axis=1)
        rest_csf_mean = rest_csf_signal.mean(axis=1)

    # creating dataframe for WM & CSF signals combined with other confounds
    if for_task:
        task_confounds_dfs = []
        for i, (task_wm_mean, task_csf_mean) in enumerate(zip(task_wm_means, task_csf_means)):
            task_wm_csf_df = pd.DataFrame({
                'wm_signal': task_wm_mean,
                'csf_signal': task_csf_mean
            })
            if len(other_task_confounds_dfs)>0:
                if len(other_task_confounds_dfs)==len(task_wm_means):
                    task_confounds_df = pd.concat(
                        [other_task_confounds_dfs[i], task_wm_csf_df], axis=1
                    )
                else:
                    task_confounds_df = pd.concat(
                        [other_task_confounds_dfs[0], task_wm_csf_df], axis=1
                    )
                task_confounds_dfs.append(task_confounds_df)
            else:
                task_confounds_dfs.append(task_wm_csf_df)
    if for_rest:
        rest_wm_csf_df = pd.DataFrame({
            'wm_signal': rest_wm_mean,
            'csf_signal': rest_csf_mean
        })
        if other_rest_confounds_df is not None:
            rest_confounds_df = pd.concat(
                [other_rest_confounds_df, rest_wm_csf_df], axis=1
            )
        else:
            rest_confounds_df = rest_wm_csf_df
        
    if for_task and not for_rest:
        return task_confounds_dfs
    if not for_task and for_rest:
        return rest_confounds_df
    if for_task and for_rest:
        return task_confounds_df, rest_confounds_df      

def fit_fl_glm(runs, events_tables, confounds_tables):
    # Define First-level GLM
    fl_glm = FirstLevelModel(
        t_r=TASK_TR,
        slice_time_ref=SLICE_TIME_REF,
        smoothing_fwhm=SMOOTHING,
        hrf_model=HRF_MODEL,
        drift_model=DRIFT_MODEL,
        high_pass=HIGH_PASS,
        standardize=STANDARDIZE,
        n_jobs=-1,
    )

    # fit fl glm
    fl_glm.fit(runs, 
               events=events_tables,
               confounds=confounds_tables
              )
    
    return fl_glm

def fl_contrast_testing(glm, contrast_name, output_type='all'):
    contrast_outputs = glm.compute_contrast(
        contrast_def=CONTRASTS[contrast_name],
        stat_type='t',
        output_type=output_type
    )
    
    return contrast_outputs

def plot_map(map_type, input_map, map_title, mean_fmri, threshold=None, coords=[0,0,0], figsize=(8,4)):
    if map_type=='stat':
        plotting_config = {
            'bg_img': mean_fmri,
            'black_bg': True,
            'cut_coords':coords
        }
        if threshold:
            plot_stat_map(
                input_map,
                threshold=threshold,
                title=map_title,
                figure=plt.figure(figsize=figsize),
                **plotting_config
            )
        else:
            plot_stat_map(
                input_map,
                title=map_title,
                figure=plt.figure(figsize=figsize),
                **plotting_config
            )
    elif map_type=='roi':
        plotting_config = {
            'bg_img': mean_fmri,
            'black_bg': False,
            'cut_coords':coords
        }
        plot_roi(
            input_map,
            title=map_title,
            figure=plt.figure(figsize=figsize),
            **plotting_config
        )
    show()

def get_fl_contrast_map(sub_info, contrast_name, verbose=0):
    # get subject task scans & tables
    run1, run2, mean_fmri, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, wm_probmap, csf_probmap = get_maps_tables(sub_info, get_task=True, get_task_mean=True, get_rest=False, get_probmaps=True, get_anat=False)

    # get complete confounds tables
    task_confounds_tables = get_confounds(
        wm_probmap, csf_probmap, mean_fmri, 
        for_task=True, for_rest=False, 
        other_task_confounds_dfs=[run1_confounds_table, run2_confounds_table],
        other_rest_confounds_df=None, 
        task_fmri_runs=[run1, run2], rest_fmri=None, 
        threshold=WM_CSF_PROB_THRESH
    )
    
    if verbose==1:
        print("Fetched task-fMRI scans and events & confounds tables")

    # perform first-level GLM fitting
    fl_glm = fit_fl_glm(
        [run1, run2], 
        [run1_events_table, run2_events_table], 
        task_confounds_tables
    )
    if verbose==1:
        print("Fitted First-level GLM")

    # perform contrast computation
    contrast_map = fl_contrast_testing(fl_glm, contrast_name, 'effect_size')
    if verbose==1:
        print(f"Obtained contrast map for testing {contrast_name}")

    return contrast_map
    
def create_sl_design_matrix(sub_list, phenotype_table, show_matrix=False, save_matrix=False, save_path='..'):
    # filter phenotype table for given subjects
    filtered_pheno = phenotype_table[phenotype_table['participant_id'].isin(
        ['sub-0'+str(sub) for sub in sub_list]
    )]
    
    # create design matrix
    sl_design_matrix = pd.DataFrame()
    sl_design_matrix['age_centered'] = filtered_pheno['age_baseline'] - filtered_pheno['age_baseline'].mean()
    sl_design_matrix['IQ_centered'] = filtered_pheno['WASI_FULL_2_IQ'] - filtered_pheno['WASI_FULL_2_IQ'].mean()
    sl_design_matrix['sex'] = filtered_pheno['sex']
    sl_design_matrix['race_white'] = filtered_pheno['RACE_WHITE']
    sl_design_matrix['race_black'] = filtered_pheno['RACE_BLACK']
    sl_design_matrix['race_asian'] = filtered_pheno['RACE_ASIAN']
    sl_design_matrix['race_multiple'] = filtered_pheno['RACE_MULTIPLE']
    sl_design_matrix['hispanic'] = filtered_pheno['ETHNICITY']
    sl_design_matrix['cohort'] = filtered_pheno['COHORT']
    sl_design_matrix['scanner'] = filtered_pheno['SCANNER']
    for x in range(1,9):
        conditions = [
            filtered_pheno['INCOME'] == x,
            filtered_pheno['INCOME'] == 9
        ]
        choices = [1, -1]
        sl_design_matrix[f'income_level_{x}'] = np.select(conditions, choices, default=0)
    sl_design_matrix['missing_IQ'] = filtered_pheno['IQ_IS_MISSING']
    sl_design_matrix['missing_income'] = filtered_pheno['INCOME_IS_MISSING']
    sl_design_matrix['intercept'] = np.ones(len(sub_list))

    # show design matrix
    if show_matrix:
        plot_design_matrix(sl_design_matrix)
        show()

    # save design matrix
    if save_matrix:
        save_output_file(
            file=sl_design_matrix, 
            path=save_path, 
            filename="group_sl_design_matrix.csv"
        )
        
    return sl_design_matrix

def get_group_mean_task_fmri(subs_info, save_img=False, save_path='..'):
    # fetch mean tfmris of all subjects
    mean_fmris = []
    for sub_info in subs_info:
        mean_fmri = get_maps_tables(
            sub_info, get_task=False, get_task_mean=True, 
            get_rest=False, get_probmaps=False, get_anat=False
        )
        mean_fmris.append(mean_fmri)

    # obtain mean of list of mean tfmris
    group_mean_fmri = mean_img(mean_fmris)

    # save group mean tfmri
    if save_img:
        save_output_file(
            file=group_mean_fmri, 
            path=save_path, 
            filename="group_mean_tfmri.nii"
        )
        
    return group_mean_fmri

def get_group_mean_task_zmap(group_mean_fmri, subs_contrast_map, sl_design_matrix, sl_verbose=0, save_img=False, save_path='..'):
    # get group-mean brain-mask
    brain_mask = compute_epi_mask(group_mean_fmri)
    # plot_map('roi', brain_mask, 
    #              "Brain Mask", 
    #              group_mean_fmri, coords=[0,0,0], figsize=(8,9))

    # fit sl glm
    sl_glm = SecondLevelModel(
        mask_img=brain_mask,
        n_jobs=-1,
        verbose=sl_verbose
    )
    sl_glm.fit(
        second_level_input=subs_contrast_map,
        design_matrix=sl_design_matrix
    )
    
    # test intercept for sl glm
    mean_task_zmap = sl_glm.compute_contrast(
        'intercept', output_type='z_score'
    )

    # save zmap
    if save_img:
        save_output_file(
            file=mean_task_zmap, 
            path=save_path, 
            filename="group_mean_task_zmap.nii"
        )

    return mean_task_zmap

def get_group_mean_maps(sub_list, phenotype_table, fl_contrast_name='incongruent_vs_congruent', fl_contrast_map_verbose=0, sl_verbose=0, show_sl_design_matrix=False, save_group_mean_tfmri=False, save_group_mean_tfmri_path=HOLDOUTS_PATH, save_fl_contrasts=False, save_fl_contrasts_path=HOLDOUTS_PATH, save_sl_design_matrix=False, save_sl_design_matrix_path=HOLDOUTS_PATH, save_mean_task_zmap=False, save_mean_task_zmap_path=HOLDOUTS_PATH):

    # create list of small dictionaries with important subject attributes
    subs_info = get_sub_info_list(sub_list, phenotype_table)

    # get group's mean-task-fmri
    group_mean_tfmri = get_group_mean_task_fmri(
        subs_info, save_group_mean_tfmri, save_group_mean_tfmri_path
    )
    print("Computed group mean task-fmri")

    # get subject fl-glm contrast maps
    try:
        subs_contrast_map = list(iter_img(load_file(
            path=save_fl_contrasts_path, 
            filename=f"fl_{fl_contrast_name}-contrast_maps_{len(sub_list)}subs.nii.gz"
        )))
    except:
        subs_contrast_map = []
        if fl_contrast_map_verbose==2:
            print("Subjects completed:")
        for sub_info in subs_info:
            if fl_contrast_map_verbose==1:
                print(f"\nSubject {sub_info['subject_id']} - {sub_info['diagnosis']} - Cohort {sub_info['cohort']}")
            sub_contrast_map = get_fl_contrast_map(sub_info, fl_contrast_name, fl_contrast_map_verbose)
            subs_contrast_map.append(sub_contrast_map)
            if fl_contrast_map_verbose==2:
                print(f"{sub_info['subject_id']}", end=", ")
        # save 4d contrast maps
        if save_fl_contrasts:
            save_output_file(
                file=concat_imgs(subs_contrast_map), 
                path=save_fl_contrasts_path, 
                filename=f"fl_{fl_contrast_name}-contrast_maps_{len(sub_list)}subs.nii.gz"
            )
        print()
    print("Obtained first-level contrast maps")

    # get sl design matrix
    sl_design_matrix = create_sl_design_matrix(
        sub_list, phenotype_table, show_sl_design_matrix, 
        save_sl_design_matrix, save_sl_design_matrix_path
    )
    print("Created second-level design matrix")
    
    # get group's mean-task-zmap
    group_mean_task_zmap = get_group_mean_task_zmap(
        group_mean_tfmri, subs_contrast_map, sl_design_matrix, 
        sl_verbose, save_mean_task_zmap, save_mean_task_zmap_path
    )
    print("Obtained Group mean task maps")
    
    return group_mean_tfmri, group_mean_task_zmap

def stat_map_thresholding(stat_map, threshold=3.1, alpha=None, height_control=None, two_sided=False, cluster_threshold=0, save_img=False, save_name='thresholded', save_path='..'):
    sided = "two-sided" if two_sided else "one-sided"
    if alpha:
        if cluster_threshold==0:
            clean_map, clean_threshold = threshold_stats_img(
                stat_map,
                alpha=alpha,
                height_control=height_control,
                two_sided=two_sided
            )
        else:
            clean_map, clean_threshold = threshold_stats_img(
                stat_map,
                alpha=alpha,
                height_control=height_control,
                cluster_threshold=cluster_threshold,
                two_sided=two_sided
            )
        print(f"Threshold: {clean_threshold:.3f}")
        corrected = "uncorrected" if height_control=='fpr' else f"{height_control}-corrected"
        if cluster_threshold==0:
            sub_directory = f"{corrected}_alpha-{alpha}_{sided}"
        else:
            sub_directory = f"{corrected}_alpha-{alpha}_{sided}_cluster-thresh-{cluster_threshold}"
        save_path = os.path.join(save_path, sub_directory)
        if save_img:
            filename = f"{save_name}_map.nii"
            json_file = {
                'clean_map': save_output_file(clean_map, save_path, filename),
                'threshold': clean_threshold
            }
            json_filename = f"{save_name}_info.json"
            save_output_file(json_file, save_path, json_filename)
    else:
        if cluster_threshold==0:
            clean_map, clean_threshold = threshold_stats_img(
                stat_map,
                threshold=threshold,
                height_control=None,
                two_sided=two_sided
            )
        else:
            clean_map, clean_threshold = threshold_stats_img(
                stat_map,
                threshold=threshold,
                height_control=None,
                cluster_threshold=cluster_threshold,
                two_sided=two_sided
            )
        if cluster_threshold==0:
            sub_directory = f"threshold-{threshold}_{sided}"
        else:
            sub_directory = f"threshold-{threshold}_{sided}_cluster-thresh-{cluster_threshold}"
        save_path = os.path.join(save_path, sub_directory)
        if save_img:
            filename = f"{save_name}_map.nii"
            save_output_file(clean_map, save_path, filename)
            
    return clean_map, clean_threshold, sub_directory

def get_clusters(clean_z_map, threshold_z, cluster_threshold=0, two_sided=False):
    cluster_table, cluster_label_map = get_clusters_table(
        stat_img=clean_z_map,
        stat_threshold=threshold_z,
        cluster_threshold=cluster_threshold,
        two_sided=two_sided,
        return_label_maps=True
    )
    cluster_map = cluster_label_map[0]
    cluster_coordinates = find_parcellation_cut_coords(cluster_map)
    
    return cluster_table, cluster_map, cluster_coordinates

def get_resampled_atlas(mean_fmri, atlas_name='schaefer', n_rois=100, yeo_networks=7, resolution_mm=PREP_RESOLUTION):
    # check if exists
    if atlas_name == 'schaefer':
        atlas_info = f"{n_rois}_{yeo_networks}"
    atlas_filename = f"resampled_{atlas_name}_{atlas_info}.nii"
    atlas_table_filename = f"resampled_{atlas_name}_{atlas_info}_label_coords.csv"
    try:
        resampled_atlas = load_file(RESAMPLED_ATLAS_PATH, atlas_filename)
        atlas_table = load_file(RESAMPLED_ATLAS_PATH, atlas_table_filename)
        atlas_table['coordinates'] = atlas_table['coordinates'].apply(safe_parse)
    except:
        # fetch atlas
        if atlas_name == 'schaefer':
            atlas = datasets.fetch_atlas_schaefer_2018(
                n_rois=n_rois,
                yeo_networks=yeo_networks,
                resolution_mm=resolution_mm,
                verbose=0
            )
            atlas_labels = [l[11:] if l.startswith(f'{yeo_networks}Networks') else l for l in atlas.labels]
            
        # resample atlas
        resampled_atlas = resample_to_img(atlas.maps, mean_fmri, interpolation='nearest')
    
        # get atlas ROI coordinates
        atlas_coordinates = [list(c) for c in find_parcellation_cut_coords(resampled_atlas)]

        # create dataframe for labels & coordinates
        atlas_table = pd.DataFrame(data={
            'label': atlas_labels[1:],
            'coordinates': atlas_coordinates
        })
        # save resampled atlas & coordinates
        save_output_file(resampled_atlas, RESAMPLED_ATLAS_PATH, atlas_filename)
        save_output_file(atlas_table, RESAMPLED_ATLAS_PATH, atlas_table_filename)
        
    return resampled_atlas, atlas_table, atlas_info

def get_masked_z_scores(clean_z_map, mask, return_masker=False):
    task_masker = NiftiLabelsMasker(
        labels_img=mask,
        standardize=STANDARDIZE
    )
    masked_z_scores = task_masker.fit_transform(clean_z_map).flatten()
    
    if return_masker:
        return masked_z_scores, task_masker
        
    return masked_z_scores

def get_active_rois_mask(masked_z_scores, resampled_atlas, atlas_table, atlas_name='schaefer', atlas_info='', save_img=False, save_path='..'):
    # Identify indices of ROIs where z-score!=0
    active_rois_indices = np.where(masked_z_scores!=0)[0]
    
    # Convert ROI indices to Atlas indices
    # atlas labels start at 1 (0 is bg) : masker output index 0 => label 1
    active_atlas_indices = active_rois_indices + 1

    # Create spatial mask
    atlas_data = resampled_atlas.get_fdata()
    active_mask_data = np.where(np.isin(atlas_data, active_atlas_indices), atlas_data, 0)
    active_mask_img = new_img_like(resampled_atlas, active_mask_data)

    # Get active atlas coordinates
    active_coordinates = atlas_table['coordinates'].values[active_rois_indices]

    # Get active atlas labels
    active_labels = atlas_table['label'].values[active_rois_indices]

    # Create dataframe for labels & coordinates
    active_atlas_table = pd.DataFrame(data={
        'label': active_labels,
        'coordinates': active_coordinates
    })

    # save files
    if save_img:
        filename = f"{len(active_coordinates)}ROIs_{atlas_name}{atlas_info}"
        save_output_file(active_mask_img, save_path, f"Mask-{filename}.nii")
        save_output_file(active_atlas_table, save_path, f"Table-{filename}.csv")
    
    return active_mask_img, active_atlas_table

def get_group_task_ROIs(group_mean_zmap, group_mean_fmri, alpha=ALPHA, height_control=HEIGHT_CONTROL, two_sided=False, cluster_threshold=0, atlas_name='schaefer', n_rois=100, yeo_networks=7, save_thresholded_zmap=False, save_thresholded_zmap_path=HOLDOUTS_PATH, save_ROIs_mask=False, save_ROIs_mask_path=HOLDOUTS_PATH):
    # thresholding of group-mean-zmap
    clean_zmap, zmap_threshold, sub_directory = stat_map_thresholding(
        group_mean_zmap, 
        alpha=alpha, 
        height_control=height_control, 
        two_sided=two_sided, 
        cluster_threshold=cluster_threshold,
        save_img=save_thresholded_zmap, 
        save_name="thresholded_group_mean_task",
        save_path=save_thresholded_zmap_path
    )
    save_ROIs_mask_path = os.path.join(save_ROIs_mask_path, sub_directory)
    print("Thresholded the group-task-mean Z-map")
    
    # resampling atlas
    resampled_atlas, atlas_table, atlas_info = get_resampled_atlas(
        group_mean_fmri, atlas_name=atlas_name, n_rois=n_rois, yeo_networks=yeo_networks
    )
    print("Fetched resampled atlas map, labels, and coordinates")

    # getting atlas-masked group-mean-zscores
    masked_mean_zscores = get_masked_z_scores(clean_zmap, resampled_atlas)
    print("Applied atlas on Z-map")

    # getting group-task ROIs mask
    group_ROIs_mask, group_ROIs_table = get_active_rois_mask(
        masked_mean_zscores,
        resampled_atlas, 
        atlas_table,
        atlas_name=atlas_name,
        atlas_info=atlas_info,
        save_img=save_ROIs_mask, 
        save_path=save_ROIs_mask_path
    )
    print("Obtained group ROIs mask, labels, and coordinates")
    
    return group_ROIs_mask, group_ROIs_table, sub_directory

def get_masked_t_fc(runs, runs_confounds, mask):
    task_masker = NiftiLabelsMasker(
        labels_img=mask,
        t_r=TASK_TR,
        detrend=True,
        high_pass=HIGH_PASS,
        low_pass=LOW_PASS,
        standardize=STANDARDIZE
    )
    masked_task_ts = task_masker.fit_transform(
        runs,
        confounds=runs_confounds
    )
    # Partial correlation
    partial_corr_measure = ConnectivityMeasure(
        kind='partial correlation', standardize=STANDARDIZE
    )
    partial_corr_matrix = partial_corr_measure.fit_transform(masked_task_ts)[0]

    return partial_corr_matrix

def get_masked_rs_fc(rest, rest_confounds, mask):
    rest_masker = NiftiLabelsMasker(
        labels_img=mask,
        t_r=REST_TR,
        detrend=True,
        high_pass=HIGH_PASS,
        low_pass=LOW_PASS,
        standardize=STANDARDIZE
    )
    masked_rest_ts = rest_masker.fit_transform(
        rest,
        confounds=rest_confounds
    )
    # Partial correlation
    partial_corr_measure = ConnectivityMeasure(
        kind='partial correlation', standardize=STANDARDIZE
    )
    partial_corr_matrix = partial_corr_measure.fit_transform(masked_rest_ts)[0]

    return partial_corr_matrix
    
def get_task_rest_features(sub_info, group_ROIs_mask, group_ROIs_coords, fl_contrast_name, fl_alt_test=None, task_zmap=None, ROIs_task_zmap=None, ROIs_rest_FC_matrix=None, show_zmap=False, show_ROIs_zmap=False, show_ROIs_connectome=False, coords=[0,0,0], figsize=(8,4)):
    # get subject scans
    run1, run2, mean_fmri, rest, run1_events_table, run2_events_table, run1_confounds_table, run2_confounds_table, rest_confounds_table, wm_probmap, csf_probmap = get_maps_tables(sub_info, get_task=True, get_task_mean=True, get_rest=True, get_probmaps=True, get_anat=False)

    # get complete confounds tables
    task_confounds_tables, rest_confounds_table_final = get_confounds(
        wm_probmap, csf_probmap, mean_fmri, 
        for_task=True, for_rest=True, 
        other_task_confounds_dfs=[run1_confounds_table, run2_confounds_table],
        other_rest_confounds_df=rest_confounds_table, 
        task_fmri_runs=[run1, run2], rest_fmri=rest, 
        threshold=WM_CSF_PROB_THRESH
    )

    # Task-features
    fl_glm = fit_fl_glm(
        [run1, run2], 
        [run1_events_table, run2_events_table], 
        task_confounds_tables
    )
    if task_zmap is None:
        # print("computing task zmap")
        task_zmap = fl_contrast_testing(fl_glm, fl_contrast_name, 'z_score')
    if show_zmap:
        plot_map('stat', task_zmap, 
                 f"Subject {sub_info['subject_id']}: {fl_alt_test}, Z-score Map", 
                 mean_fmri, coords=coords, figsize=figsize
        )
    ROIs_task_zscores, fitted_task_masker = get_masked_z_scores(
        task_zmap, group_ROIs_mask, return_masker=True
    )
    if ROIs_task_zmap is None:
        # print("computing ROIs task zmap")
        ROIs_task_zmap = fitted_task_masker.inverse_transform(ROIs_task_zscores)
    if show_ROIs_zmap:
        plot_map('stat', ROIs_task_zmap, 
                 f"Subject {sub_info['subject_id']}: {fl_alt_test}, ROIs Z-score Map", 
                 mean_fmri, coords=coords, figsize=figsize
        )

    # ROIs_task_FC_matrix1 = get_masked_t_fc(run1, run1_confounds_table, group_ROIs_mask)
    # plot_connectome(
    #     ROIs_task_FC_matrix1,
    #     group_ROIs_coords,
    #     figure=plt.figure(figsize=figsize),
    #     edge_vmin=-1.0,
    #     edge_vmax=1.0,
    #     colorbar=True,
    #     title=f"Subject {sub_info['subject_id']} ROIs Task-Run1 Functional Connectivity Map"
    # )
    # show()
    # ROIs_task_FC_matrix2 = get_masked_t_fc(run2, run2_confounds_table, group_ROIs_mask)
    # plot_connectome(
    #     ROIs_task_FC_matrix2,
    #     group_ROIs_coords,
    #     figure=plt.figure(figsize=figsize),
    #     edge_vmin=-1.0,
    #     edge_vmax=1.0,
    #     colorbar=True,
    #     title=f"Subject {sub_info['subject_id']} ROIs Task-Run2 Functional Connectivity Map"
    # )
    # show()
    
    # Rest-features
    if ROIs_rest_FC_matrix is None:
        # print("computing ROIs rest FC matrix")
        ROIs_rest_FC_matrix = get_masked_rs_fc(
            rest, rest_confounds_table_final, group_ROIs_mask
        )
    if show_ROIs_connectome:
        plot_connectome(
            ROIs_rest_FC_matrix,
            group_ROIs_coords,
            figure=plt.figure(figsize=figsize),
            edge_vmin=-1.0,
            edge_vmax=1.0,
            colorbar=True,
            title=f"Subject {sub_info['subject_id']} ROIs Resting-state Functional Connectivity Map"
        )
        show()
    upper_tri_indices = np.triu_indices(len(group_ROIs_coords), k=1) #exclude diagonal
    vectorized_ROIs_rest_fc = ROIs_rest_FC_matrix[upper_tri_indices]
    ROIs_rest_zscores = np.arctanh(vectorized_ROIs_rest_fc)
    
    return task_zmap, ROIs_task_zmap, ROIs_rest_FC_matrix, ROIs_task_zscores, ROIs_rest_zscores

def combine_task_rest_features(task_zscores, rest_zscores):
    print(f"Min. task z-score: {min(task_zscores):.3f}, Max. task z-score: {max(task_zscores):.3f}")
    print(f"Length of task features: {len(task_zscores)}")
    print(f"Min. rest z-score: {min(rest_zscores):.3f}, Max. rest z-score: {max(rest_zscores):.3f}")
    print(f"Length of rest features: {len(rest_zscores)}")
    feature_vector = np.concatenate([task_zscores, rest_zscores])
    print(f"Length of feature vector: {len(feature_vector)}")
    return feature_vector

def create_ml_dataset(sub_list, phenotype_table, group_ROIs_mask, group_ROIs_table, fl_contrast_name, fl_alt_test, show_zmaps=False, show_ROIs_zmaps=False, show_ROIs_connectomes=False, coords=[0,0,0], figsize=(8,4), feature_extraction_verbose=2, save_zmaps=False, save_zmaps_path=HOLDOUTS_PATH, save_ROIs_zmaps=False, save_ROIs_zmaps_path=HOLDOUTS_PATH, save_ROIs_rest_FC_matrices=False, save_ROIs_rest_FC_matrices_path=HOLDOUTS_PATH, save_table=False, save_name="ml_dataset.csv", save_path=HOLDOUTS_PATH):
    
    # create list of small dictionaries with important subject attributes
    subs_info = get_sub_info_list(sub_list, phenotype_table)

    # create feature table
    task_feature_names = list(group_ROIs_table['label'].values)
    rest_feature_names = [f"{a}_to_{b}" 
                          for a,b in itertools.combinations(
                              task_feature_names, 2
                          )]
    ml_dataset = pd.DataFrame(
        columns=['participant_id']+task_feature_names+rest_feature_names
    )
    
    # Get Task & Rest features
    group_ROIs_table['coordinates'] = group_ROIs_table['coordinates'].apply(safe_parse)
    group_ROIs_coords = np.array([np.array(c) for c in group_ROIs_table['coordinates'].values.tolist()])
    
    # pre-output flags
    zmaps_filename = f"fl_{fl_contrast_name}-contrast_zmaps_{len(sub_list)}subs.nii.gz"
    ROIs_zmaps_filename = f"fl_{fl_contrast_name}-contrast_zmaps_{len(group_ROIs_table)}ROIs_{len(sub_list)}subs.nii.gz"
    ROIs_rest_FC_matrices_filename = f"rs-FunctionalConnectivity-matrices_{len(group_ROIs_table)}ROIs_{len(sub_list)}subs.npy"
    compute_zmap, compute_ROIs_zmap, compute_ROIs_rest_FC_matrix = True, True, True
    try:
        task_zmaps = list(iter_img(
            load_file(save_zmaps_path, zmaps_filename)
        ))
        compute_zmap = False
    except:
        task_zmaps = []

    try:
        ROIs_task_zmaps = list(iter_img(
            load_file(save_ROIs_zmaps_path, ROIs_zmaps_filename)
        ))
        compute_ROIs_zmap = False
    except:
        ROIs_task_zmaps = []

    try:
        ROIs_rest_FC_matrices = load_file(
            save_ROIs_rest_FC_matrices_path, ROIs_rest_FC_matrices_filename
        )
        compute_ROIs_rest_FC_matrix = False
    except:
        ROIs_rest_FC_matrices = []

    if feature_extraction_verbose==2:
        print("Subjects completed:")
    
    for i,sub_info in enumerate(subs_info):
        try:
            # print(len(task_zmaps), type(task_zmaps))
            task_zmap = task_zmaps[i]
            # plot_stat_map(task_zmap)
            # show()
        except:
            task_zmap = None

        try:
            # print(len(ROIs_task_zmaps), type(ROIs_task_zmaps))
            ROIs_task_zmap = ROIs_task_zmaps[i]
            # plot_stat_map(ROIs_task_zmap)
            # show()
        except:
            ROIs_task_zmap = None

        try:
            ROIs_rest_FC_matrix = ROIs_rest_FC_matrices[i]
        except:
            ROIs_rest_FC_matrix = None
            
        task_zmap, ROIs_task_zmap, ROIs_rest_FC_matrix, task_features, rest_features = get_task_rest_features(
            sub_info, 
            group_ROIs_mask, 
            group_ROIs_coords, 
            fl_contrast_name, 
            fl_alt_test=fl_alt_test, 
            task_zmap=task_zmap,
            ROIs_task_zmap=ROIs_task_zmap,
            ROIs_rest_FC_matrix=ROIs_rest_FC_matrix,
            show_zmap=show_zmaps, 
            show_ROIs_zmap=show_ROIs_zmaps, 
            show_ROIs_connectome=show_ROIs_connectomes, 
            coords=coords, figsize=figsize
        )
        
        ml_dataset.loc[len(ml_dataset)] = np.concatenate(
            [['sub-0'+str(sub_info['subject_id'])], task_features, rest_features]
        )
        if compute_zmap:
            task_zmaps.append(task_zmap)
        if compute_ROIs_zmap:
            ROIs_task_zmaps.append(ROIs_task_zmap)
        if compute_ROIs_rest_FC_matrix:
            ROIs_rest_FC_matrices.append(ROIs_rest_FC_matrix)

        if feature_extraction_verbose==2:
            print(f"{sub_info['subject_id']}", end=", ")
    print("\nCreated ML Features table") 

    # save results
    if save_table:
        save_output_file(ml_dataset, save_path, save_name)
    
    if compute_zmap and save_zmaps:
        save_output_file(concat_imgs(task_zmaps),
                         save_zmaps_path, 
                         zmaps_filename
                        )
        
    if compute_ROIs_zmap and save_ROIs_zmaps:
        save_output_file(concat_imgs(ROIs_task_zmaps),
                         save_ROIs_zmaps_path, 
                         ROIs_zmaps_filename
                        )
        
    if compute_ROIs_rest_FC_matrix and save_ROIs_rest_FC_matrices:
        save_output_file(np.array(ROIs_rest_FC_matrices),
                         save_ROIs_rest_FC_matrices_path, 
                         ROIs_rest_FC_matrices_filename
                        )

    return ml_dataset