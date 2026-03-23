%-----------------------------------------------------------------------
% Job saved on 03-Feb-2026 14:07:41 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

subjects = [22172];
tau = '1';
for sub=subjects

    subject = num2str(sub);
    disp(['Unzipping MRI scans of Cohort 2 subject-' subject]);

    if isfile(['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-TAU' tau '_run-1_bold.nii.gz']) == 0
        tau = '2';
    end

    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = {
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\anat\sub-' subject '_ses-1_T1w.nii.gz']
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-TAU' tau '_run-1_bold.nii.gz']
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-TAU' tau '_run-2_bold.nii.gz']
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-rest_echo-1_bold.nii.gz']
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-rest_echo-2_bold.nii.gz']
                                                                           ['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' subject '\ses-1\func\sub-' subject '_ses-1_task-rest_echo-3_bold.nii.gz']
                                                                       };
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {''};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep = true;

    spm_jobman('run', matlabbatch);
end