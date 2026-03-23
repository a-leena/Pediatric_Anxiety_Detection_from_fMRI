%-----------------------------------------------------------------------
% Job saved on 07-Feb-2026 15:39:04 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

subjects = [24046, 24067, 24089, 24091, 24106, 24115, 24118, 24139, 24146, 24148, 24151, 24182, 24217, 24236, 24240, 24260, 24300, 24353];

cohort = '2';

if cohort=='1'
    nslices = 41;
    ta = 2.24390243902439;
    so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40];
else
    nslices = 42;
    ta = 2.2452380952381;
    so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42];
end

for subject=subjects
    sub = num2str(subject);

    disp(['Preprocessing Tasks fMRI for cohort ' cohort ' subject-' sub]);

    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'task_run1run2_files';
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {
                                                                         {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\func\sub-' sub '_ses-1_task-TAU' cohort '_run-1_bold.nii']}
                                                                         {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\func\sub-' sub '_ses-1_task-TAU' cohort '_run-2_bold.nii']}
                                                                         }';
    matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Named File Selector: task_run1run2_files(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
    matlabbatch{2}.spm.temporal.st.scans{2}(1) = cfg_dep('Named File Selector: task_run1run2_files(2) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{2}));
    matlabbatch{2}.spm.temporal.st.nslices = nslices;
    matlabbatch{2}.spm.temporal.st.tr = 2.3;
    matlabbatch{2}.spm.temporal.st.ta = ta;
    matlabbatch{2}.spm.temporal.st.so = so;
    matlabbatch{2}.spm.temporal.st.refslice = 1;
    matlabbatch{2}.spm.temporal.st.prefix = 'a';
    matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
    matlabbatch{3}.spm.spatial.realign.estwrite.data{2}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 2)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{2}, '.','files'));
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.quality = 0.95;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.sep = 1.5;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.fwhm = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{3}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    matlabbatch{4}.spm.spatial.preproc.channel.vols = {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\anat\sub-' sub '_ses-1_T1w.nii,1']};
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg = 0.0001;
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{4}.spm.spatial.preproc.channel.write = [0 1];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,1'};
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,2'};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,3'};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,4'};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,5'};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm = {'C:\Program Files\MATLAB\spm\tpm\TPM.nii,6'};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{4}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.reg = [0 0 0.1 0.01 0.04];
    matlabbatch{4}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{4}.spm.spatial.preproc.warp.samp = 3;
    matlabbatch{4}.spm.spatial.preproc.warp.write = [0 1];
    matlabbatch{4}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{4}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                                  NaN NaN NaN];
    matlabbatch{5}.spm.spatial.coreg.estimate.ref(1) = cfg_dep('Segment: INU corrected (1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','channel', '()',{1}, '.','biascorr', '()',{':'}));
    matlabbatch{5}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
    matlabbatch{5}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
    matlabbatch{5}.spm.spatial.coreg.estimate.other(2) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 2)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{2}, '.','rfiles'));
    matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{5}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    matlabbatch{6}.spm.spatial.normalise.write.subj.def(1) = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
    matlabbatch{6}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
    matlabbatch{6}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                              78 76 85];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{6}.spm.spatial.normalise.write.woptions.prefix = 'w';

    spm_jobman('run', matlabbatch);

end
