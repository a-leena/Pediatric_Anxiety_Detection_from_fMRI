%-----------------------------------------------------------------------
% Job saved on 07-Feb-2026 18:29:28 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

subjects = [24106, 24115, 24118, 24139, 24146, 24148, 24151, 24182, 24217, 24236, 24240, 24260, 24300, 24353];

cohort = '2';

for subject=subjects    
    sub = num2str(subject);

    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.name = 'rest_file';
    if cohort=='1'
        disp(['Preprocessing Rest fMRI for cohort ' cohort ' subject-' sub]);
        matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\func\sub-' sub '_ses-1_task-rest_bold.nii']}};

        matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Named File Selector: rest_file(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
        matlabbatch{2}.spm.temporal.st.nslices = 36;
        matlabbatch{2}.spm.temporal.st.tr = 2;
        matlabbatch{2}.spm.temporal.st.ta = 1.94444444444444;
        matlabbatch{2}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36];
        matlabbatch{2}.spm.temporal.st.refslice = 1;
        matlabbatch{2}.spm.temporal.st.prefix = 'a';
        matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
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
        matlabbatch{4}.spm.spatial.coreg.estimate.ref = {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\anat\msub-' sub '_ses-1_T1w.nii,1']};
        matlabbatch{4}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
        matlabbatch{4}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
        matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
        matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
        matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
        matlabbatch{5}.spm.spatial.normalise.write.subj.def = {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\anat\y_sub-' sub '_ses-1_T1w.nii']};
        matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
        matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                  78 76 85];
        matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
        matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
        matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
    
        spm_jobman('run', matlabbatch);

    else
        echos = ['1', '2', '3'];
        for echo=echos
            disp(['Preprocessing Rest fMRI echo ' echo ' for cohort ' cohort ' subject-' sub]);
            matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_named_file.files = {{['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\func\sub-' sub '_ses-1_task-rest_echo-' echo '_bold.nii']}};
    
            matlabbatch{2}.spm.temporal.st.scans{1}(1) = cfg_dep('Named File Selector: rest_file(1) - Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files', '{}',{1}));
            matlabbatch{2}.spm.temporal.st.nslices = 34;
            matlabbatch{2}.spm.temporal.st.tr = 2;
            matlabbatch{2}.spm.temporal.st.ta = 1.94117647058824;
            matlabbatch{2}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34];
            matlabbatch{2}.spm.temporal.st.refslice = 1;
            matlabbatch{2}.spm.temporal.st.prefix = 'a';
            matlabbatch{3}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
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
            matlabbatch{4}.spm.spatial.coreg.estimate.ref = {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\anat\msub-' sub '_ses-1_T1w.nii,1']};
            matlabbatch{4}.spm.spatial.coreg.estimate.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
            matlabbatch{4}.spm.spatial.coreg.estimate.other(1) = cfg_dep('Realign: Estimate & Reslice: Resliced Images (Sess 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rfiles'));
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{4}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
            matlabbatch{5}.spm.spatial.normalise.write.subj.def = {['F:\MSc Capstone\Pediatric_Anxiety_Disorder\sub-' sub '\ses-1\anat\y_sub-' sub '_ses-1_T1w.nii']};
            matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
            matlabbatch{5}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                                      78 76 85];
            matlabbatch{5}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
            matlabbatch{5}.spm.spatial.normalise.write.woptions.interp = 4;
            matlabbatch{5}.spm.spatial.normalise.write.woptions.prefix = 'w';
        
            spm_jobman('run', matlabbatch);
        end
    end
end
