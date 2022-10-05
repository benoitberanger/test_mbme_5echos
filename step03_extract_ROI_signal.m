clear
clc

e = regroup_series();

vol = e.getSerie('run').getVolume('^s5wts').removeEmpty();
tags = {vol.getSerie.tag};

mdl_dir = gdir(e.getPath(),'mdl','clean','.*');
spm_file  = fullfile(mdl_dir,'SPM.mat');
mask_file = fullfile(mdl_dir,'mask.nii');

load aal3.mat

nRun = length(mdl_dir);

nRegion = size(aal3,1);

for iRun = 1 : nRun
    
    fprintf('run %d/%d : %s \n', iRun, nRun, mdl_dir{iRun})
    
    matlabbatch = cell(nRegion,1);
    for iRegion = 1 : nRegion
        
        ROIabbr = aal3.ROIabbr(iRegion);
        ROIid   = aal3.ROIid  (iRegion);
        
        matlabbatch{iRegion}.spm.util.voi.spmmat = spm_file(iRun);
        matlabbatch{iRegion}.spm.util.voi.adjust = NaN; % NaN means "adjust with everythong" => in RS, all regressors are noise, so we want to regress them out
        matlabbatch{iRegion}.spm.util.voi.session = 1;
        matlabbatch{iRegion}.spm.util.voi.name = sprintf('region_%03d_%s', ROIid, ROIabbr);
        matlabbatch{iRegion}.spm.util.voi.roi{1}.label.image = {'./aal3.nii'};
        matlabbatch{iRegion}.spm.util.voi.roi{1}.label.list = ROIid;
        matlabbatch{iRegion}.spm.util.voi.roi{2}.mask.image = {mask_file{iRun}};
        matlabbatch{iRegion}.spm.util.voi.roi{2}.mask.threshold = 0.5;
        matlabbatch{iRegion}.spm.util.voi.expression = 'i1 & i2';
        
    end
    
    spm_jobman('run',matlabbatch)
    
end
