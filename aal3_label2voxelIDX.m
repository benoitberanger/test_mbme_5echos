clear
clc

load aal3.mat
nRegion = size(aal3,1);

Vaal = spm_vol('aal3.nii');
Yaal = spm_read_vols(Vaal);

for iRegion = 1 : nRegion
    ROIabbr = aal3.ROIabbr(iRegion);
    ROIid   = aal3.ROIid  (iRegion);
    
    mask_roi = Yaal == ROIid;
    
    aal3.voxel_list{iRegion} = find(mask_roi);
    
end

save aal3_voxel_list.mat aal3
