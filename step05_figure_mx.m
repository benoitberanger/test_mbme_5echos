clear
clc
close all

e = regroup_series();

vol = e.getSerie('run').getVolume('^s5wts').removeEmpty();
tags = {vol.getSerie.tag};

mdl_dir = gdir(e.getPath(),'mdl','clean','.*');
mx_file  = fullfile(mdl_dir,'correlation_matrix.mat');

load aal3.mat

nRun = length(mdl_dir);

nRegion = size(aal3,1);

f = figure('Name','Measured Data','NumberTitle','off');
tg = uitabgroup(f);

for iRun = 1 : nRun
    
    fprintf('run %d/%d : %s \n', iRun, nRun, mdl_dir{iRun})
        
    load(mx_file{iRun})
    
    t = uitab('Title', tags{iRun});
    a(iRun) = axes(t);
    
    % mx( eye(size(mx))==1 ) = NaN;
    
    imagesc(a(iRun), mx)
    caxis([-1 +1])
    colormap(jet)
    colorbar
    axis equal
    
    xticks(1:nRegion)
    xticklabels(aal3.ROIabbr)
    xtickangle(90)
    yticks(1:nRegion)
    yticklabels(aal3.ROIname)
    
    
end

linkaxes(a, 'xy')