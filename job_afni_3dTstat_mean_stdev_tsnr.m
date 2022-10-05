function job = job_afni_3dTstat_mean_stdev_tsnr( img , par )
% JOB_AFNI_REMOVE_NAN - AFNI:3dcalc
%
% INPUT : img can be 'char' of volume(file), multi-level 'cellstr' of volume(file), '@volume' array
%
% To build the image list easily, use get_subdir_regex & get_subdir_regex_files
%
% See also get_subdir_regex get_subdir_regex_files exam exam.AddSerie exam.addVolume


%% Check input arguments

if ~exist('par','var')
    par = ''; % for defpar
end

if nargin < 1
    help(mfilename)
    error('[%s]: not enough input arguments - img is required',mfilename)
end

obj = 0;
if isa(img,'volume')
    obj = 1;
    volumeArray  = img;
    img = volumeArray.toJob(1);
end


%% defpar

defpar.OMP_NUM_THREADS = 0; % number pf CPU threads : 0 means all CPUs available

defpar.sge      = 0;
defpar.jobname  = 'job_afni_3dTstat_mean_stdev_tsnr';
defpar.walltime = '00:30:00';

defpar.auto_add_obj = 1;

defpar.pct      = 0;
defpar.redo     = 0;
defpar.run      = 1;
defpar.display  = 0;
defpar.verbose  = 1;

par = complet_struct(par,defpar);

if par.sge || par.pct
    par.OMP_NUM_THREADS = 1; % in case of parallelization, only use 1 thread per job
end


%% 3dTstat

nSubj = length(img);
job  = cell(0);

cjob = 0;
for iSubj = 1 : nSubj
    
    for vol = 1 : length(img{iSubj})
        
        src = deblank( img{iSubj}{vol} ) ;
        
        for fcn = {'mean', 'stdev', 'tsnr'}
            
            prefix = [    char(fcn) '_'];
            arg    = ['-' char(fcn)    ];
            
            % output exists ?
            dst = addprefixtofilenames( src , prefix );
            if ~par.redo   &&   exist(dst,'file')
                fprintf('[%s]: skiping subj/vol %d/%d because %s exist \n',mfilename,iSubj,vol,dst);
            else
                
                cjob = cjob + 1;
                
                cmd = sprintf('export OMP_NUM_THREADS=%d; 3dTstat -prefix %s %s %s; \n', par.OMP_NUM_THREADS, dst, arg, src);
                
                job{cjob,1} = char(cmd);
                
            end
            
        end
        
    end % vol
    
end % iSubj


%% Run the jobs

% Run CPU, run !
job = do_cmd_sge(job, par);


%% Add outputs objects

if obj && par.auto_add_obj && (par.run || par.sge)
    
    for iVol = 1 : length(volumeArray)
        
        % Shortcut
        vol = volumeArray(iVol);
        ser = vol.serie;
        tag = vol.tag;
        sub = vol.subdir;
        
        for fcn = {'mean', 'stdev', 'tsnr'}
            
            prefix = [ char(fcn) '_'];
            
            if par.run
                
                ext  = '.*.nii';
                
                ser.addVolume(sub, ['^' prefix tag ext],[prefix tag],1)
                
            elseif par.sge
                
                ser.addVolume('root', addprefixtofilenames(vol.path,prefix),[prefix tag])
                
            end
            
        end
        
    end % iVol
    
end % obj


end % function
