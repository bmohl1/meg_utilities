function Oddball_batch_tft_bmm(list,coords,labels,ep)
%batching epoch extractor to split oddball files into conditions (10s
%separate from 20s)
addpath('path_to_your_MEG_tools_dir/tools/megcode_v2_Aug2013/megtools');
% Select working directories

if exist('list')
    pth_subjdirs = list;
else

    cwd             = spm_select(1,'dir','Select root directory for studies',...
        '',pwd);
    cd(cwd);
    pth_subjdirs    = spm_select([1,Inf],'dir','Select subject directories to process',...
        '',pwd);
end
nsub = size(pth_subjdirs,1);
fprintf('The following %d subject(s) will be preprocessed:\n',nsub);
disp(pth_subjdirs);


% loop through directories
for sub=1:nsub
    fprintf('working on subject %d of %d/n',sub,nsub);

    % change working directory
    cd(pth_subjdirs(sub,:));

    %get mat name from dat name to avoid other mat files
    if strcmp(pth_subjdirs(sub,end),filesep)
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end-1)); % get subject id from path
    else
        [pth id ext]     = fileparts(pth_subjdirs(sub,1:end));
    end

    if exist('ep','var')
        conds = {ep}
    else
        conds = {'ep10' 'ep20' 'correp20'};
    end

    for i = 1:length(conds)
        cond = conds{i};
        %Kristina's coordinates
        %mni_coords = [52 -10 4;
        %              -56 -12 4;
        %              40 18 20;
        %              -40 18 20;
        %              40 8 36;
        %              -40 8 36;
        %              56 -44 32;
        %              -56 -44 32;
        %              28 -62 60;
        %             -28 -62 60;

        %              54 -4 0;
        %                -48 -2 -14;
        %                46 -56 -8;
        %                54 -2 0;
        %                -54 2 10];

        %labels     = {'R_Aud',...
        %              'L_Aud',...
        %              'R_Finf',...
        %              'L_Finf'...
        %              'R_Fmid'...
        %              'L_Fmid'...
        %              'R_smg'...
        %              'L_smg'...
        %              'R_Spar'...
        %              'L_Spar'...
        %              'R_Audb',...
        %              'L_Audb',...
        %              'R_infTL',...
        %              'R_Aud20',...
        %              'L_Aud20'};

        %Brianne's coordinates
        mni_coords = coords;
        rois = labels{1,1};


        for j = 1:length(rois)
            roi=rois{j}
            % load file left aud ep20
            saveFile = [pth filesep id filesep id '_' cond '_' roi '_tft.mat'];
            if ~exist(saveFile,'file')
                file = ([id '_' roi '_' cond '_ssp.mat']);
                load(file);
                tf = tft(file,[5 80]);
                fprintf('time-frequency analysis...');

                %save tft file
                save ([id '_' cond '_' roi '_tft.mat'],'tf');
                clear tf;
            end

        end
    end
end
