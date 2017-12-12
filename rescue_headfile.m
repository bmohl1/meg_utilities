function rescue_headfile(list, ep)

disp('Rescuing headshape')
if exist('list')
    pth_subjdirs = list;
else
    cwd     = spm_select(1,'dir','Select root directory for studies',...
        '',pwd);
    cd(cwd);
    list    = spm_select([1,Inf],'dir','Select subject directories to process',...
        '',pwd);
end

for s = 1:size(list,1)
    subj = list(s, 1:end);
    subjName = textscan(subj, '%s', 'Delimiter', '/');
    subjName = char(subjName{1,1}{end});
    subjName = textscan(subjName, '%s', 'Delimiter', '_');
    subjName = char(subjName{1,1}{end});
    if exist('ep','var')
        options = {['_' ep '.mat']};
    else
        %options = {'_neweps.mat', '_spm.mat' '_cnt.mat'};
        options = {'_ep20.mat' '_cnt.mat' '_orig_cnt.mat' '_ep10.mat' '_correp20.mat' '_spm.mat'};
    end

    for i = 1: length(options)
        suffix = char(options{i});
        pth =[subj,subjName,suffix];
        dest=[subj,'fiducials.mat'];
        touchFile = [subj,'touch_hs.txt'];
        spmMat=rdir([subj, subjName,'_spm.mat']);
        if ~isempty(spmMat())
            spmMat = spmMat.name;
        else
            disp('Did not find *_spm.mat');
            continue
        end

        if exist(pth, 'file') %&& ~exist (touchFile, 'file')

            fprintf('Loading: %s\n',pth);

            if contains('20',{suffix})
                fName = 'ep20';
            elseif contains('10',{suffix});
                fName = 'ep10';
            elseif contains('orig',{suffix});
                fName = 'cnt_orig';
            elseif contains('cnt',{suffix});
                fName = 'cnt';
            elseif contains('spm',{suffix});
                fName = 'D';
            else
                fName = 'eps';
            end

            struct1 = load(pth);
            if ~isfield(struct1.(fName), 'fiducials')
                if ~exist('fiducials','var')
                    load(spmMat,'D')
                    fiducials = D.fiducials;
                    save(dest,'fiducials');
                end
                struct1.(fName).fiducials = fiducials;
                disp('>>>Saving fiducials<<<');
                save(pth,'-struct','struct1');
            else
                fprintf ('found fiducials - %s\n',fName);
            end
            clear cnt

        end

    end
    fclose(fopen(touchFile,'w+'));
end
