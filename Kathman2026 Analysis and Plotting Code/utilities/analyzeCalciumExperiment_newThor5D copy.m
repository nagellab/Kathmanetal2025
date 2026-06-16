function analyzeCalciumExperiment_newThor5D(headpath,fileexpnum,correction,fluor,new_rois)
%give this function the experimental path you want to start in
%and the filename and experiment number (ie vt019532_f3e4)
%correction if the motion correction needs to be done, fluor if the dff
%hasn't been extracted yet.

eeees = strfind(fileexpnum, 'e');
us = strfind(fileexpnum, '_');
if us(end) > eeees(end)
    expt = str2double(fileexpnum(eeees(end)+1:us(end)-1)); %pybmt expt #
else
    expt = str2double(fileexpnum(eeees(end)+1:end)); %pybmt expt #
end

if correction
    trlcat = 1; zsplit = 0;
    [trl_dirs_corr, FOV] = NDK_Motion_Correction5D(fileexpnum, trlcat, zsplit);

    save(['E' num2str(expt) '_regroi_vars.mat'], 'trl_dirs_corr', 'FOV')
    %     load regroi_test_vars.mat
%     FOV = [96, 160];
%     Correct_all5D(headpath,fileexpnum); %motion correct all data
    disp('Confirm ROI set is saved to first trial folder. Hit any key when ready...')
    pause

    base = trl_dirs_corr{1};
    other = trl_dirs_corr(2:end);
    other = other';
    register_rois2(base,other, FOV);

elseif new_rois
    load(['E' num2str(expt) '_regroi_vars.mat'])
    disp('Confirm ***NEW*** ROI set is saved to first trial folder. Hit any key when ready...')
    pause
    base = trl_dirs_corr{1};
    other = trl_dirs_corr(2:end);
    other = other';
    register_rois2(base,other, FOV);
end


%analyze the data for each
[Xpix,Ypix,FrameRate,Zframes,flybackframes]=findkeyparams_newThor5D(fileexpnum);
FrameRate=str2double(FrameRate)/(Zframes+flybackframes);

if new_rois; 
    fluor = 1; end
if fluor
    copy=0; %only do this if can use same roi for every trial (will overwrite .zips if already in directory)
    if copy
        disp('Please select the master ROIset');
        copyROIs(headpath,fileexpnum);%copy the ROIdataset to all folders necessary
    end
    cd(headpath);
    filelist=dir;
    for k=1:numel(filelist)%should do something if it should skip a trial? why hasn't this been a problem before 
        dataname=filelist(k).name;
        if (contains(dataname,fileexpnum)&& contains(dataname,'corr'))
            %should pass the full path here I guess?
            dataname
            fullpath=strcat(headpath,'/');
            fullpath=strcat(fullpath,dataname);
            try %writes Cin (raw fluor) and Ain (roi?) * not actually dff (that done in sync2struct)
                %does no processing of image (except corr image?), just extracts
                dff2(fullpath,Ypix,Xpix);%runs CDNF_noparpoolsAM which imports ROIs, analyzes for each one, and saves extracted values
            catch
                disp(['well we had to skip ' fileexpnum]);
            end
            cd(headpath);%pop back up to the top so you can go down to the next trial
        end
        
    end
end

%% andrews code to plot pooled directional responses for each ROI (good to see quick summary for each ROI)

%CURRENTLY SET TO USE ROISET2 -> IF YOU WANT ORIGINALS WILL HAVE TO MAKE
%CHANGES IN dff2 and COPYROIS 
%SAVE IN SEPARATE EXTRACTED FOLDER!!! 
%MAKE PLOTS IN SEPARATE FOLDER NOT TO OVERWRITE 
%HOPEFULLY THIS HARDDRIVE HAS ENOUGH SPACE 
% % function analyzeCalciumExperiment_newThor5D(headpath,fileexpnum,imskip,expskip,correction,fluor,copy,savemode)

% plot_rois = 0;
% if plot_rois
% deltamode = 0;
% [fluordata,fluordata_red,fluordata_diff,ROIs,redtemps]=pooldirections(headpath,fileexpnum,imskip,expskip,deltamode,FrameRate);
% % need something here to pool airspeeds
% %ROIs=reshape(ROIs,[Ypix,Xpix]);%make it an easy 3-D matrix so each plane is a ROI image
% ROIs=full(ROIs);
% %meanBaseGreen=
% %avggreen=mean(CinGreen(:,1:5*framerate),2);
% %divGreen=CinGreen./avggreen
% 
% 
% %[dgG,drR,diffRG]=RGdiff(fluordata,fluordata_red)
% [avgfluor,stdfluor]=avgROIs(fluordata);
% [avgred,stdred]=avgROIs(fluordata_red);%get the red averages? -> probably want to convert to deltaF/F for each first
% % [avgdiff,stddiff]=avgROIs(fluordata_diff);
% % plotfluor_dirspd5D(avgfluor,redtemps,ROIs,Xpix,Ypix,FrameRate);
% 
% 
% us = strfind(fileexpnum, '_f');
% slsh = strfind(headpath, '/');
% outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum(us(end)+1:end) '_extracted'];
% 
% % us = strfind(fileexpnum, '_');
% % outdirname = [fileexpnum(us(end)+1:end) '_extracted'];
% if ~exist(outdirname,'dir')
%     mkdir(outdirname);
% end
% cd(outdirname)
% plotfluordirections_45only(avgfluor,redtemps,ROIs,Xpix,Ypix,FrameRate);
% cd(headpath)
% % %make a version of this that takes them all and then plots them together
% % if savemode
% %     us = strfind(fileexpnum, '_');
% %     outdirname = [fileexpnum(us(end)+1:end) '_extracted'];
% %     if ~exist(outdirname,'dir')
% %         mkdir(outdirname);
% %     end
% %     cd (outdirname);
% %     savefilename=strcat(fileexpnum,'_analysis2');
% %     save(savefilename,'fluordata','fluordata_red','fluordata_diff','ROIs','redtemps','avgfluor','stdfluor','avgred','stdred','avgdiff','stddiff','Xpix','Ypix','FrameRate');
% %     cd(headpath);
% % end
% 
% end



end