function[syncstruct]= sync2struct(headpath,fileexpnum, kill_column)
%[syncstruct]=
%sync2struct(headpath,fileexpnum,imskip,expskip,deltamode,savemode)

%analyze the data for each
[Xpix,Ypix,FrameRate,Zframes,flybackframes]=findkeyparams_newThor5D(fileexpnum);
FrameRate=str2double(FrameRate)/(Zframes+flybackframes);
trig = 3; %sec after pybmt starts to trig thor
% FrameRate = FrameRate - 1; 

cd(headpath);
filelist=dir;

eeees = strfind(fileexpnum, 'e');
us = strfind(fileexpnum, '_');
if us(end) > eeees(end)
    load(['E' fileexpnum(eeees(end)+1:us(end)-1) '.mat']); %pybmt expt #
else
    load(['E' fileexpnum(eeees(end)+1:end) '.mat']) %pybmt expt #
end

%% Andrews clean trials, but just returns sorted thor trial  names
%add some errorchecking for indices you want to remove or something
imskip = []; expskip = [];
[sortednames,data]=cleantrials(data,fileexpnum,filelist,imskip,expskip); %i think just returns list of thor trial names sorted by trial number
syncstruct = data;
for t = 1:size(syncstruct,2)
    balltime(t) = datenum(syncstruct(t).filename(end-9:end-4), 'HHMMSS');
end

%% find "all trial" baseline (avg of lowest 3% of values) for each ROI
all_CinGreen = [];
all_CinRed = [];
for k = 1:numel(sortednames)
    %find xml timestamp
    cd([headpath '/' sortednames{k}(1:end-5)])
    xmlexp=parseXML5D('Experiment.xml');
    kids=xmlexp.Children;
    trl_time{k} = kids(4).Attributes(1).Value;
    twoptime = datenum(trl_time{k}(end-7:end), 'HH:MM:SS');
    [m,synctrial] = min(abs(twoptime-balltime));
    cd(headpath)
    if m < 0.0004%approx a 40s
        cd(sortednames{k})
        if isfile('MC_Videocdf2.mat')
            load('MC_Videocdf2.mat');
            all_CinGreen=horzcat(all_CinGreen, Cin);
            load('MC_Video_Redcdf2.mat');%load in the red channel 
            all_CinRed=horzcat(all_CinRed, Cin);
        end
    end
end
sorted_allCinGreen = sort(all_CinGreen, 2);
green_bases = mean(sorted_allCinGreen(:,1:round(size(all_CinGreen,2)*0.03)),2);
sorted_allCinRed = sort(all_CinRed, 2);
red_bases = mean(sorted_allCinRed(:,1:round(size(all_CinRed,2)*0.03)),2);

%% load fluor data, find dff, add to struct
for k=1:numel(sortednames)
    cd([headpath '/' sortednames{k}(1:end-5)])
    xmlexp=parseXML5D('Experiment.xml');
    kids=xmlexp.Children;
    trl_time{k} = kids(4).Attributes(1).Value;
    twoptime = datenum(trl_time{k}(end-7:end), 'HH:MM:SS');
    [m,synctrial] = min(abs(twoptime-balltime));
    cd(headpath)
    if m < 0.0004 %approx a 40s
        cd(sortednames{k})
        if isfile('MC_Videocdf2.mat')
%         try
            %load fluor data and find df/f
            load('MC_Parameters.mat');
            load('MC_Videocdf2.mat');
            dfG=(Cin - green_bases)./green_bases;
            
            if ~isempty(kill_column)
                for kc = 1:length(kill_column) %replace "kill_columns" with nans
                    dfG(kill_column(kc), :) = nan(1,size(dfG,2));
                end
            end

            load('MC_Video_Redcdf2.mat');%load in the red channel 
            dfR=(Cin - red_bases)./red_bases;
%             dfR = dfR.*(trapz(dfG,2)./trapz(dfR,2)); %normalize red channel to green (no good bc change each roi seperates,losing point of red channel)
            m = mean((trapz(dfG,2)./trapz(dfR,2)));

            if ~isempty(kill_column)
                for kc = 1:length(kill_column) %replace "kill_columns" with nans
                    dfR(kill_column(kc), :) = nan(1,size(dfR,2));
                end
            end

            syncstruct(synctrial).('dff_green') = dfG;
            syncstruct(synctrial).('dff_red') = dfR;
            syncstruct(synctrial).('dff_diff') = dfG-dfR*m*.02; %from ratio of integrals    (.1; %scale found empirically, related to intensity diff of 2 fluors?)
            [numROIs,triallength]=size(Cin);
            %makes ts by framerate, not real ts, and add trigger time at end (this is how andrew did it, prob not best way, but have no thor ts)
            syncstruct(synctrial).('dff_ts') =(linspace(1,triallength,triallength)/(FrameRate)) + trig; 
            syncstruct(synctrial).('thorname') = sortednames{k};
            syncstruct(synctrial).('thortime') = [trl_time{k}];
%         catch
%             disp('fuck');
%             keyboard
%         end
        end
        cd(headpath);
    end
end

%% interpolate and add bump analysis
if ~isempty(syncstruct(1).dff_green) && size(syncstruct(1).dff_green,1)>15 %look for two layers of rois
    use_cols = [1:8; 9:16];
elseif ~isempty(syncstruct(2).dff_green) && size(syncstruct(2).dff_green,1)>15 %in case first trial empty
    use_cols = [1:8; 9:16];
elseif ~isempty(syncstruct(3).dff_green) && size(syncstruct(3).dff_green,1)>15 %in case second trial empty
    use_cols = [1:8; 9:16];
else
    use_cols = [1:8]; %2:7; %1:8
end
% use_cols = 1:8;
syncstruct = bumpstats(syncstruct, use_cols);

%% save stuct

%remove unused fields
fields = {'light_on', 'light', 'ai', 'gain', 'plume_dx', 'plume_dy', 'plume_odor', 'plume_codor', 'plume_A', 'plume_beta'};
existing_fields = intersect(fields, fieldnames(syncstruct));
syncstruct = rmfield(syncstruct, existing_fields);
% syncstruct = rmfield(syncstruct,fields);

savemode = 1;
if savemode
    us = strfind(fileexpnum, '_');
    slsh = strfind(headpath, '/');
    outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum(us(end)+1:end) '_extracted'];
    
    if ~exist(outdirname,'dir')
        mkdir(outdirname);
    end
    cd (outdirname);
%     savefilename=strcat(fileexpnum,'_analysis2');
%     save(savefilename,'fluordata','fluordata_red','fluordata_diff','ROIs','redtemps','avgfluor','stdfluor','avgred','stdred','avgdiff','stddiff','Xpix','Ypix','FrameRate');
    savefilename=strcat(fileexpnum,'_syncstruct');
    save(savefilename,'syncstruct');
    cd(headpath);
end

%% old (with other fluor data used to save)
%     if m < 0.0004 %approx a 40s
%         cd(sortednames{k})
%         try
%             %load('MC_VidROI.mat');
%             load('MC_Parameters.mat');
%             load('MC_Videocdf2.mat');
%             CinGreen=Cin;
%             ROIs=Ain;
%             load('MC_Video_Redcdf2.mat');%load in the red channel 
%             CinRed=Cin;
%             redtemps{k}=template_red;
% 
%             [dCinGreen,dCinRed]=takeDelta(CinGreen,CinRed,FrameRate, green_bases, red_bases);
%             
% 
% %             syncstruct(synctrial).('green') = CinGreen;
% %             syncstruct(synctrial).('red') = CinRed;
% %             syncstruct(synctrial).('greenbases') = green_bases;
% %             syncstruct(synctrial).('redbases') = red_bases;
% %             syncstruct(synctrial).('allgreen') = all_CinGreen;
% %             syncstruct(synctrial).('allred') = all_CinRed;            
% 
%             syncstruct(synctrial).('dff_green') = dCinGreen;
%             syncstruct(synctrial).('dff_red') = dCinRed;
%             syncstruct(synctrial).('dff_diff') = dCinGreen-dCinRed.*.1; %scale found empirically, related to intensity diff of 2 fluors?
%             [numROIs,triallength]=size(CinGreen);
%             %makes ts by framerate, not real ts, and add trigger time at end (this is how andrew did it, prob not best way, but have no thor ts)
%             syncstruct(synctrial).('dff_ts') =(linspace(1,triallength,triallength)/(FrameRate)) + trig; 
%             syncstruct(synctrial).('thorname') = sortednames{k};
%             syncstruct(synctrial).('thortime') = [trl_time{k}];
%         catch
%             disp('fuck');
%         end
%         cd(headpath);
%     end