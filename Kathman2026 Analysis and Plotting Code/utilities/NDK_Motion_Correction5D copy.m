function [trl_dirs_corr, FOV] = NDK_Motion_Correction5D(fileexpnum, trialcat, zsplit)
% clear all
% zsplit = 0;
% trialcat = 1;
% fileexpnum = '62617-g7_Mar2322_f3e23'; %thor expt #

%set headpath (and reset if in "extracted" output folder)
tic
headpath = pwd; %'/Volumes/Samsung_T5/2p rotaball/21D07/210722'; %directory with all data (should be one for the entire day)
if contains(headpath,'extracted')
    hslsh = strfind(headpath, '/');
    headpath = headpath(1:hslsh(end)-1);
    cd(headpath)
end

%make stacks and concat all trials or mot corr each trial
cd(headpath);
filelist=dir;

if trialcat
    f = waitbar(0,'Concatenating trials...');
else
    f = waitbar(0,'Making stacks and motion correcting...');
end
n = 0; i = 1;
for k=1:numel(filelist)
    dataname=filelist(k).name;
    if (contains(dataname,fileexpnum)&& ~contains(dataname,'corr')) && ~strcmp(dataname(1), '.')%is it a data directory and not a corr directory  
        files = jdir(dataname,'ChanB_0*.tif');
        if ~isempty(files)
            n = n + 1;
            trl_dirs{i} = dataname;
            trl_dirs_corr{i} = fullfile(pwd,[dataname, '_corr']);
            i = i + 1;
        end
    end
end

FOV = [];

trl_dirs = cell(1,n);
Ycat = []; Yredcat = []; 
i = 1; ii = 1;
for k=1:numel(filelist)
    dataname=filelist(k).name;
    if (contains(dataname,fileexpnum)&& ~contains(dataname,'corr')) && ~strcmp(dataname(1), '.')%is it a data directory and not a corr directory        
        files = jdir(dataname,'ChanB_0*.tif');
        if ~isempty(files)
            [Y, Yred] = makeStacks5D_ndk(dataname); %trial directory, progress interval
    
%             maximum z projection for all timepoints
            Y_maxp = squeeze(max(Y,[],3));
            Yred_maxp = squeeze(max(Yred,[],3));

%             % mean z projection for all timepoints
%             Y_maxp = squeeze(mean(Y,3));
%             Yred_maxp = squeeze(mean(Yred,3));
            
            
            if trialcat %concat trials to single stack
                if zsplit
                    Ycat = cat(4, Ycat, Y); %for splits
                    Yredcat = cat(4, Yredcat, Yred);
                else
                    Ycat = cat(3, Ycat, Y_maxp); %for max proj        
                    Yredcat = cat(3, Yredcat, Yred_maxp);
                end
    
                trl_t{i}  = ii:ii+size(Y_maxp, 3)-1;
                ii = ii+size(Y_maxp, 3);
                trl_dirs{i} = dataname;
                trl_dirs_corr{i} = fullfile(pwd,[dataname, '_corr']);
                i = i + 1;
                waitbar(i/n,f,'Concatenating trials...');
            else %motion correct single trial (only doing rigid MC bc non-rigid didn't seem to do any better, but not well tested)
                
                FOV = [size(Yred_maxp,1),size(Yred_maxp,2)];
                options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'init_batch', 200); %, 'iter', 10
                [Mred,shifts_red,template_red] = normcorre_batch(Yred_maxp,options);
                M_final = apply_shifts(Y_maxp,shifts_red,options);
                M_final_red = apply_shifts(Yred_maxp,shifts_red,options);
    
                %find corr folder or create it
                outFolder=strcat(dataname,'_corr');
                if(~exist(outFolder,'dir'))
                    mkdir(outFolder);
                end
    
                %write tiffs
                save(fullfile(outFolder,'MC_Parameters'), 'options', 'shifts_red', 'template_red');
                file_name = fullfile(outFolder,'MC_Video.tif');
                file_name_red = fullfile(outFolder,'MC_Video_Red.tif');
                options.overwrite = true;
                saveastiff(int16(M_final), file_name, options);
                saveastiff(int16(M_final_red), file_name_red, options);
    
                trl_t{i}  = ii:ii+size(Y_maxp, 3)-1;
                ii = ii+size(Y_maxp, 3);
                trl_dirs{i} = dataname;
                trl_dirs_corr{i} = fullfile(pwd,[dataname, '_corr']);
                i = i + 1;
                waitbar(i/n,f,'Making stacks and motion correcting...');
            end
        end
    end
end
% save('220323_cat', 'Ycat', 'Yredcat', 'trl_t', 'trl_dirs')
close(f)

if trialcat
    f = waitbar(0,'Motion correcting...');

    %run mot corr on alltrial stack
    FOV = [size(Ycat,1),size(Ycat,2)];
    options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'init_batch', 200); %, 'iter', 10
    [~,shifts_red_cat,template_red_cat] = normcorre_batch(Yredcat,options);
    Mcat = apply_shifts(Ycat,shifts_red_cat,options);
    Mcat_red = apply_shifts(Yredcat,shifts_red_cat,options);
    
    %split back to trials, run mot corr again, and write tiffs
    for i = 1:n
        for k=1:numel(filelist)
            dataname=filelist(k).name;
            if strcmp(dataname, trl_dirs{i}) %(contains(dataname,fileexpnum)&& ~contains(dataname,'corr'))%is it a data directory and not a corr directory
                Y = Mcat(:,:,trl_t{i});
                Yred = Mcat_red(:,:,trl_t{i});
%                 shifts_red = shifts_red(trl_t{i});
%                 template_red = template_red(trl_t{i});
                FOV = [size(Yred,1),size(Yred,2)];
                options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'init_batch', 200); %, 'iter', 10
                [~,shifts_red,template_red] = normcorre_batch(Yred,options);
                M_final = apply_shifts(Y,shifts_red,options);
                M_final_red = apply_shifts(Yred,shifts_red,options);


                %find corr folder or create it
                outFolder=strcat(dataname,'_corr');
                if(~exist(outFolder,'dir'))
                    mkdir(outFolder);
                end
        
                %write tiffs
                save(fullfile(outFolder,'MC_Parameters'), 'options', 'shifts_red', 'template_red'); %these arent correct bc don't account for first round of normcorre (but close enough for roi demo currently used for?)
                file_name = fullfile(outFolder,'MC_Video.tif');
                file_name_red = fullfile(outFolder,'MC_Video_Red.tif');
                options.overwrite = true;
                saveastiff(int16(M_final), file_name, options);
                saveastiff(int16(M_final_red), file_name_red, options);
        
            end
        end
        waitbar(i/n,f,'Motion correcting...');
    end
    close(f)
end

toc





% Run on each individual z plane (no zprojection)
% load 220323_cat_zsplit.mat
% %run normcorre
% for z = 1:size(Ycat, 3)
%     FOV = [size(Ycat,1),size(Ycat,2)];
%     options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'init_batch', 200, 'iter', 10); %are these parameters really what I want for the motion correction?
%     % % options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'bin_width',200,'max_shift',30,'us_fac',50, 'init_batch', 200, 'iter', 10); %rigid
%     % % options = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',30,'max_dev',3,'us_fac',50,'init_batch',200,'shifts_method','cubic','iter',5);
%     % [M,shifts,template] = normcorre_batch(Ycat,options);
%     this_z = squeeze(Yredcat(:,:,z,:));
%     [Mred,shifts_red,template_red] = normcorre_batch(this_z,options);
%     
%     M_final = apply_shifts(Ycat(:,:,z,:),shifts_red,options);
%     M_final_red = apply_shifts(Yredcat(:,:,z,:),shifts_red,options);
%     
%     outFolder = pwd;
%     save(fullfile(outFolder,['MCcat10x_Parameters_' num2str(z)]), 'options', 'shifts_red', 'template_red');
%     file_name = fullfile(outFolder,['MCcat10x_Video_' num2str(z) '.tif']);
%     file_name_red = fullfile(outFolder,['MC10x_redcat_Video_' num2str(z) '.tif']);
%     
%     options.overwrite = true;
%     saveastiff(int16(M_final), file_name, options);
%     saveastiff(int16(M_final_red), file_name_red, options); 
% end