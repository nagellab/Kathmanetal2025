function dff2(foldername,Ypix,Xpix)
%% load files 
% Folder with stack images of diferent sessions
%foldername='Z:\Andrew\2-photon imaging\vt019532\vt019532_f3e4_010_corr'; %path to image folder
cd(foldername);
listfiles = dir('*.tif'); %will look for all the tif files
% !!! Files are listed in alphabetical order !!!
for i= 1:length(listfiles)
files{i}=fullfile(foldername, listfiles(i).name);
end
% Import ROI
FOV = [Ypix,Xpix]; % Image resolution (change with image)%should figure this out from the image since it will change

% ROI_file='RoiSet.zip'; %path to ROI file (.zip)
ROI_file = dir('*.zip');
% Need ReadImageJROI.m script
if length(ROI_file) < 2
[a,ROI] = ReadImageJROI(ROI_file.name,[Ypix,Xpix]);%need to reformat a using poly2mask
else
[a,ROI] = ReadImageJROI(ROI_file(2).name,[Ypix,Xpix]);%need to reformat a using poly2mask
end
[a,ROI] = fixROIs(a,ROI);
% Create Structure
input.foldername=foldername;
input.image=files;
input.list=listfiles;
input.ROI.file=ROI;
input.ROI.a=a;
input.param.FOV=FOV;

%% Set parameters soma 
K = 1;                                           % number of components to be found
tau = 4;                                          % std of gaussian kernel (size of neuron)
p = 2;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;                                  % merging threshold
options = CNMFSetParms(...
                       'd1',Ypix,'d2',Xpix,...                        % dimensions of datasets
                       'search_method','dilate','dist',3,...       % search locations when updating spatial components
                       'deconv_method','constrained_foopsi',...    % activity deconvolution method
                       'method','dual',...                         %hopeless prayer AM %prevents the code from breaking cvx or whatever doesn't work
                       'temporal_iter',2,...                       % number of block-coordinate descent steps
                       'fudge_factor',0.98,...                     % bias correction for AR coefficients
                       'merge_thr',merge_thr,...                    % merging threshold
                       'gSig',tau...
                       );
% Create Structure                
input.param.options=options;
input.param.tau=tau;
input.param.p=p;
input.param.merge_thr=merge_thr;
                  
%% CNMF (constrained nonnegative matrix factorization)
%ndk: just pulls raw pixels (nnz from corr max proj) within ROI, no other
%part of functions used (eg tau, merg, deconv, etc)
input.refine=0; %refine components 
%Session to process
%session=1; 
%OR
%all session
session=1:size(files,2);

for i=1:length(session)
% [output]= CNMF_noparpool(input, session(i)); 
[output]= CNMF_noparpoolAM(input, session(i)); 
end





