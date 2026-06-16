function [output]= CNMF_noparpoolAM(input, session)  
%% load structure
foldername = input.foldername;
files = input.image;
listhdf5 = input.list;
ROI = input.ROI.file;
a = input.ROI.a;
FOV = input.param.FOV;
options = input.param.options;
tau = input.param.tau;
p = input.param.p;
merge_thr =input.param.merge_thr;
%Session to process
nam=files{session}; 
%sframe=1;
%num2read=5000;
%Y = bigread2(nam,sframe,num2read);
Y = bigread2(nam);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double
[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;
%% Data pre-processing
[P,Y] = preprocess_data(Y,p); %interpolate missing pixels
%% fast initialization of spatial components using greedyROI and HALS
%Y and ROI are my friends (and hence Yr)
%Ain is to ROI as Yr is to Y
%what is Cin


Cn_max =  correlation_image_max(Y);
Cn = Cn_max;%basically a max project but using correlations
%Use ROI imported ROI for Ain
Ain= sparse(reshape(ROI,FOV(1)*FOV(2),[]));
for k=1:size(ROI,3)%Cin is just the raw fluorescence DUMB NAMING 
    ROI_bin{k}=ROI(:,:,k);
    Cin(k,:) = ROI_bin{k}(:).'*reshape(Y,[],T)/nnz(ROI_bin{k});%nnz of ROI bin returns a scalar which is the number of pixels in the ROI (same as sum in both directions)
end 
% Yr = reshape(Y,d,T);
% res = reshape(Y,[],T) - Ain*Cin; %subtracting the spatial components from the background I think
% fin = mean(res); %what is this
% for nmfiter = 1:20 %what does this do
% bin = max((res*fin')/(fin*fin'),0);%sounds like B is the background    
% fin = max((bin'*bin)\(bin'*res),0);%right division you hate to see it
% end
%% update spatial components
%[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);
%%turn off because it remaps it all weird and I'm not a fan. Just take my
%%badly drawn hand ROIs please 
%% update temporal components
% P.p = 2;    
% [C,f,P,S,YrA] = update_temporal_components(Yr,Ain,bin,Cin,fin,P,options); %have to update to get the p component right (I do not like this)
%% classify components
%[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components_noparpool(Y,A,C,b,f,YrA,options);
% run GUI for modifying component selection (optional, close twice to save values)
% run_GUI = true;
% if run_GUI
%    Coor = plot_contours(A,Cn,options,1); close;
%    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
%    options = GUIout{2};
%    keep = GUIout{3};    
% end

%% refine estimates excluding rejected components
% if input.refine==true %not true by default so wont enter AM
% %[A2,b2,C2] = update_spatial_components(YrA,C,f,[A,b],P,options);
% %[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,P,options);
% elseif input.refine==false
% A2=Ain;
% b2=bin;
% C2=C;
% P2=P;
% f2=f;
% S2=S;
% YrA2=YrA;
% end
%% Extract DF/F 
%[C_df,~] = extract_DF_F(Yr,A2,C2,P2,options); %would be great if this ran without P2
%New df/f extract
alpha=0.05;
%[expDffMedZeroed, expDff,dff,F,bf,dfc] = dff_extract_3(YrA2, A2,C2, b2,f2,alpha);
%% detrend fluorescence and extract DF/F values %also skip? AM
%options.df_window = 1000; 
%[F_dff,F0] = detrend_df_f(A2,b,C2,f2,YrA2,options);
%% deconvolve data %skip this??? AM
%nNeurons = size(F_dff,1);
%C_dec = zeros(size(F_dff));
%S = zeros(size(F_dff));
%kernels = cell(nNeurons,1);
%min_sp = 3;    % find spikes resulting in transients above min_sp x noise level
%for i = 1:nNeurons
%    [C_dec(i,:),S(i,:),kernels{i}] = deconvCa(F_dff(i,:), [], min_sp, true, false, [], 20, [], 0);
%end
%% display components
%[Coor,json_file] = plot_contours(A2,Cn_max,options,1);
%plot_components_GUI(Yr,A2,C2,b2,f2,Cn,options)
%center for ROI
%center= com(A2,d1,d2);
%Save
disp(['Saving : ' listhdf5(session).name]);
tic;
%trying to go down a directory when it can't?
%save(fullfile([nam(1:end-6),'Cdf']),'C_df','expDffMedZeroed','dff'); %skip saving S AM ,'S','C_dec','F_dff'
save(fullfile([nam(1:end-4),'cdf2']),'Cin','Ain'); %Ain is a lot smaller than ROI so save that - reshape later
%save(fullfile([nam(1:end-6),'ROI']),'center','A2','C2','b2','f2','Cn_max', 'options'); %take out ROIvars 'ROIvars','Coor' ,'json_file',, 'keep'
disp(['Done saving : ' listhdf5(session).name]);
output=[];
%output.expDffMedZeroed= expDffMedZeroed;
%output.C_df= C_df;
%output.center=center;
%output.Coor=Coor;
%output.keep=keep;
%output.Cn_max=Cn_max;
%output.A2=A2;





