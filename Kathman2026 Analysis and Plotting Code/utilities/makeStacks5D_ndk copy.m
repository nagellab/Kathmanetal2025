function [Ygreen, Yred] = makeStacks5D_ndk(folder)
% This function will make 4-dimensional stacks from 5-dimensional data.
% 5D data: x, y, z, t, c
% Volume imaging time series with two color channels.

%make 5D stacks
Ygreen = []; Yred = [];
chan = {'ChanA_0*.tif', 'ChanB_0*.tif'};
for c = 1:2
    files = jdir(folder,chan{c});
    N = length(files);
    
    %make blank
    r = imread(fullfile(folder,files(1).name));
    t = str2num(files(end).name(end-6:end-4));
    z = str2num(files(end).name(end-10:end-8));
    Y{c} = zeros(size(r,1),size(r,2),z,t,'uint16');

    %fill maxtrix
    for i=1:N    
        r = imread(fullfile(folder,files(i).name));
        t = str2num(files(i).name(end-6:end-4));
        z = str2num(files(i).name(end-10:end-8));
        Y{c}(:,:,z,t) = r;
    end
end
Ygreen = Y{1};
Yred = Y{2};



% 
%     if(nargin<2)
%         pctUpdate = [];
%         flag = false;
%     else
%         flag = true;
%     end
% 
%     
%     files = jdir(folder,'ChanA_0*.tif');
%     
%     %make chan A blank
%     if(~isempty(files))
%         r = imread(fullfile(folder,files(1).name));
%         t = str2num(files(end).name(end-6:end-4));
%         z = str2num(files(end).name(end-10:end-8));
%     else
%         r = NaN(512,512);
%         t = NaN(1,1);
%         z = NaN(1,1);
%     end
%     N = length(files);
%     Ygreen = zeros(size(r,1),size(r,2),z,t,'uint16');
%     
%     %set progress markers
%     if(flag==1)
%         npoints = ceil(100/pctUpdate);
%         checkpoints = ceil(linspace(0,N,npoints+1));
%         checkpoints = checkpoints(2:end);
%     else
%         checkpoints = [];
%     end
% 
%     %make chan A stack
%     count = 1;
%     for i=1:N
%         if(flag && i>=checkpoints(count))
%             count = count+1;
%             fprintf('%d%%...',ceil(100*i/N));
%         end
%         r = imread(fullfile(folder,files(i).name));
%         t = str2num(files(i).name(end-6:end-4));
%         z = str2num(files(i).name(end-10:end-8));
%         Ygreen(:,:,z,t) = r;
%     end
%     
%     %make chan B blank
%     files = jdir(folder,'ChanB_0*.tif');
%     if(~isempty(files))
%         r = imread(fullfile(folder,files(1).name));        
%         t = str2num(files(end).name(end-6:end-4));
%         z = str2num(files(end).name(end-10:end-8));
%     else
%         r = NaN(512,512);
%         t = NaN(1,1);
%         z = NaN(1,1);
%     end
%     N = length(files);
%     Yred = zeros(size(r,1),size(r,2),z,t,'uint16');
%     
%     %set progress markers
%     if(flag==1)
%         npoints = ceil(100/pctUpdate);
%         checkpoints = ceil(linspace(0,N,npoints+1));
%         checkpoints = checkpoints(2:end);
%     else
%         checkpoints = [];
%     end
% 
%     %make chan B stack
%     count = 1;
%     for i=1:N
%         if(flag && i>=checkpoints(count))
%             count = count+1;
%             fprintf('%d%%...',ceil(100*i/N));
%         end
%         r = imread(fullfile(folder,files(i).name));        
%         t = str2num(files(i).name(end-6:end-4));
%         z = str2num(files(i).name(end-10:end-8));
%         Yred(:,:,z,t) = r;
%     end
%     
%     % maximum z projection for all timepoints
%     Ygreen_split = Ygreen;
%     Yred_split = Yred;
%     Ygreen = max(Ygreen,[],3);
%     Yred = max(Yred,[],3);
%     
%     Ygreen_split = squeeze(Ygreen_split);
%     Yred_split = squeeze(Yred_split);
%     Ygreen = squeeze(Ygreen);
%     Yred = squeeze(Yred);
end
