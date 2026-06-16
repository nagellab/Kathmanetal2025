clear all
% close all

headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/all flies/'];
cd(headpath)
allflies_filenames

savepath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/';
savepath = '/Users/kathmn01/Desktop/';
% use_flies = [21, 22, 25, 33, 35, 38]; %can be vector
% fly_yl = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]; %matches use_flies (vector or value)
% dff_color = [.05 .5; .05 .7; .05 .5; .05 .5; .05 1; .05 1];
% sdms = [.7, .5, .7, .5, .7, .7];
% which_trials{1} = [23];
% which_trials{2} = [8];
% which_trials{3} = [17, 34];
% which_trials{4} = [16, 17];
% which_trials{5} = [20];
% which_trials{6} = [34, 35];

%bump doesn't move trial examples
% use_flies = [21, 22, 25, 33, 35, 38]; %can be vector
% fly_yl = [1.2, 1.2, 1.2, 1.2, 1.2, 1.2]; %matches use_flies (vector or value)
% dff_color = [.05 .5; .05 .7; .05 .5; .05 .5; .05 1; .05 1];
% sdms = [.7, .5, .7, .5, .7, .7];
% which_trials{1} = [23];
% which_trials{2} = [8];
% which_trials{3} = [17, 34];
% which_trials{4} = [16, 17];
% which_trials{5} = [20];
% which_trials{6} = [34, 35];

% % add2txt = '_fig2ex';

%no corr with pos and heading example
% use_flies = [25];
% which_trials{1} = [1:4, 7, 11:2:17, 18:2:48, 61, 64];
% use_flies = [21];
% which_trials{1} = [1:9, 16:25, 42];
% which_trials{1} = [17, 30, 32];
% fly_yl = 1.2;
fig_no = 100;
% dff_coloraxis = [0.05 .7];

%bump doesn't shift examples
% use_flies = 33;
% which_trials{1} = 16;
% fly_yl = .75;
% dff_color = [.05 .5];
% sdms = .65;

use_flies = 25;
which_trials{1} = 34;
fly_yl = .75;
dff_color = [.05 .5];
sdms = .7;

% %bump does shift examples
% use_flies = 32;
% % use_flies = 154;
% % which_trials{1} = 42;
% which_trials{1} = [1];
% fly_yl = 1.3;
% dff_color = [0.1 .8];
% sdms = .5;

% %bump persists examples
% use_flies = 100;
% which_trials{1} = [8, 16, 24];
% fly_yl = 1;
% dff_color = [0.2 .6];

% %bump persists examples
% use_flies = 38;
% which_trials{1} = [39];
% fly_yl = 1.1;
% dff_color = [.1 .7];
% sdms = .5;
% 
%bump persists examples (paper example!!!)
use_flies = 25;
% which_trials{1} = [15, 18, 22, 38, 61, 63];
which_trials{1} = [22];
fly_yl = .8;
dff_color = [.1 .5];

% % owoo
% use_flies = [25];
% which_trials{1} = [29, 37]; %schematic/first plume examples
% fly_yl = .5;

% %plume fig examples
% use_flies = [27];
% which_trials{1} = [15, 13]; %schematic/first plume examples
% fly_yl = 1.2;

%more plume examples
% use_flies = [27, 36, 48, 49];
% fly_yl = [1.2, 1.2, 2.2, 2.5];
% which_trials{1} = [15, 13];
% which_trials{2} = [14];
% which_trials{3} = [9, 11];
% which_trials{4} = [9];


% %wind shift trial %remember to add lines at 30 and use fictrac heading
% use_flies = 117; %[57];
% fly_yl = 2;
% which_trials{1} = 36; %15; %[6, 18, 26];
% % which_trials{2} = 15;
% % use_flies = [25];
% % fly_yl = 0.8;
% % which_trials{1} = 15;

add2txt = '';
print_flag = 0;
file_type = 'eps';

sdm = .7;

for f = 1:length(use_flies)
    data = [];
    load(filename{use_flies(f)})
    data = syncstruct;

    sdm = sdms(f);
    
    %find fly stat ranges for axes and bump threshold
    allamp_thr = []; use_trials = [];
    for t = 1:length(data)
        nonan = ~isnan(data(t).bump_amp(:,1));
    
        allamp(t) = quantile(data(t).bump_amp(~isnan(data(t).bump_amp(:,1))), .90); 
        allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) %plume
            use_trials = [use_trials, t];
%         end
    end
    if isempty(which_trials{f})
        which_trials{f} = use_trials;
    end
    
    %     ylm = [-.05 nanmean(allamp)*1.65];
    ylm = [0 nanmax(allamp)];
    %     ylm = [0 0.52];
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;
    
    
    dff_yaxis = [0 fly_yl(f)];
%     dff_coloraxis = [0 fly_yl(f)*.85];
%     dff_coloraxis = [0.05 .5];
%     dff_coloraxis = [0 1.4];
    dff_coloraxis = dff_color(f,:);
    
    
    fig_no = plottrials_flex(data, savepath, use_flies(f), which_trials{f}, dff_yaxis, dff_coloraxis, base_thr, add2txt, print_flag, file_type, 'jet', fig_no);
    fig_no = fig_no + 1;
end



