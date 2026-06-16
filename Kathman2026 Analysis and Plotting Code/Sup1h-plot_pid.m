clear all
close all


allflyinfo

l = 9; %fly line
e = 2; %2 = 15s block, 7 = 10ms pulses; 
print_figs = 0;
trunc4meandff = 1;

if e == 2
    add2txt = '15sblock_pid_n5';
else
    add2txt = '10mspulses_pid_n10';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fus = strfind(fly{l}{e}, '_');
% headpath = ['/Volumes/T7/' fly{l}{e}(1:fus-1)];
% headpath = ['/Users/kathmn01/Desktop/preprocessed/' fly{l}{e}(1:fus-1)];
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/Data/PID data/' fly{l}{e}(1:fus-1)];
headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/pid-anemom/' fly{l}{e}(1:fus-1)];
cd(headpath)

% headpath = pwd; %'/Volumes/Samsung_T5/2p rotaball/21D07/210722'; %directory with all data (should be one for the entire day)
if contains(headpath,'extracted')
    hslsh = strfind(headpath, '/');
    headpath = headpath(1:hslsh(end)-1);
    cd(headpath)
end
eeees = strfind(fileexpnum{l}{e}, 'e');
us = strfind(fileexpnum{l}{e}, '_');
if us(end) > eeees(end)
    expt = str2double(fileexpnum{l}{e}(eeees(end)+1:us(end)-1)); %pybmt expt #
else
    expt = str2double(fileexpnum{l}{e}(eeees(end)+1:end)); %pybmt expt #
end

[data, filenames] = rotaball_importer('date', expt); %don't use date input anymore
balldata_struct(data, filenames);
% clearvars -except expt print_figs add2txt
load(['E' num2str(expt) '.mat'])


j = 1;
maxframes = 0;
% use_trials = [3:11, 18];
use_trials = 2:7; %[6 7 8 9 10 11 12 13 14 15]; %[1:length(data)];

% maxfrm = 0;
% for i = 1:length(use_trials)
%     maxfrm = max(maxfrm, max(data(use_trials(i).frame)));
% end
st = 28.05; %stim on time

figure(1); clf; 
for i =  use_trials%3:length(data)
    hold on
    x = downsample(data(i).calc_ts, 5);
    plot(x-st, downsample(data(i).ai, 5), 'k')
    plot(x-st, downsample(data(i).wind*.03+.55, 5), 'k')
    plot(x-st, downsample(data(i).odor*.035+.5, 5), 'k')

%     this_aligned = nan(1,maxfrms)
%     alltrl = [alltrl, nan()]
%     keyboard
%     j = j + 1;
    if max(data(i).frame) > maxframes
        maxframes = max(data(i).frame);
    end
    
end
% ylim([.08 .2])
% xlim([16 55.19])


pid = nan(length(use_trials), maxframes);
% 
figure(1);% clf; hold on
j = 1;
for i = use_trials
    pid(j,data(i).frame) = data(i).ai;
    j = j + 1;
end
mnpid = nanmean(pid);
ts = .01:.01:maxframes*.01;
plot(ts-st, mnpid, 'r')
if trunc4meandff
    xlim([-12 27])
end

subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/';
if print_figs
    exportgraphics(gcf, [subpath add2txt '.png'])
    print_fig(1, gcf, [subpath add2txt], 'eps', 1);   
    
end