clear all
%%

allflyinfo

l = 13; %fly line
for e = 20
    %experiment number

    process = 1; %parent flag for motion correct, register rois, sync2struct (set next 4 flags below)
    %fluor and correction off will only run sync2struct
    correction = 1; %run motion correction
    fluor = 1; %calc dff and makes syncstruct 
    new_rois = 0; %new_rois calculates new dff for all trials based on NEW RoiSet.zip in trial 1
    kill_col = []; %nans out each column in vector (1-8), use when columns are coming in and out of stack

plot_trials = 0; %plot time series and trajectory for each trial

%dead subfunctions (can fix plot_rois easily though)
plot_rois = 0;
plot_respsummary = 0;
plot_means = 0;
look_at_dips = 0;
use_roi = 0; %can be a set
windodor_calibration = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fus = strfind(fly{l}{e}, '_');
% headpath = ['/Volumes/T7/' fly{l}{e}(1:fus-1)];
% headpath = ['/Volumes/T7 Shield/preprocessed/' fly{l}{e}(1:fus-1)]; %my black durable ssd
headpath = ['/Volumes/T7/preprocessed/' fly{l}{e}(1:fus-1)]; %lab black ssd
% headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/Data/hAck imaging data/' fly{l}{e}(1:fus-1)]; %lab black ssd


% headpath = ['/Users/kathmn01/Desktop/preprocessed/' fly{l}{e}(1:fus-1)];
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/Data/PID data/' fly{l}{e}(1:fus-1)];
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
%% process and plot trials
if process
    % make stimdata struct
    [data, filenames] = rotaball_importer('date', expt); %don't use date input anymore
    balldata_struct(data, filenames);

    % run motion correct, caiman, import ROIs, calc fluor, plot pooled means for each ROI
%     correction = 0; fluor = 1; new_rois = 1;

%%%%%%%for some reason, RoiSeet.zip for first trial is not working with
%%%%%%%ReadImageJROI function (within dff2 function). cheap fix is to just replace it
%%%%%%%with one from another trial (which all work), but should figure out
%%%%%%%why. Is it getting screwed up while running
%%%%%%%RegisterROI2??
    analyzeCalciumExperiment_newThor5D(headpath,fileexpnum{l}{e},correction,fluor,new_rois);

    % sync with stimdata, interpolate and add bump stats
    cd(headpath)
    close all
%     kill_col = [];
    [syncstruct] = sync2struct(headpath,fileexpnum{l}{e}, kill_col);
end

if plot_rois
    plot_roi(headpath,fileexpnum{l}{e})
end

if plot_trials
    % plot timeseries for each trial
%     clearvars -except fileexpnum{l}{e} headpath use_roi
    cd(headpath)
    us = strfind(fileexpnum{l}{e}, '_');
    slsh = strfind(headpath, '/');
    outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end)+1:end) '_extracted'];
%     outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end-1)+1:us(end)-1) '_extracted'];
    cd(outdirname)
    load([fileexpnum{l}{e} '_syncstruct'])
%     for r = use_roi 
%         pref_roi = r;
%         plottrials_rotaball2p(syncstruct, pref_roi)
%     end
    layers = 0;
    if size(syncstruct(1).bump_amp, 2) > 1
        layers = 1;end
    plottrials_flex_rotaball2p(syncstruct, layers)
%     plottrials_multiroi_rotaball2p(syncstruct, pref_roi{l}(e), columns{l}{e}, prefroi_txt{l}{e}, 1, 1)
    cd(headpath)
end

if plot_respsummary
    
    cd(headpath)
    us = strfind(fileexpnum{l}{e}, '_');
    slsh = strfind(headpath, '/');
    outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end)+1:end) '_extracted'];
%     outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end-1)+1:us(end)-1) '_extracted'];
    cd(outdirname)
    load([fileexpnum{l}{e} '_syncstruct'])
    for r = use_roi 
        pref_roi = r;
        resp_summary_rotaball2p(syncstruct, pref_roi)
    end
    cd(headpath)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot behavior trial means
if plot_means
%     cd(headpath)
%     us = strfind(fileexpnum{l}{e}, '_');
%     slsh = strfind(headpath, '/');
%     outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end)+1:end) '_extracted'];
% %     outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end-1)+1:us(end)-1) '_extracted']; %for 210602 bad naming
%     cd(outdirname)

    cd(headpath)
    us = strfind(fileexpnum{l}{e}, '_');
    slsh = strfind(headpath, '/');
    outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end)+1:end) '_extracted'];
%     outdirname = [headpath(slsh(end)+1:end) '_' fileexpnum{l}{e}(us(end-1)+1:us(end)-1) '_extracted'];
    cd(outdirname)

    load([fileexpnum{l}{e} '_syncstruct'])
%         for zm = [1]
%             for cl = [0, 1]
%                 for dr = [-1, 1, 0]
%                     for spd = [21, 38, inf]
    for zm = [1]
        for cl = [0,1]
            for dr = [0]
                for spd = [inf]                    
                    closed_loop = cl; odor_on = 1; light_on = 0; wind_on = 1; 
                    wind_speed = spd; wind_dir = dr; wind_pos = inf; 
                    add2txt = ''; printfigs = 1; fig_no = 1;
                    sorttrials_rotaball2p(syncstruct, use_roi, closed_loop, odor_on, light_on, wind_on,...
                        wind_speed, wind_dir, wind_pos, add2txt, printfigs, fig_no, zm);
                end
            end
        end
    end
end 
%%
if look_at_dips
    thresh_anlys
end

%% compare max intensities
% CompMaxIntensityDist

%% look at behavior (pin tether)
%     [data, filenames] = rotaball_importer('date', expt);
%     balldata_struct(data, filenames);
%     % balldata_struct_4pooldirs(data, filenames);
%     % rotoball_anlys_v1(noc, [1 1], -1, 'NOC4 OL 75cmsRt', 1)
%     % rotoball_anlys_v1(noc, [1 1], 1, 'NOC4 OL 75cmsLft', 1)
% %     closed_loop = 0; light_on = 0; wind_on = 1; wind_speed = inf; wind_dir = 0; wind_pos = 45; add2txt = 'test'; printfigs = 1; fig_no = 1;
% %     [out, mns] = rotoball_anlys(data, closed_loop, light_on, wind_on, wind_speed, wind_dir, wind_pos, add2txt, printfigs, fig_no);
%     plottrials_rotaballonly(data)
    
%% plot analog in with wind speed and position
if windodor_calibration
    % c(1,:) = [0, 0.4470, 0.7410];
    % c(2,:) = [0.8500, 0.3250, 0.0980];
    % c(3,:) = [0.9290, 0.6940, 0.1250];
    % c(4,:) = [0.4940, 0.1840, 0.5560];
    % c(5,:) = [0.4660, 0.6740, 0.1880];
    % c(6,:) = [0.3010, 0.7450, 0.9330];
    % c(7,:) = [0.6350, 0.0780, 0.1840];
    fly = data;
    figure(1); clf, hold on
    set(gcf, 'Position', [100 100 1000 700])
    for f = 1:length(fly)
%         plot(fly(f).calc_ts, fly(f).ai)
        plot(fly(f).ai)
%         mn(f) = (mean(fly(f).ai(300:600))-0.19)*1+0; %for wind test with [1, 9] time window
% %         mn(f) = (mean(fly(f).ai(420:820))-0.19)*1+0; %for odor test with [5, 13] time window [420:820]
%         plot(fly(f).winddir + randn(1), mn(f), 'ok', 'MarkerSize', 6)
    end 
%     fly = data;
%     figure(1); clf, hold on
%     set(gcf, 'Position', [100 100 1000 700])
%     plot(fly(1).calc_ts, fly(1).ai, fly(2).calc_ts, fly(2).ai, fly(3).calc_ts, fly(3).ai, fly(4).calc_ts,...
%         fly(4).ai, fly(5).calc_ts, fly(5).ai, fly(6).calc_ts, fly(6).ai, fly(7).calc_ts, fly(7).ai)
%     legend('-135deg', '-90deg', '90deg', '135deg', '0deg', '45deg', '-45deg', 'Location', 'northwest')
    base = mean(fly(f).ai(10:50));
    plot([-180 180], [base, base], '--k')
    xlabel('wind pos (deg, -/+ = rt/lft)'), ylabel('volts')

    alldirs = [fly.winddir];
    udirs = unique(alldirs);
    i = 1;
    for d = udirs(udirs < 0)
        bias(i) = (mean(mn(alldirs == d))-base)/(mean(mn(alldirs == -d))-base);
        i = i + 1;
    end
    title(['AT FLY, COMP ON, VAC ON; Biases: ' num2str(bias) ', wind speed = 0.45 L/m, odor speed = 0.4 L/m, n = ' num2str(mode(histc(alldirs, udirs))) ', -- = windoff'])
    % legend('0.2 L/m odor', '0.4 L/m odor', '0.6 L/m odor')
    % print_fig(1, gcf, 'Windbalance-atfly-compoff-speed0d45', 'png', 1);
end
% %test legend error (this is from the help page, still doesn't work, wtf)
% figure;
% x = 0:.2:12;
% plot(x,besselj(1,x),x,besselj(2,x),x,besselj(3,x));
% legend('First','Second','Third','Location','NorthEastOutside')

%% check ball calibration
% ttl = {'YawTurn-LftSlow', 'YawTurn-LftFast', 'YawTurn-RightSlow', 'YawTurn-RightFast',...
%     'PitchTurn-BkwdSlow', 'PitchTurn-BkwdFast', 'PitchTurn-FwdSlow', 'PitchTurn-FwdFast', ...
%     'RollTurn-LftSlow', 'RollTurn-LftFast', 'RollTurn-RightSlow', 'RollTurn-RightFast'};
% for t = 1:size(fly,2)
%     figure(t); clf; smth = 15;
%     subplot(3,1,1); plot(smoothdata(data(t).calc_deltaz, 'movmean', smth)); ylim([-10 10]); ylabel('deltaZ'); 
%     title([data(t).filename(1:end-4), ' ' ttl{t}], 'Interpreter', 'none')
%     subplot(3,1,2); plot(smoothdata(data(t).calc_deltapitch, 'movmean', smth)); ylim([-10 10]); ylabel('deltaPitch');
%     subplot(3,1,3); plot(smoothdata(data(t).calc_deltaroll, 'movmean', smth)); ylim([-10 10]); ylabel('deltaRoll');    
%     print_fig(1, gcf, [data(t).filename(1:end-4), '_' ttl{t}], 'png', 1); 
% end
end
%%
%%
% %MB052B
% fileexpnum{1}{1} = 'MB052B_Feb2822_f3e8';
% fileexpnum{1}{2} = 'MB052B_Feb2822_f3e11';
% fileexpnum{1}{3} = 'MB052B_Feb2822_f4e12';
% fileexpnum{1}{4} = 'MB052B-g7_Mar2222_f1e4';
% fileexpnum{1}{5} = 'MB052B-g7_Mar2922_f1e6';
% fileexpnum{1}{6} = '';
% fileexpnum{1}{7} = '';
% %FB5AB
% fileexpnum{2}{1} = '21D07_Mar0322_f1e3'; %actually 220301
% % fileexpnum{2}{1} = '21D07_Mar0322_f2e17'; %actually 220301 %ol only, poor gcamp expression 
% fileexpnum{2}{2} = '21D07_Mar0122_f3e23';
% fileexpnum{2}{3} = '21D07-g7_Mar1422_f3e23';
% fileexpnum{2}{4} = 'FB5AB-g7_Mar2222_f2e14';
% fileexpnum{2}{5} = 'FB5AB-g7_Mar2322_f1e3';
% fileexpnum{2}{6} = '21D07-g7_Mar2822_f1e3';
% fileexpnum{2}{7} = '21D07-g7_Mar2822_f1e4';
% fileexpnum{2}{8} = ''; %plume
% fileexpnum{2}{9} = ''; %plume
% fileexpnum{2}{10} = ''; %plume
% %HDELTAC
% % fileexpnum{3}{1} = '62617_Mar0122_f5e44'; %saline leak
% % fileexpnum{3}{1} = '62617-g6_Mar1422_f1e4'; %poor health, behavior, expression
% fileexpnum{3}{1} = '62617-g7_Mar1522_f1e3';
% fileexpnum{3}{2} = '62617-g6_Mar1522_f2e13';
% fileexpnum{3}{3} = '62617-g6_Mar2122_f1e3'; %bad expression???
% fileexpnum{3}{4} = '62617-g6_Mar2122_f2e13'; %bad expression???
% fileexpnum{3}{5} = '62617-g7_Mar1522_f3e25'; %actually 220321
% fileexpnum{3}{6} = '62617-g7_Mar2322_f3e23';
% fileexpnum{3}{7} = '62617-g7-Apl1822_f1e4';
% fileexpnum{3}{8} = '62617-g7-Apl1922_f1e3'; %bad expression???
% 
% fileexpnum{3}{9} = '62617-g7-Jun1422_f1e5'; %cwoo, plume
% 
% fileexpnum{3}{10} = '0818_5'; %cwoo (blowup)
% fileexpnum{3}{11} = '0829_13'; %cwoo, cwoo/owoo (weak)
% fileexpnum{3}{12} = 'vt062617-gc7f-tdtom(II)-Sept2922_f1e3'; %cwoo, heading, cwoo, heading
% fileexpnum{3}{13} = '1004_13'; %plume, cwoo, heading, cwoo
% fileexpnum{3}{14} = '1011_13'; %(weak)
% fileexpnum{3}{15} = '1012_3'; %(no resp?)
% fileexpnum{3}{16} = '1012_13'; %cwoo, plume, blanks, cwoo, cwoo/owoo
% fileexpnum{3}{17} = '1017_3'; %cwoo
% fileexpnum{3}{18} = '1017_13-14'; %cwoo, cwoo/owoo
% fileexpnum{3}{19} = '1018_13'; %cwoo, cwoo/owoo
% fileexpnum{3}{20} = '1018_23'; %cwoo, cwoo/owoo
% fileexpnum{3}{21} = '1018_34'; %plume, cwoo
% fileexpnum{3}{22} = '1019';
% fileexpnum{3}{23} = '';
% fileexpnum{3}{24} = '';
% fileexpnum{3}{25} = '';
% fileexpnum{3}{26} = '';
% fileexpnum{3}{27} = '';
% fileexpnum{3}{28} = '';
% fileexpnum{3}{29} = '';
% fileexpnum{3}{30} = '';
% fileexpnum{3}{31} = '';
% 
% %FC
% fileexpnum{4}{1} = '52g12-gc6f-tdtom-Oct1022_f1e3';  %cwoo, cwoo/owoo
% fileexpnum{4}{2} = '1011_3'; %cwoo, cwoo/owoo
% 
% %HDELTAK
% fileexpnum{5}{1} = '';  %cwoo, cwoo/owoo
% fileexpnum{5}{2} = ''; %cwoo, cwoo/owoo
% 
% 
% %RGECO/62617; GCaMP7f/21D07
% fileexpnum{6}{1} = '1003_3,5'; %cwoo (no rgeco seen)
% 
% % fly{3}{6} = '220323_23'; age{3}(6) = 13; zt{3}(6) = 10; pref_roi{3}(6) = 7; columns{3}{6} = [4, 2, 1, 3, 3, 5, 7, 8];  prefroi_txt{3}{6} = 'rtFB'; bump_thresh(3,6) = 3; %12='cntr collat?';

%% testbed
% fig_no = 1; i = 1;
% for s = [4 8 25 50 76 101]
%     for p = [0 45 90 135 180]
%         for l = 0
%             for d = -1   % wind_dir: -1=only plot wind right, 1=only plot wind left, 0=plot all dirs
%                 closed_loop = l;
%                 wind_speed = s; %[8, 25, 50, 76, 101 cm/s] = [0.1, 0.3, 0.6, 0.9, 1.2 L/m]
%                 wind_pos = p; %[0 45 90 135 180]
%                 wind_dir = d; %dir_select: -1=only plot wind right, 1=only plot wind left, 0=plot all dirs0
%                 light = 1;
%                 wind = 1; 
%                 add2txt = '';
%                 printfigs = 1;
% 
%                 out{i} = rotoball_mns(fly, closed_loop, light, wind, wind_speed, wind_dir, wind_pos,...
%                     add2txt, printfigs, fig_no);
%                 fig_no = fig_no + 1;
%                 i = i + 1;
%             end
%         end
%     end
% end

% ws = [];
% for i = 1:length(eval(fly))
%     eval(['ws = [ws; ' fly '{i}.wind_speed];'])
% end
% unique(ws)
% 
% wp = [];
% for i = 1:length(eval(fly))
%     eval(['wp = [wp; ' fly '{i}.wind_pos];'])
% end
% unique(wp)
