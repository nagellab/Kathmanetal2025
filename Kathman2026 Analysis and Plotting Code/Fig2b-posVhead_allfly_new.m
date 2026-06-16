%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2024 paper/sandbox/exploratory figs/bump postion/';
cd(headpath)

allflies_filenames

%%
% use_flies = [21:27, 30:40]; %[21:40, 91:93];%[21:26, 33, 37]; %[21:26,30, 32:35, 37:39]; %21:39;%[6,8:10];%1:17; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging


%%
% use_flies = [21:26,30, 32:35, 37:39]; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging
% use_flies = [21:27, 30:40]; sdm = 0.5;%[21:40, 91:93];%[21:26, 33, 37]; %[21:26,30, 32:35, 37:39]; %21:39;%[6,8:10];%1:17; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging
% use_flies = 41:51; sdm = 0; %sdm not used in colornorm
% use_flies = [33]; sdm = 1;
% use_flies = [22, 23, 25, 33]; sdm = 1;
use_flies = [21:26,30, 32:35, 37:39]; sdm = .5;

% use_trials{1} = [1:3,5,6,8,9];
% use_trials{2} = [1,2,15,21,25,28,29,33,47,49,51];
% use_trials{3} = [2:4, 7, 11, 13, 15, 17 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 42, 54, 64];
% use_trials{4} = [16, 21:26, 28, 30];

print_figs = 0;
add2txt = '';

all_pos = [];
all_head = [];
all_cc = [];
allfly_cc_mns = [];

c = []; 
j = 1;
fig_no = 1;
figure(fig_no); clf; hold on
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;


    use_trials{f} = [];
    allamp_thr = [];
    for t = 1:length(data)
        nonan = ~isnan(data(t).bump_amp(:,1));
        odor = data(t).odor(nonan);
        idyl = odor >= 0.5; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        oncrossi = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)        

        allamp(t) = quantile(data(t).bump_amp(~isnan(data(t).bump_amp(:,1))), .95); 
        allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
        if data(t).closed_loop && length(oncrossi) < 3 && max(odor) > 0 %data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ~isempty(data(t).bump_amp) && length(oncrossi) < 2 && max(odor) > 0%cwoo
%         if length(oncrossi) < 2 && max(odor) > 0
            use_trials{f} = [use_trials{f} t];
        end
    end
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;


    all_t = 0;
    all_tt = 0;
    all_dff = [];
    all_pos = [];
    all_head = [];
    all_odor = [];
    allfly_cc = [];
    if ~isempty(use_trials{f})
        for t = use_trials{f}
            nonan = ~isnan(data(t).bump_amp(:,1));
    
            %find stim times before trial sorting bc used to rule out plume trials
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
    

    %         nonan = ~isnan(data(t).bump_amp(:,1));
            odor = data(t).odor(nonan);
            wind = data(t).wind(nonan);
            ampthresh = .5;%0.125;
            idyl = odor >= ampthresh; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(odor(idy-1)<ampthresh);
    
    
            dff = data(t).bump_amp(nonan, 1);
            speed = data(t).calc_speed(nonan);
            fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
            fps_im = median(diff(data(t).dff_ts));
            if data(t).closed_loop
                heading = data(t).calc_windpos(nonan);
            else
                heading = data(t).calc_heading(nonan);
            end
            heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
            uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 100);
            tt = data(t).calc_ts(nonan);
            col = data(t).bump_col(nonan, 1);
            com = data(t).bump_com(nonan, 1);
            pos = smoothdata(data(t).bump_col(nonan), 'movmean', 100); %200);
            avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 50);
            abs_avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 50);
            uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
            fvel = data(t).calc_deltapitch(nonan)*9.52/2;
            fvel = smoothdata(fvel, 'movmean', 30); 

%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 5, base_thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
%             bumpi = find(bump);

            im = data(t).dff_green(1:8,:);
            ts = data(t).dff_ts;
            fs = size(im, 2)/max(ts);            
            [pos, uwpos, ~, thr] = bumpposfinder(im, ts, tt, base_thr, sdm); %thr is min of expt thresh (base_thr) and trial thresh

%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 3, thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
            [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 5, base_thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
            bumpi = find(bump);

%                 imagesc(data(t).dff_ts, (1:8), im, [0 1])
%                 colormap('gray')
%                 scatter(tt(bumpi), uwpos(bumpi), 2, 'g', 'filled')
%                 scatter(tt(bumpi), wrapTo2Pi(uwheading(bumpi))*1.1+1, 2, 'r', 'filled')


            if ~isempty(bumpi)% && length(bumpi) > 200 %no need for length limit, already in bumpfinder
                all_pos = [all_pos, mean(pos(bumpi))];
                all_head = [all_head, mean(heading(bumpi))];
                r = corrcoef(uwpos(bumpi), uwheading(bumpi));
                all_cc = [all_cc, r(2)];
                allfly_cc = [allfly_cc, r(2)];
            end

% 
%             im = data(t).dff_green(1:8,:);
%  
%             ttt = data(t).dff_ts;
% %             all_dff = [all_dff, [im, ones(8, 10)]];
% %             all_t = [all_t, [all_t(end)+ttt]]; 
% %             all_t = [all_t, (all_t(end)+1/fps_im):1/fps_im:(all_t(end)+10/fps_im)];
% %             all_pos = [all_pos, pos(bumpi)];
% %             all_pos_t = [all_pos_t, tt(bumpi)+all_pos_t(end)+11];
%             all_head = [all_head, heading'];
%             all_tt = [all_tt, tt'+all_tt(end)+1/fps];
% 
%             all_dff = [all_dff, im];
%             all_t = [all_t, [all_t(end)+ttt]]; 
% 
%             all_odor = [all_odor, odor'];
        end
    end



if length(all_pos) > 4
    figure(100);
    subplot(4,4,j)
    scatter(all_head, all_pos, 80, 'filled', 'k'); 
    r = corrcoef(all_head, all_pos);
    title(['fly' num2str(use_flies(f)) ', cc=' num2str(r(2))])
    % mdl = fitlm(all_head,all_pos);
    % plot(mdl)
    ylim([0 8])
    xlim([-180, 180])

    c = [c, r(2)];
    j = j + 1;
end



%     if ~strcmp(plot_what, 'dff')
%         plot(tt, eval(plot_what), 'k'); hold on
%         if plot_stimlines
%             plot(tt, wind.*.25*abs(diff(ylm))+ylm(2)*1.3, 'LineW', 2, 'Color', [.5 .5 .5]) %wind.*30+160
%             plot(tt, odor.*.25*abs(diff(ylm))+ylm(2)*1.3, 'k', 'LineW', 2)
%             ylim([ylm(1)-abs(diff(ylm))*.1 ylm(2)+abs(diff(ylm))*.6]) %%%%%change these scalers to fix ylim; fvel: .5, .6; dff: .2, .6
%         end
%         
% %                     ylim([ylm(1)-abs(diff(ylm))*.5 ylm(2)+abs(diff(ylm))*.6])
%         ylabel(plot_what)
%         oi = find(odor);
%         title(['t=' num2str(t), ',thr=' num2str(thrsh, 2), 'avel=' num2str(mean(avel(oi)))])
%     else
%         ylm = [0 nanmax(nanmax(im)) ];
%         imagesc(data(t).dff_ts, (1:8), im, ylm)
% 
%         colormap('gray')
% %             sdff = smoothdata(dff, 'sgolay', 300);
% %             plot(tt, sdff, 'k')
%         plot(tt, wind.*1.25+9, 'LineW', 2, 'Color', [.5 .5 .5])
%         plot(tt, odor.*1.25+9, 'k', 'LineW', 2)
%         ylim([1 11])
%         colorbar
%     end
% 
%     if plot_bumppos
%         scatter(tt(bumpi), pos(bumpi)*360/8-202.5, 2, 'g', 'filled')
%     end  



%     sgtitle(['Fly #' num2str(use_flies(f))])


% figure(200); hold on;
% scatter(1+randn(1)*.05, nanmean(allfly_cc), 80, 'filled', 'k')
% ylim([-1 1])
% xlim([0 2])
allfly_cc_mns = [allfly_cc_mns, nanmean(allfly_cc)];


end

if print_figs
    print_fig(1, gcf, [subpath 'bumpposVheading_cwoo_allfliesw5trials', add2txt], 'eps', 1); 
    exportgraphics(gcf, [subpath 'bumpposVheading_cwoo_allfliesw5trials', add2txt, '.png'])
end
fig_no = fig_no + 1;

% c = [0.21, 0.01, .42, .39, .7, .05, .143, .22, .345, .33, .54, .2, .188, .34];
figure(5); clf;
h = histogram(c, [-1:.1:1]);
set(h, 'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', [0.5 0.5 0.5])
xlim([-1 1])

figure(6);
scatter(ones(1, length(c))+randn(1, length(c))*.01, c, 50, 'filled', 'k')
xlim([0 2])
ylim([-1 1])
[~, p] = ttest(c);
title(['p=' num2str(p, '%.3f')])


if print_figs
    print_fig(1, gcf, [subpath 'posVhead_corrdist_cwoo_allfliesw5trials', add2txt], 'eps', 1); 
    exportgraphics(gcf, [subpath 'posVhead_corrdist_cwoo_allfliesw5trials', add2txt, '.png'])
end
fig_no = fig_no + 1;


figure(8); clf; 
h = histogram(all_cc, 21); %, 'Normalization', 'probability')
[~, p] = ttest(all_cc);
set(h, 'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', [0.5 0.5 0.5])
title(['N=' num2str(length(all_cc)), ', p=' num2str(p)])
xlabel('Correlation Coefficient')
ylabel('Trials')
% ylabel('Probability')

figure(9); clf; hold on;
scatter(ones(1,length(allfly_cc_mns))+randn(1, length(allfly_cc_mns))*.05, allfly_cc_mns, 80, 'filled', 'k')
boxplot(allfly_cc_mns)
[~, p] = ttest(allfly_cc_mns);
plot([0 2],[0 0], '--k')
ylim([-1 1])
xlim([0 2])
title(['N=' num2str(length(allfly_cc_mns)), ', p=' num2str(p)])
ylabel(['corr coef'])