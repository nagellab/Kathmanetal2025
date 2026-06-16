%% all fly
%%
clear all
close all


use_flies = [27, 32, 36, 38, 109:112];
sdm = .5; line_label = 'hAck';
% use_flies = 35;
% use_flies = [25]; sdm = .5; line_label = 'hAck';
% use_flies = [41:51]; line_label = 'FC1';

print_figs = 1;
add2txt = '_windowcomp';
add2txt = '_bmplngth_plume3_longavelwin';
bmp1 = 3;
bmp2 = 3;
plot_trials = 0;

fig_no = 1;

use_mets = [1, 2, 3, 11];

all_mets{1} = 'dff'; 
all_mets{2} = 'abs(avel)';%'pos';
all_mets{3} = 'fvel';
all_mets{4} = 'abs(heading)';%'abs(heading)';
all_mets{5} = 'uwv';
all_mets{6} = 'straight'; %'straight'; 'vect_str';
% all_mets{7} = 'pos';
all_mets{7} = 'odor';
all_mets{8} = 'wind';
all_mets{9} = 'cumsumdtheta_tail';
all_mets{10} = 'heading_tail';
all_mets{11} = 'goalhead';



%% load data and start figs
% headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
headpath = ['/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames
subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/exploratory figs/bump trig/';
subpath = '/Users/kathmn01/Desktop/';
subpath = '/Users/kathmn01/Desktop/revision figs/';

ylab = cell(1);
ylab{1} = 'dff';
ylab{2} = 'abs(avel)(deg/s)';%'bump pos'; %abs(avel)(deg/s)';
ylab{3} = 'fvel(mm/s)';
ylab{4} = 'abs(heading)(deg)';%'abs(heading)(dg)';
ylab{5} = 'uwv(deg/s)';
ylab{6} = 'straight';
% ylab{7} = 'pos (col)'; %'straight'; 'vect str'; [.65 1]
ylab{7} = 'odor';
ylab{8} = 'wind';
% ylab{9} = 'odor_heading';
ylab{9} = 'cumsum delta theta(tailonly)';
ylab{10} = 'abs(heading)-tailonly';
ylab{11} = 'abs(head-goal)';

ylm(1,:) = [0 .7];
ylm(2,:) = [-.1 1.5];
ylm(3,:) = [1 5];
ylm(4,:) = [-10 190];
ylm(5,:) = [2 8];
ylm(6,:) = [0 .2];
ylm(7,:) = [0 0];
ylm(8,:) = [0 0];
ylm(9,:) = [-100 2500];
ylm(10,:) = [-10 190];
ylm(11,:) = [-10 150];
%% metric windows
mwins = [0 3; -1 2; -1 2; 0 3];
mwins = [0 3; 0 1; -1 2; 0 5];
mwins = [0 3; 0 2.5; -1 2; 0 5];

mwin = [-6 -3 3 6;
        -6 -3 0 3;
        -6 -3 -1.5 1.5;
        -6 -3 0 3;];
% 
% mwin = [-6 -3 3 6;
%         -6 -3 0 3;
%         -6 -3 -1.5 1.5;
%         -6 -3 0 3;];

% mwin = [-4 -1 5 8;
%         -5 -2 3 6;
%         -4 -1 5 8;
%         -4 -1 5 8;];

%% collect data
evpair = [];
figure(1); clf;
set(gcf, 'position', [1000, 0, 275, 110*length(use_mets)+200])
nm = length(use_mets); %number of metrics to plot
i = 1; N = 0; n = 0;
allfly_mnoffs = cell(nm,1);
allfly_mnons = cell(nm,1);
comps = cell(nm, 1);

A_trials = cell(nm,1);
B_trials = cell(nm,1);

for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;

    %find use_trials and alltrial_dff data to set bump threshold
    use_trials = []; allamp_thr = [];
    for t = 1:length(data)
        nonan = ~isnan(data(t).bump_amp(:,1));

        %find stim times before trial sorting bc used to rule out plume trials
        odori = data(t).odor(nonan);
        idyl = odori >= 0.5; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        odoroni = idy(odori(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
%         if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ...
%                 ~isempty(data(t).bump_amp) && length(odoroni) == 1 && max(odori) > 0 && odori(end) == 0%cwoo
        if data(t).closed_loop && data(t).odor_on && data(t).wind_on && any(~isnan(data(t).plume_t)) && ...
                ~isempty(data(t).bump_amp) && length(odoroni) > 1 && max(odori) > 0 && odori(end) == 0%plume
                use_trials = [use_trials, t];
%                 allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
                allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
        end
    end
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;


%     allamp_thr = [];
%     for i = 1:length(data)
%         allamp(i) = quantile(data(i).bump_amp(~isnan(data(i).bump_amp(:,1))), .95); 
%         allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
%     end
% %     ylm = [-.05 nanmean(allamp)*1.5];
% %     sdm = 1;
%     base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;


    offs = cell(nm,1);
    ons = cell(nm,1);
    trackm = [];
    bmp = [];
    if ~isempty(use_trials)
    for t = use_trials
        nonan = ~isnan(data(t).bump_amp(:,1));

        tt = data(t).calc_ts(nonan);
        odor = data(t).odor(nonan);
        h = exp(-tt/2); %exponential function, tau = 2s;
        ofreq = conv(h,odor)/10;
        ofreq = ofreq(1:length(odor));


        odorevent = find(odor); 
        odorevent = [odorevent(1)+10 odorevent(end)+10];

        wind = data(t).wind(nonan);
        ampthresh = .5;%0.125;
        idyl = odor >= ampthresh; %all superthresh pts
        idy = find(idyl); %superthresh indices
        idy = idy(idy>1); %all subthresh after 1
        oncrossi = idy(odor(idy-1)<ampthresh);
        offcrossi = idy(odor(idy(idy < length(odor))+1)<ampthresh); %all pts where next pt is subthresh (aka downward crossing)
        offcrossi = offcrossi + 10; %compesate for odor delay (100ms)

        dff = data(t).bump_amp(nonan, 1);
        speed = data(t).calc_speed(nonan);
        fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
        if data(t).closed_loop
            heading = data(t).calc_windpos(nonan);
        else
            heading = data(t).calc_heading(nonan);
        end
        heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
        uwheading = smoothdata(unwrap(deg2rad(heading)), 'movmean', 1);

        col = data(t).bump_col(nonan, 1);
        com = data(t).bump_com(nonan, 1);
%                 uwpos = unwrap(data(t).bump_col(nonan);
        oldpos = smoothdata(data(t).bump_col(nonan), 'movmean', 100); %200);
%                 newpos = data(t).bump_uwspos(nonan);
%                 pos = data(t).bump_col(nonan);
        avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 200);
        abs_avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 50);
        uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
        fvel = data(t).calc_deltapitch(nonan)*9.52/2;
        fvel = smoothdata(fvel, 'movmean', 200);

        x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
        y = cumsum(speed.*cos(heading*pi/180))/fps;
        pmove = zeros(length(fvel), 1);
        pmove(fvel > 2) = 1;
%         strait = data(t).straight(nonan);        
%         straight = calc_straightness(x, y, tt, 100);
%         straight = smoothdata(straight, 'movmean', 10);     


        oni = []; offi = []; 
%         [bump, oni, offi] = bumpfinder2(dff, tt, 5, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
%         [bump, oni, offi] = bumpfinder2(dff, tt, 3, 5, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
% % % % %         [bump, oni, offi] = bumpfinder2(dff, tt, 5, 5, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
        [bump, oni, offi] = bumpfinder2(dff, tt, bmp1, bmp2, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration


        win = [-10, 10]; 
        for c = 1:length(offi) %offs
            bumpn = offi(c)-oni(c);
            if max(odor(oni(c):offi(c))) > 0 && any(offcrossi > oni(c) & offcrossi < offi(c))% && max(ofreq(oni(c):offi(c)) > 0.4) && bumpn < 6000 %longer bumps are prob errors

                these_offcrosses = offcrossi(offcrossi > oni(c) & offcrossi < offi(c));
                this_odor = [oni(c):these_offcrosses(end)];

%                     goal = mean(heading(find(this_odor)));
%                     goalv = cos((heading-goal)*pi/180).*speed;
        
                vect_str = calc_vect_str(heading(this_odor), tt(this_odor), 150);
                cumsumdtheta_tail = nan(length(heading), 1);
                cumsumdtheta_tail(this_odor(end):end) = cumsum(abs(avel(this_odor(end):end)));
                heading_tail = nan(length(uwheading), 1);
%                 heading_tail(this_odor(end):end) = abs(uwheading(this_odor(end):end)-mean(uwheading(this_odor)));
                heading_tail(this_odor(end):end) = abs(heading(this_odor(end):end));
                goaldhead = [];
                goalhead = abs(uwheading - mean(uwheading(this_odor)));

                for m = 1:nm

                    met = []; win1 = []; wint1= []; 
                    met = eval(all_mets{use_mets(m)}); 
                    win1 = (offi(c)+win(1)*fps):(offi(c)+win(2)*fps); %window indices around bump off index
                    win1 = win1(win1 < length(tt) & win1 > 0); %truncate if extend beyond trial
                    wint1 = tt(win1)-tt(offi(c)); %show in trial time too

                    if length(win1) ~= diff(win)*fps+1 %if truncated, add nans
                        this_met = []; this_off = [];
                        this_met = met(win1);

                        this_off = nan(diff(win)*fps+1, 1); 

                        this_off(1:length(this_met)) = this_met;
                    else
                        offt = wint1;
                        this_off = met(win1);
                    end
                    offs{m} = [offs{m}, this_off]; %all centered traces around bump offset

                    if bumpn > 9000
                        trim4comp = bumpn - 800;
                    else
                        trim4comp = 1000;
                    end
                    trim4comp = 500;
                    this_on = nan(60*fps+1, 1);
                    starti = ((60*fps)/2)-round((bumpn-trim4comp)/2);
                    this_on(starti:starti+(bumpn-trim4comp)) = met((oni(c)+trim4comp):offi(c));
                    ons{m} = [ons{m}, this_on]; %all centered traces during bump on period
                    ont = -(60/2):(1/fps):(60/2);

                    if plot_trials
                        figure(1);
                        subplot(size(use_mets,2),1,m); hold on; 
                        plot(wint1, met(win1), 'k') 
                        xlim([-12 12])
                    end

                    if m == 1
                        n = n + 1;
                    end
                end
            end
        end
    end
    end

    %plot trial mean for this fly for each met
    for m = 1:nm
        if ~isempty(offs{m})

            these_ons = ons{m}(2001:end-2000,:);
            ont = -10:.01:10;


            A_trials{m} = [A_trials{m}; offs{m}'];
            B_trials{m} = [B_trials{m}; these_ons'];            

            mw1 = nanmean(offs{m}, 2);
            mw2 = nanmean(these_ons, 2);

            allfly_mnoffs{m} = [allfly_mnoffs{m}, mw1];
            allfly_mnons{m} = [allfly_mnons{m}, mw2];


% % %             subplot(size(use_mets,2),1,m); hold on; 
% % % %             subplot(size(use_mets,2)+1,2,m*2+2); hold on; 
% % %             plot(plot_wint, mw, 'k', 'LineW', 1)
% % %             if plot_trials
% % %                 plot(plot_wint, mw, 'r', 'LineW', 2)
% % %             end
% % % %             plot([0 0], ylim, 'k', 'LineW', 0.5)
% % % %             xlim([win(1)-1 win(2)+1]); %ylim(ylm(m,:))
% % % %             ylabel(ylab{use_mets(m)})
% %             
% %             %for stat comparisons
% % %             comps{m} = [comps{m}; [mean(mw1(offt > mwin(m,1) & offt < mwin(m,2))), mean(mw1(offt > mwin(m,3) & offt < mwin(m,4)))]];

%             %find difference of beginning 5s and ending 5s for this fly (for off condition and on condition)
%             startmn = nanmean(mw1(1:500,:));
%             endmn = nanmean(mw1(end-1000+1:end-500,:));
%             diff1 = endmn-startmn;
% 
%             startmn = nanmean(mw2(1:500,:));
%             endmn = nanmean(mw2(end-1000+1:end-500,:));
%             diff2 = endmn-startmn;
% 
%             comps{m} = [comps{m}; [diff1, diff2]];

            %find means of a time window for both offs and ons
%             mwins = [0 3; -1 2; -1 2; 0 3];
% %             mwins = [1 4; 1 4; 1 4; 1 4];
            win_off = nanmean(mw1(offt > mwins(m,1) & offt < mwins(m,2)));
            win_on = nanmean(mw2(offt > mwins(m,1) & offt < mwins(m,2))); 

            comps{m} = [comps{m}; [win_off, win_on]];

% if f == 8 && m == 1
%     keyboard
% end

            if m == 1
                N = N + 1;
            end
        end
    end


end

%% plot fly mean for each met
for m = 1:nm
    if ~isempty(allfly_mnoffs{m})
        d1 = allfly_mnoffs{m}';
        mn1 = nanmean(d1,1);
        se1 = nanstd(d1,0,1)./sqrt(sum(~isnan(d1),1));
        ts = offt';

        subplot(size(use_mets,2),1,m); hold on;
        plot(ts, mn1, 'r', 'LineW', 3)
        fill([ts, fliplr(ts)], [mn1 + se1, fliplr(mn1 - se1)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none')

        d2 = allfly_mnons{m}';
        mn2 = nanmean(d2,1);
        se2 = nanstd(d2,0,1)./sqrt(sum(~isnan(d2),1));
        ts = ont;

        subplot(size(use_mets,2),1,m); hold on;
        plot(ts, mn2, 'k', 'LineW', 3)
        fill([ts, fliplr(ts)], [mn2 + se2, fliplr(mn2 - se2)], 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none')        


%         for ii = 1:4
%             plot([mwin(m,ii) mwin(m,ii)], ylim, 'b')
%         end
        yl = ylim;
%         for ii = [1,3]
%             g = fill([mwin(m,ii) mwin(m,ii+1), mwin(m,ii+1) mwin(m,ii)], [yl(2) yl(2) yl(1) yl(1)], 'b');
%             set( g, 'edgecolor', 'none', 'facealpha', .1 );
%         end

            g = fill([mwins(m,1) mwins(m,2), mwins(m,2) mwins(m,1)], [yl(2) yl(2) yl(1) yl(1)], 'b');
            set( g, 'edgecolor', 'none', 'facealpha', .1 );


%             ts = ts(~isnan(ts));
%             binarized = binarized(~isnan(binarized));
%             binarized = make_row(binarized);
%             on = find(diff(binarized)>0);
%             off = find(diff(binarized)<0);
%             if binarized(1) 
%                 on = [1 on]; end
%             if binarized(end)
%                 off = [off length(binarized)]; end
%             for s = 1:length(on)
%                 g = fill(ax_handle, [ts(on(s)), ts(off(s)), ts(off(s)), ts(on(s))], [ylim(2) ylim(2) ylim(1) ylim(1)], colorstr);
%                 set( g, 'edgecolor', 'none', 'facealpha', facealpha );
%             end

%         plot([-2.5 -2.5], ylim, 'k')
        xlim([win(1)-1 win(2)+1]); 
%         ylim(ylm(use_mets(m),:))
        ylabel(ylab{use_mets(m)})
        plot([0 0], ylim, 'k')



        %sliding window sig diff test
        % Parameters
        fs = 100;                     % Sampling rate (Hz)
%         win_size_sec = .15;           % Sliding window size (in seconds)
%         step_size_sec = 0.01;          % Step size (in seconds)
        win_size_sec = 0.5;           % Sliding window size (in seconds)
        step_size_sec = 0.1;          % Step size (in seconds)
        
        win_size = round(win_size_sec * fs);  % Samples per window
        step_size = round(step_size_sec * fs);
        n_timepoints = size(A_trials{m}, 2);
        
        % Preallocate
        n_windows = floor((n_timepoints - win_size) / step_size) + 1;
        p_values = NaN(1, n_windows);
        h_values = NaN(1, n_windows);
        t_values = NaN(1, n_windows);
        time_centers = NaN(1, n_windows);
        
        % Sliding window t-test
        for i = 1:n_windows
            start_idx = (i-1)*step_size + 1;
            end_idx = start_idx + win_size - 1;
            
            A_win = nanmean(A_trials{m}(:, start_idx:end_idx), 2);
            B_win = nanmean(B_trials{m}(:, start_idx:end_idx), 2);
            
%             [h_ttest, p, ~, stats] = ttest2(A_win, B_win);
            [h_ttest, p, ~, stats] = ttest2(A_win, B_win, 'alpha', 0.03);
            if p < 0.05
                h_test = 1; end
            h_values(i) = h_ttest;
            p_values(i) = p;
            t_values(i) = stats.tstat;
            time_centers(i) = (start_idx + end_idx) / 2 / fs;  % in seconds
        end
        
%       correct for multiple comparisons (e.g., FDR)
        [h_fdr, crit_p, adj_p] = fdr_bh(p_values);  % Benjamini-Hochberg FDR

        % Plot sig windows
%         time_to_shade = time_centers(find(h_fdr)) - 10;
        time_to_shade = time_centers(find(h_values)) - 10;
        yl = ylim;
        
        % Group nearby times into contiguous segments (define a gap threshold)
        dt = diff(time_to_shade);
        gap_thresh = 0.15;  % seconds between groups
        breaks = [0, find(dt > gap_thresh), length(time_to_shade)];
        segments = arrayfun(@(i) time_to_shade(breaks(i)+1:breaks(i+1)), 1:length(breaks)-1, 'UniformOutput', false);
        
        % Plot the vertical shading (or line) for each segment
        if ~isempty(segments{1})
            for i = 1:length(segments)
                seg = segments{i};
                x1 = seg(1);
                x2 = seg(end) + 0.01;  % small width beyond the last point
                plot([x1 x2], [yl(2) yl(2)], 'k', 'LineW', 2)
%                 fill([x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
%                      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            end
        end
        ylim([yl(1) yl(2)*1.1])




    end
end
% sgtitle(['fly#' num2str(use_flies(f)) ', N = ' num2str(N) ', n(bumps)=' num2str(n)])
sgtitle(['N = ' num2str(N) ', n(bumps)=' num2str(n)])


%% plot bump/nobump means as paired scatters
fig_no = fig_no + 1;
figure(fig_no); clf;
set(gcf, 'position', [100, 0, 210*4+200, 400])
for m = 1:4
    subplot(1,4,m); hold on
    for mm = 1:size(comps{m},1)
        plot([1 2], comps{m}(mm,:), 'k')
        plot([1 2], comps{m}(mm,:), '.k')
    end
    mncmp(1) = nanmean(comps{m}(:,1));
    nmcmp(2) = nanmean(comps{m}(:,2));
    plot([.75 1.25], [mncmp(1) mncmp(1)], 'r', 'LineW', 3)
    plot([1.75 2.25], [nmcmp(2) nmcmp(2)], 'r', 'LineW', 3)
    [~, p] = ttest2(comps{m}(:,1), comps{m}(:,2));
    sigstar({[1,2]}, p)
    title(['m=[', num2str(nmcmp(1), '%.2f'), ',', num2str(nmcmp(2),'%.2f'), '], p=' num2str(p, '%.4f')])
    xlim([0 3])
    ylabel(ylab{use_mets(m)})
end

figure(1);
if print_figs
%     print_fig(1, gcf, [subpath 'bumpofftrigVon_plume_slidwin_' line_label add2txt], 'eps', 1); 
%     exportgraphics(gcf, [subpath 'bumpofftrigVon_plume_slidwin_' line_label add2txt, '.png'])
    savefig([subpath 'bumpofftrigVon_plume_' line_label add2txt])
end

figure(2);
if print_figs
%     print_fig(1, gcf, [subpath 'bumpofftrigVonStats_plume_slidwin_' line_label add2txt], 'eps', 1); 
%     exportgraphics(gcf, [subpath 'bumpofftrigVonStats_plume_slidwin_' line_label add2txt, '.png'])
    savefig([subpath 'bumpofftrigVonStats_plume_' line_label add2txt])
end



%% calc_straightness
function straight = calc_straightness(x, y, tt, win)
    pth = []; strait = [];
    for ii = 1:length(x)-1
        pth(ii) = sqrt((x(ii+1)-x(ii))^2 + (y(ii+1)-y(ii))^2);
    end
    pth = [pth, pth(end)];
    
    % pth = sqrt(x.^2 + y.^2);
    i = 1;
    for j = 1:length(x)-win
        strait(i) = sqrt((x(j+win)-x(j))^2 + (y(j+win)-y(j))^2)/sum(pth(j:j+win));
        strait(i) = strait(i)^2;
        i = i + 1;
    end
    straight = nan(length(tt),1);
    straiti = round(win/2):(round(win/2)+length(strait)-1);
    straight(straiti) = strait;
    
end

%% calc_vect_str
function vect_str = calc_vect_str(heading, tt, win)
    vector_strengths = [];
    vector_angles = [];
    bin_num = 150;

    i = 1;
    for j = 1:length(tt)-win
        [counts,centers] = hist(heading(j:j+win), win); 
        counts = counts / sum(~isnan(heading));
        centers = centers * pi/180;
    
        % Convert polar plot histogram into data x,y cartesian data points
        data_x = counts.* cos(centers);
        data_y = counts.* sin(centers);
    
        vec_x = sum(data_x);
        vec_y = sum(data_y);
    
        vec_strength = sqrt(vec_x^2 + vec_y^2);
    
        vector_strengths = [vector_strengths, vec_strength];
%                 vector_angle = atan(vec_y/vec_x); % Find angle
%     
%                 % Correct for angle if x component of PC vector is negative
%                 if vec_x < 0
%                     vector_angle = vector_angle + pi;
%                 end
%                 vector_angles = [vector_angles, vector_angle];
    end
    vect_str = nan(length(tt),1);
    vsi = round(win/2):(round(win/2)+length(vector_strengths)-1);
    vect_str(vsi) = vector_strengths;
end    