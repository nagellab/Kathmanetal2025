%% all fly
%%
clear all
close all


% use_flies = [21:27, 30, 32:34, 35, 36:39]; sdm = .5; line_label = 'hAck';
use_flies = [21:26, 30, 32:35, 37:39]; sdm = 0.5; line_label = 'hAck';
% use_flies = 35;
% use_flies = [25]; sdm = .5; line_label = 'hAck';
% use_flies = [41:51]; line_label = 'FC1';

print_figs = 0;
add2txt = '_stimfill_newwin';
plot_trials = 0;

fig_no = 1;

use_mets = [1, 2, 3, 4];

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

ylm(1,:) = [0 .9];
ylm(2,:) = [-.1 2.2];
ylm(3,:) = [-2 12];
ylm(4,:) = [0 110];
ylm(5,:) = [0 14];
ylm(6,:) = [0.75 1];
ylm(7,:) = [0 0];
ylm(8,:) = [0 0];
ylm(9,:) = [-100 2500];
ylm(10,:) = [-10 190];
ylm(11,:) = [-10 150];
%% metric windows
mwin = [-4 -1 5 8;
        -5 -2 0 3;
        -7 -4 -1.5 1.5;
        -5 -2 0 3;];

%% collect data
evpair = [];
figure(1); clf;
set(gcf, 'position', [1000, 0, 275, 110*length(use_mets)+200])
nm = length(use_mets); %number of metrics to plot
i = 1; N = 0; n = 0;
allfly_mnons = cell(nm,1);
comps = cell(nm, 1);
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;

    allamp_thr = [];
    for i = 1:length(data)
        allamp(i) = quantile(data(i).bump_amp(~isnan(data(i).bump_amp(:,1))), .95); 
        allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
    end
%     ylm = [-.05 nanmean(allamp)*1.5];
%     sdm = 1;
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;

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
        if data(t).closed_loop && data(t).odor_on && data(t).wind_on && ~any(~isnan(data(t).plume_t)) && ...
                ~isempty(data(t).bump_amp) && length(odoroni) == 1 && max(odori) > 0 && odori(end) == 0%cwoo        

            odoroffi = idy(odori(idy+1) < 0.5); %all pts where next pt is subthresh (aka downward crossing)
%             if length(odoroffi) == 1 
            if length(odoroni) == 1 %&& all([use_flies(f) ~= 35, t ~= 6]) %kill this trial (f35t6) bc spinning whole trial (make spinning trial filter)
                use_trials = [use_trials, t];
%                 allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
%                 allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
            end
        end
    end
%     base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;


%     allamp_thr = [];
%     for i = 1:length(data)
%         allamp(i) = quantile(data(i).bump_amp(~isnan(data(i).bump_amp(:,1))), .95); 
%         allamp_thr = [allamp_thr; data(i).bump_amp(~isnan(data(i).bump_amp(:,1)))];
%     end
% %     ylm = [-.05 nanmean(allamp)*1.5];
% %     sdm = 1;
%     base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;

    
    ons = cell(nm,1);
    trackm = [];
    bmp = [];
    for t = use_trials 

        nonan = ~isnan(data(t).bump_amp(:,1));
        odor = data(t).odor(nonan);
        odorevent = find(odor); 
        odorevent = [odorevent(1)+10 odorevent(end)+10];


        dff= data(t).bump_amp(nonan,1);
%         if data(t).closed_loop
%             heading = data(t).calc_windpos(nonan);
%         else
%             heading = data(t).calc_heading(nonan);
%         end
%         heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 100)));
%         speed = data(t).calc_speed(nonan);
%         fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
%         x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
%         y = cumsum(speed.*cos(heading*pi/180))/fps;
%         tt = data(t).calc_ts(nonan);
%         strait = data(t).straight(nonan);
%         avel = data(t).calc_deltaz(nonan);
%         avel = smoothdata(avel, 'movmean', 250);
%         fvel = abs(data(t).calc_deltapitch(nonan))*9.52/2;
%         fvel = smoothdata(fvel, 'movmean', 30); 
%         uwv = cos(heading*pi/180).*speed;
%         uwv = smoothdata(uwv, 'movmean', 30);

        %settings from bumptrig.m
        if data(t).closed_loop
            heading = smoothdata(data(t).calc_windpos(nonan), 'movmean', 1);
        else
            heading = smoothdata(data(t).calc_heading(nonan), 'movmean', 1);
        end
        heading = wrapTo180(rad2deg(smoothdata(unwrap(deg2rad(heading)), 'movmean', 10)));
        speed = data(t).calc_speed(nonan);
        fps = cat(1,data.fps);fps = median(fps(~isnan(fps)));
        x = -cumsum(speed.*sin(heading*pi/180))/fps;%+90;
        y = cumsum(speed.*cos(heading*pi/180))/fps;
        tt = data(t).calc_ts(nonan);
        avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 100);
        uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
        fvel = data(t).calc_deltapitch(nonan)*9.52/2;
        fvel = smoothdata(fvel, 'movmean', 30);   


        pmove = zeros(length(fvel), 1);
        pmove(fvel > 2) = 1;
        goal = mean(heading(find(odor)));
        goalv = cos((heading-goal)*pi/180).*speed;
        straight = calc_straightness(x, y, tt, 100);
        straight = smoothdata(straight, 'movmean', 10);
        vect_str = calc_vect_str(heading(find(odor)), tt(find(odor)), 150);
        cumsumdtheta_tail = nan(length(heading), 1);
%         cumsumdtheta(odorevent(2)+1:end) = cumsum(abs(diff(heading(odorevent(2):end))));
        cumsumdtheta_tail(odorevent(2):end) = cumsum(abs(avel(odorevent(2):end)));
        heading_tail = nan(length(heading), 1);
        heading_tail(odorevent(2):end) = abs(heading(odorevent(2):end));
        goalhead = abs(heading - mean(heading(odorevent(1):odorevent(2))));


        oni = []; offi = []; 
% %         [bump, oni, offi] = bumpfinder2(dff, tt, 5, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
%         [bump, oni, offi] = bumpfinder2(dff, tt, 3, 5, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
% %         [bump, oni, offi] = bumpfinder2(dff, tt, 5, 5, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration
        [bump, oni, offi] = bumpfinder2(dff, tt, 3, 3, base_thr);% maxoff = fill in gaps < "maxoff" sec wide, minon = gets rid of bumps less than "minon" duration


        win = [-10, 10]; 
        mets = [];
        for m = 1:nm
            mets(:,m) = eval(all_mets{use_mets(m)});
            met = []; 
            met = mets(:,m); 
            win1 = []; wint1= []; over = 0;
            good_c = [];

            for c = 1:length(oni) %offs
                if oni(c) > odorevent(1)-0% && offi(c) > odorevent(2)
    %                     if odor(oni(c))% ~odor(offi(c)) && offi(c) > odorevent(2) %~odor(offi(c))
                    win1 = (oni(c)+win(1)*fps):(oni(c)+win(2)*fps);
                    win1 = win1(win1 < length(tt) & win1 > 0);
                    wint1 = tt(win1)-tt(oni(c));
    %                         if ~isempty(offs{m}) && length(win1) ~= diff(win)*fps+1
                    if length(win1) ~= diff(win)*fps+1
                        this_met = []; this_on = [];
                        this_met = met(win1);
    %                             this_off = nan(size(offs{m},1),1);
                        this_on = nan(diff(win)*fps+1, 1);
                        if over
                            this_on(end-length(this_met)+1:end) = this_met;
                        else
                            this_on(1:length(this_met)) = this_met;
                        end
                        ons{m} = [ons{m}, this_on];
                    else
                        plot_wint = wint1;
                        ons{m} = [ons{m}, met(win1)];
                    end
%                     compwin1 = 
%                     comp = [comp, [mean(offs)]]
                    if plot_trials
                        subplot(size(use_mets,2),1,m); hold on; 
                        plot(wint1, met(win1), 'k')  
                    end

%                     figure(100); hold on; 
%                     plot(tt-tt(offi(c)), cumsumdtheta); 
%                     plot([tt(odorevent(2))-tt(offi(c)) tt(odorevent(2))-tt(offi(c))], ylim, 'k'); 
%                     plot([tt(offi(c))-tt(offi(c)) tt(offi(c))-tt(offi(c))], ylim, 'k') 

%                 taili = odorevent(2):offi(i);
%                 posti = offi(i)+1:offi(i)+length(taili);
%                 figure; hold on; 
%                 plot(heading); 
%                 plot([oni oni], ylim, 'r'); 
%                 plot([offi offi], ylim, 'r'); 
%                 plot([odorevent(1) odorevent(1)], ylim, 'k'); 
%                 plot([odorevent(2) odorevent(2)], ylim, 'k');
%                 plot(taili, ones(length(taili),1), 'r', 'LineW', 4)
%                 plot(posti, ones(length(taili),1), 'k', 'LineW', 4)

                    if m == 1
                        n = n + 1;
                    end
                end
            end
        end
    end

    for m = 1:nm
        if ~isempty(ons{m})
            allfly_mnons{m} = [allfly_mnons{m}, nanmean(ons{m},2)];
            subplot(size(use_mets,2),1,m); hold on; 
%             subplot(size(use_mets,2)+1,2,m*2+2); hold on; 
            mw = [];
            mw = nanmean(ons{m}, 2);
            plot(plot_wint, mw, 'k', 'LineW', 1)
%             plot([0 0], ylim, 'k', 'LineW', 0.5)
%             xlim([win(1)-1 win(2)+1]); %ylim(ylm(m,:))
%             ylabel(ylab{use_mets(m)})

            comps{m} = [comps{m}; [mean(mw(plot_wint > mwin(m,1) & plot_wint < mwin(m,2))), mean(mw(plot_wint > mwin(m,3) & plot_wint < mwin(m,4)))]];

            if m == 1
                N = N + 1;
            end
        end
    end


end

for m = 1:nm
    if ~isempty(ons{m})
        subplot(size(use_mets,2),1,m); hold on;
        plot(plot_wint, nanmean(allfly_mnons{m},2), 'r', 'LineW', 3)
%         for ii = 1:4
%             plot([mwin(m,ii) mwin(m,ii)], ylim, 'b')
%         end
        yl = ylim;
        for ii = [1,3]
            g = fill([mwin(m,ii) mwin(m,ii+1), mwin(m,ii+1) mwin(m,ii)], [yl(2) yl(2) yl(1) yl(1)], 'b');
            set( g, 'edgecolor', 'none', 'facealpha', .1 );
        end

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
        ylim(ylm(use_mets(m),:))
        ylabel(ylab{use_mets(m)})
        plot([0 0], ylim, 'k')
    end
end
% sgtitle(['fly#' num2str(use_flies(f)) ', N = ' num2str(N) ', n(bumps)=' num2str(n)])
sgtitle(['N = ' num2str(N) ', n(bumps)=' num2str(n)])

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
    mncmp(2) = nanmean(comps{m}(:,2));
    plot([.75 1.25], [mncmp(1) mncmp(1)], 'r', 'LineW', 3)
    plot([1.75 2.25], [mncmp(2) mncmp(2)], 'r', 'LineW', 3)
    [~, p] = ttest2(comps{m}(:,1), comps{m}(:,2));
    sigstar({[1,2]}, p)
    title(['m=[', num2str(mncmp(1), '%.2f'), ',', num2str(mncmp(2),'%.2f'), '], p=' num2str(p, '%.4f')])
    xlim([0 3])
    ylabel(ylab{use_mets(m)})
end

figure(1);
if print_figs
    print_fig(1, gcf, [subpath 'bumpontrig_' line_label add2txt], 'eps', 1); 
    exportgraphics(gcf, [subpath 'bumpontrig_' line_label add2txt, '.png'])

%     print_fig(1, gcf, [subpath 'fly' num2str(use_flies(f)) '_bumptrig_cumsumdthetaATodoroff' line_label add2txt], 'eps', 1); 
%     exportgraphics(gcf, [subpath 'fly' num2str(use_flies(f)) '_bumptrig_cumsumdthetaATodoroff' line_label add2txt, '.png'])
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