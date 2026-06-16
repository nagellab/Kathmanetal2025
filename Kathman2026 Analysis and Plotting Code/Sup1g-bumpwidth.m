%% all fly
%%
clear all
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/2024 paper/sandbox/exploratory figs/bump postion/';
cd(headpath)

allflies_filenames

%%
% use_flies = [21:26,30, 33:35, 37:39]; sdm = 1;
use_flies = [21:27, 30, 32:39]; sdm = .25;
% use_flies = [41:51]; sdm = .5;

print_figs = 0;

fig_no = 1;
figure(fig_no); clf; hold on
set(gcf, 'Position', [1 1 620 550]); %[100 100 500 1000])

j = 1; jj = 1;
flies_b = [];
flies_bwo = [];
ntrials = [];
for f = 1:length(use_flies)
    load(filename{use_flies(f)})
    data = syncstruct;



    %find use_trials and base_thresh for whole expt
    allamp_thr = []; use_trials = [];
    for t = 1:length(data)
        if ~isempty(data(t).bump_amp)
            nonan = ~isnan(data(t).bump_amp(:,1));
    
            %find stim times before trial sorting bc used to rule out plume trials
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)

            wind = data(t).wind(nonan);
            idyl = wind >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            windoni = idy(wind(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
    
            allamp(t) = quantile(data(t).bump_amp(~isnan(data(t).bump_amp(:,1))), .95); 
            allamp_thr = [allamp_thr; data(t).bump_amp(~isnan(data(t).bump_amp(:,1)))];
            if length(odoroni) == 1  && length(windoni) == 1 %&& ~data(t).closed_loop
                use_trials = [use_trials, t];
            end
        end
    end
    base_thr = mean(allamp_thr)+std(allamp_thr)*sdm;


    trials_b = [];
    trials_bwo = [];
    if ~isempty(use_trials)
        for t = use_trials
            nonan = ~isnan(data(t).bump_amp(:,1));

            %find stim times (again, dumb)
            odor = data(t).odor(nonan);
            idyl = odor >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            odoroni = idy(odor(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            odoroni = odoroni + 10;
            odoroffi = idy(odor(idy(idy < length(odor))+1) < .5); %all pts where next pt is subthresh (aka downward crossing)
            odoroffii = odoroffi + 10; %compesate for odor delay (100ms)

            wind = data(t).wind(nonan);
            idyl = wind >= 0.5; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            windoni = idy(wind(idy-1) < 0.5); %all pts where prev pt is subthresh (aka crossing)
            windoni = windoni + 10;
            windoffi = idy(wind(idy(idy < length(wind))+1) < .5); %all pts where next pt is subthresh (aka downward crossing)
            windoffii = windoffi + 10; %compesate for odor delay (100ms)            

            %find dff and behav metrics
            dff = data(t).bump_amp(nonan, 1);
            speed = data(t).calc_speed(nonan);
            fps = cat(1,data.fps);fps = median(fps(~isnan(fps))); %fictrac fs
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
            avel = smoothdata(data(t).calc_deltaz(nonan), 'movmean', 200);
            abs_avel = smoothdata(abs(data(t).calc_deltaz(nonan)), 'movmean', 50);
            uwv = smoothdata(cos(heading*pi/180).*speed, 'movmean', 30);%60
            fvel = data(t).calc_deltapitch(nonan)*9.52/2;
            fvel = smoothdata(fvel, 'movmean', 200); 


            im = data(t).dff_green(1:8,:);
            ts = data(t).dff_ts;
            fs = size(im, 2)/max(ts); %thor fs           
            [pos, uwpos, ~, thr] = bumpposfinder(im, ts, tt, base_thr, sdm); %thr is min of expt thresh (base_thr) and trial thresh

%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 3, thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
%             [bump, oni, offi, thrsh] = bumpfinder2(dff, tt, 5, 3, base_thr); %unsmoothed dff(nonan), tt(nonan), maxoff, minon
            [bump, oni, offi, thrsh] = bumpfinder2(max(im), ts, 5, 3, thr); 
            bumpi = find(bump);

            %find bump widths
            bwo = []; bw = []; b = []; %wind+odor, wind, none
            for i = 1:length(bumpi)
                this_bump = im(:,bumpi(i));
                [~, pkx] = max(this_bump);
                centi = mod((0:7)+(4-pkx), 8) + 1; %center bump peak to column 4
                this_cbump = nan(1,8);
                this_cbump(centi) = this_bump;
%                 figure(1); clf; hold on
%                 plot(this_bump, 'k')
%                 plot(cent_bump, 'r')

                if ts(bumpi(i)) > tt(odoroni)+2 && ts(bumpi(i)) < tt(odoroffi)
%                     bwo = [bwo; this_cbump];
                    bwo = [bwo; this_cbump-min(this_cbump)]; %min dff set to zero
                elseif ts(bumpi(i)) < tt(windoni) || ts(bumpi(i)) > tt(windoffi)
%                     b = [b; this_cbump];
                    b = [b; this_cbump-min(this_cbump)]; %min dff set to zero
                end
            end

            if size(b,1) > 2 
                trials_b = [trials_b; nanmean(b,1)];
                j = j + 1;
            end
            if size(bwo,1) > 2 
                trials_bwo = [trials_bwo; nanmean(bwo,1)];
                jj = jj + 1;
            end
%                 if f == 11
%                     keyboard; end
        end
    end
    if size(trials_b,1) > 1
        flies_b = [flies_b; nanmean(trials_b,1)];
        plot([-3:4], nanmean(trials_b,1), 'k')
%         size(trials_b,1)
    end
    if size(trials_bwo,1) > 1 
        flies_bwo = [flies_bwo; nanmean(trials_bwo,1)];
        plot([-3:4], nanmean(trials_bwo,1), 'm') 
        ntrials = [ntrials, size(trials_bwo,1)];
    end

end

plot([-3:4], nanmean(flies_b,1), 'k', 'LineW', 6)

plot([-3:4], nanmean(flies_bwo,1), 'm', 'LineW', 6)   





xlabel('columns')
ylabel('dff')
title(['magenta=during odor, black=no stim (N=' num2str(size(flies_bwo,1)), ', trials/fly=', num2str(ntrials), ', n=', num2str(jj), ', n=' num2str(j), ')'])




