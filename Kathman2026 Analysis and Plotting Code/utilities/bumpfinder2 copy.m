function [bump, oni, offi, thresh] = bumpfinder2(dff, tt, maxoff, minon, base_thr)
% [bump, oni, offi, thresh] = bumpfinder2(dff, tt, maxoff, minon, base_thr)
% dff, tt needs to be nonaned
% maxoff = fill in gaps < "maxoff" sec wide
% minon = gets rid of bumps less than minon duration

            fs = 1/median(diff(tt));
%             ipt = []; bases = [];bumpi = [];thresh = [];
            sdff = smoothdata(dff, 'movmean', round(fs)); %1 sec
%             sdff = smoothdata(dff, 'movmean', 100); %hack

            %set some boundries for trials with no bump
            thresh = base_thr;
            if thresh < 0.2  
                thresh = 0.2; end

%             thresh = 0.23; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manual theshhold if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            bumpi = find(sdff > thresh);
%             bumpi = sdff > thresh(use_flies(f));
%             figure(2); clf; hold on; plot(tt, sdff); plot([0 65], [thresh, thresh], 'm')

            %fill in gaps < "maxoff" sec wide
            gaps = diff(bumpi);idx = find(gaps > 1 & gaps < maxoff*round(fs)); 
            bump = zeros(length(sdff), 1);
            bump(bumpi) = 1;
            for g = 1:length(idx)
                this_gap = [];
                this_gap = bumpi(idx(g))+1:bumpi(idx(g)+1)-1;
                bump(this_gap) = 1;
            end
            bumpi = find(bump);
%             plot(tt(bumpi), sdff(bumpi), '*g')

        if ~isempty(bumpi)
            ampthresh = .5;%0.125;
            idyl = bump >= ampthresh; %all superthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            oncrossi = idy(bump(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)


            ampthresh = 0.5; %set thresh
            idyl = []; idy = []; 
            idyl = bump <= ampthresh; %all subthresh pts
            idy = find(idyl); %superthresh indices
            idy = idy(idy>1); %all subthresh after 1
            offcrossi = idy(bump(idy-1)>ampthresh); %all pts where prev pt is superthresh (aka crossing)

            if isempty(oncrossi) & length(offcrossi) == 1
                oncrossi = 1; end
            if isempty(offcrossi) & length(oncrossi) == 1
                offcrossi = length(bump); end
            if isempty(oncrossi) & isempty(offcrossi)
                oncrossi = 1; offcrossi = length(bump); end

            if offcrossi(1) < oncrossi(1)
                oncrossi = [1;oncrossi]; end
            if oncrossi(end) > offcrossi(end)
                offcrossi = [offcrossi; length(bump)]; end

            keep = find((offcrossi - oncrossi) > minon*round(fs)); %gets rid of bumps less than minon duration
            bump = zeros(length(sdff), 1);
%             bump(bumpi) = 1;
            for k = 1:length(keep)
                bump(oncrossi(keep(k)):offcrossi(keep(k))) = 1;
            end
            bumpi = find(bump);
            oni = oncrossi(keep);
            offi = offcrossi(keep);
            ont = tt(oni);
            offt = tt(oni);
        else
            oni = [];
            offi = [];
        end
% keyboard

%% old
%             fs = 1/median(diff(tt));
% %             ipt = []; bases = [];bumpi = [];thresh = [];
%             sdff = smoothdata(dff, 'movmean', round(fs)); %1 sec
% %             sdff = smoothdata(dff, 'movmean', 100); %hack
% 
% % %             %find different Ca states for this trial
% % % %             ipt = findchangepts(sdff,'Statistic','std','MinThreshold',round(fs*20)); %hack
% % % %             ipt = findchangepts(sdff,'Statistic','std','MinThreshold',2000);%hack!!!!!!!!!!!!!!!!!!!
% % %             ipt = findchangepts(sdff); %fc1!!!!!!!!!  %%%% i think this is now obsolete without MinThresh bc it doesn't show change pts anymore? check, it should be!!!
% % % %             ipt = [1, ipt, length(sdff)];
% % %             ipt = [1; ipt; length(sdff)];
% % %             for i = 2:length(ipt)
% % %                 bases(i-1) = mean(sdff(ipt(i-1):ipt(i)));
% % %             end
% % % 
% % %             % find separation between Ca states (midpoint of kmeans) %%% also now obsolete bc bases == 1
% % %             if length(bases) > 1   
% % %                 [idx, c] = kmeans(bases', 2);
% % %                 thresh = min(c)+abs(diff(c))/2;
% % % %                 thresh = mean(bases);
% % %             else
% % %                 thresh = mean(sdff)+std(sdff)*2.5; %set some boundries for trials with no bump
% % %             end
% 
%             %set some boundries for trials with no bump
%             thresh = base_thr;
%             if thresh < 0.2  
%                 thresh = 0.2; end
% % %             if length(bases) == 2 && abs(diff(bases)) < std(sdff)*2.5% && thresh <= 0.2
% % % %                 thresh = 0.3;
% % % %                 thresh = mean(sdff)+ std(sdff)*2.5;
% % % %                 thresh = mean(sdff)+ std(sdff)*1.5;
% % %                 thresh = base_thr;
% % %             end
% %             thresh = 0.23; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% manual theshhold if needed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             bumpi = find(sdff > thresh);
% %             bumpi = sdff > thresh(use_flies(f));
% %             figure(2); clf; hold on; plot(tt, sdff); plot([0 65], [thresh, thresh], 'm')
% 
%             %fill in gaps < "maxoff" sec wide
%             gaps = diff(bumpi);idx = find(gaps > 1 & gaps < maxoff*round(fs)); 
%             bump = zeros(length(sdff), 1);
%             bump(bumpi) = 1;
%             for g = 1:length(idx)
%                 this_gap = [];
%                 this_gap = bumpi(idx(g))+1:bumpi(idx(g)+1)-1;
%                 bump(this_gap) = 1;
%             end
%             bumpi = find(bump);
% %             plot(tt(bumpi), sdff(bumpi), '*g')
% 
%         if ~isempty(bumpi)
%             ampthresh = .5;%0.125;
%             idyl = bump >= ampthresh; %all superthresh pts
%             idy = find(idyl); %superthresh indices
%             idy = idy(idy>1); %all subthresh after 1
%             oncrossi = idy(bump(idy-1)<ampthresh); %all pts where prev pt is subthresh (aka crossing)
% %             ontimes = tt(oncrossi);
% 
% %             for c = 1:length(oncrossi)
% %                 win1 = (oncrossi(c)-win(1)*fps):(oncrossi(c)+win(2)*fps);
% %                 wint1 = tt(win1)-tt(oncrossi(c));
% %                 ons = [ons, met(win1)];
% %                 subplot(2,2,1); hold on; plot(wint1, bumpi(win1), 'k')
% %                 subplot(2,2,3); hold on; plot(wint1, met(win1), 'k')
% %             end
% 
%             ampthresh = 0.5; %set thresh
%             idyl = []; idy = []; 
%             idyl = bump <= ampthresh; %all subthresh pts
%             idy = find(idyl); %superthresh indices
%             idy = idy(idy>1); %all subthresh after 1
%             offcrossi = idy(bump(idy-1)>ampthresh); %all pts where prev pt is superthresh (aka crossing)
% %             offtimes = tt(offcrossi);
% 
% %             for c = 1:length(offcrossi)
% %                 win2 = (offcrossi(c)-win(1)*fps):(offcrossi(c)+win(2)*fps);
% %                 wint2 = tt(win2)-tt(offcrossi(c));
% %                 offs = [offs, met(win2)];
% %                 subplot(2,2,2); hold on; plot(wint2, stimi(win2), 'k')
% %                 subplot(2,2,4); hold on; plot(wint2, met(win2), 'k')
% %             end
%             if isempty(oncrossi) & length(offcrossi) == 1
%                 oncrossi = 1; end
%             if isempty(offcrossi) & length(oncrossi) == 1
%                 offcrossi = length(bump); end
%             if isempty(oncrossi) & isempty(offcrossi)
%                 oncrossi = 1; offcrossi = length(bump); end
% 
%             if offcrossi(1) < oncrossi(1)
%                 oncrossi = [1;oncrossi]; end
%             if oncrossi(end) > offcrossi(end)
%                 offcrossi = [offcrossi; length(bump)]; end
% 
%             keep = find((offcrossi - oncrossi) > minon*round(fs)); %gets rid of bumps less than minon duration
%             bump = zeros(length(sdff), 1);
% %             bump(bumpi) = 1;
%             for k = 1:length(keep)
%                 bump(oncrossi(keep(k)):offcrossi(keep(k))) = 1;
%             end
%             bumpi = find(bump);
%             oni = oncrossi(keep);
%             offi = offcrossi(keep);
%             ont = tt(oni);
%             offt = tt(oni);
%         else
%             oni = [];
%             offi = [];
%         end
% % keyboard
