clear all
close all
tic

run_model = 1;
plot_trials = 1;
sum_figs = 1;

save_res = 0;
save_filepath = '/Users/kathmn01/Desktop/datafit_fscmodel_figs/air_tauscreen_n500_strtposx45y50-300_fixed_dope';
save_filepath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/FSC modeling/datafit_fscmodel_figs/rigolli_tauscreen_n500_strtposx50y50-350_150mmgoal_15hzfs_0d5Othresh_5sAtau';
save_filepath = '/Users/kathmn01/Desktop/datafit_fscmodel_figs/air_tauscreen_n500_strtposx45y50-300_fixed_dope';

plot_travelstats = 0;
print_figs = 0;

add2txt = ''; %'_3Ps_5trl';


%% optimize pSearch
subpath = '/Users/kathmn01/Desktop/';
if run_model
    close all

    T = .02 * 2.^(1:15);
    T = .1 * 2.^(0:13);
%     T = .5 * 2.^(0:13);
%     T = [1, 6, 25, 100];
%     T = [1.6, 6.4, 25.6];
%     T = 0.8;

    % Walking
    % avvars = [5, 1.6, .5, 1]; %original values
    avvars = [1.5 0.4 0.1 0.5]; %first guess after unit conv
    avvars = [1 0.25 0.05 0.25]; %kathy fit
    avvars = [1 0.25 0.3 0.2]; %my guess
    avvars = [.7 0.3 0.1 0.2]; %narrow (paper values for walking!!!!!!!)
%     avvars = [1.1 0.5 0.2 0.4]; %wide  
%     avvars = [.7 0.3 0.05 0.3]; %narrow v, fit a

%     % Flight
%     avvars = [40 30 .9 .9]; %wide flight
% %     avvars = [30 20 .6 .6]; %best flight
% %     avvars = [15 15 .4 .4]; %narrow flight
%     avvars = [10 10 0.25 0.25]; %hybrid
%     avvars = [5 5 0.1 0.2]; %hybrid2 - used in first rigoli simulation
%     avvars = [40 30 .9 .9]; %floris flight data - used in second rigoli simulation (keeper!!)


%     start_pos = [50, 0, 50, 350]; %air

    start_pos = [50, 0, 50, 300]; %boundary

    plume_type = 'rigoliplume'; %'boundaryplume'; 'airplume'; 'waterplume'; 'rigoliplume'
    plume_type = 'boundaryplume';% 'airplume'; 'waterplume'; 'rigoliplume'
 
    ntrials = 5;%00;
    samp = 3600; %3600; %500:500:6500;
    rad = 15; %mm
%     rad = 15*(1900/160); %mm (for big plume, rigoli, scaled to proportion of plume size)
%     rad = 15*10; %for rigolli

    allpwin = [];
    
    if plot_travelstats
        allodor = [];
        allC = [];
        alls = [];
        allv = [];
        alla = [];
        all_dist = [];
    end

        for t = 1:length(T)
            disp(['Running ' num2str(t) ' of ' num2str(length(T)) ' tau values'])
            pwin = []; %success rates (time spent at source per trial)
            for i = 1:ntrials
% % %                     res(t, i) = FSCnavmodel_v2(samp(s), 'boundaryplume', '', T(t), 0.4, plot_trials); 
%                 res(t, i) = FSCnavmodel(samp, 'boundaryplume', '', T(t), avvars, plot_trials); 
                [res(t, i), mmperpix] = FSCnavmodel_altplumes(samp, plume_type, '', T(t), avvars, start_pos, plot_trials); 
%                     pwin(i) = pRad(res(e,i).x,res(e,i).y,0,0,rad);
                pwin(i) = pRad(res(t,i).x,res(t,i).y,0,0,rad,mmperpix);
                   
                if plot_travelstats

                    odor = res(t, i).odor;
                    C = res(t, i).C;
                    s = res(t, i).s;
                    v = res(t, i).v;
                    a = res(t, i).a;
                    x = res(t, i).x;
                    y = res(t, i).y;
 

                    figure(20); set(gcf,'Position',[448   363   560   420]); hold on;
                    plot(x,y,'k');
                    plot(x(s==1),y(s==1),'k.');
                    plot(x(s==2),y(s==2),'r.');
                    plot(x(s==3),y(s==3),'b.');
                    axis equal
    %                 xlim([-9000 9000]); ylim([-9000 7000])
                
                    allodor = [allodor; odor];
                    allC = [allC; C];
                    thiss = zeros(1,length(s));
                    thiss(s == 2) = 1;
                    alls = [alls; thiss];
                    allv = [allv; v];
                    alla = [alla; a];
    
                    dx = diff(x);
                    dy = diff(y);
                    dists = sqrt(dx.^2 + dy.^2);
                    all_dist = [all_dist, sum(dists)];
                end


            end
            
            if plot_travelstats
                figure(1); clf; set(gcf,'Position',[78    34   367   751]);
                subplot(5,1,1); plot(mean(allodor)); ylabel('mean odor')
                title(['avvars=' num2str(avvars)])
                subplot(5,1,2); plot(mean(allC));  ylabel('mean C') %hold on; plot(D,'r'); set(gca,'Ylim',[0 1.2]); 
                subplot(5,1,3); plot(mean(alls)); ylabel('mean goal state')
                subplot(5,1,4); plot(mean(allv)); ylabel('mean fvel')
                subplot(5,1,5); plot(mean(alla)); ylabel('mean avel')
    
                figure(3); clf; scatter(ones(1, length(all_dist)), all_dist, 20, 'k', 'filled')
    %             ylim([0 40000])
                title(['avvars=' num2str(avvars)])
            end


    %         figure(1); hold  on; scatter((ones(1,length(pwin))+(randn(1, length(pwin))*.1))*P(p), pwin, 'k', 'filled')
%             allpwin(e,:) = pwin;
            allpwin(t,:) = pwin;
        end

end
if plot_trials && print_figs
    print(gcf, '-depsc', '-painters', [subpath 'FSCmodel_traj_Tau' num2str(T(t)) '_n' num2str(ntrials) add2txt '.eps'])
end

if save_res
    save([save_filepath '.mat'], 'res', 'avvars', 'start_pos', 'allpwin', 'pwin', 'T', '-v7.3')
    save([save_filepath '_no-res.mat'], 'avvars', 'start_pos', 'allpwin', 'pwin', 'T')
end

   %% plot prad for each T
%     figure('Position', [100, 100, 200*size(allpwin,1), 600]); clf; hold  on; 
%     
%     for i = 1:size(allpwin,1)
%         jit = randn(1, length(allpwin(i,:)))*.1;
%         scatter(gca, ones(1,length(allpwin(i,:)))*i+jit, allpwin(i,:), 30, 'k', 'filled')
%         plot([i-.4 i+.4], [mean(allpwin(i,:)) mean(allpwin(i,:))], 'r', 'LineW', 4)
%     end
%     BoxByRow(gca, allpwin, 1:length(E), .5)
%     % Relabel x-ticks
%     xticks(1:length(E)); % Set the x-tick positions to match the box plot positions
%     xticklabels(E); % Set the custom x-tick labels
%     xlabel('tauE')
%     ylabel('% time at source')
%     xlim([-.05 length(E)+.5])
%     title(['n = ' num2str(ntrials)])
%     
%     if print_figs
%         % print_fig(1, gcf, [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt], 'eps', 1); 
%         exportgraphics(gcf, [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.png'])
%         print(gcf, '-depsc', '-painters', [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.eps'])
%     end   

% if sum_figs
    figure('Position', [100, 100, 200*size(allpwin,1), 600]); clf; hold  on; 
    
    for i = 1:size(allpwin,1)
        jit = randn(1, length(allpwin(i,:)))*.1;
        scatter(gca, ones(1,length(allpwin(i,:)))*i+jit, allpwin(i,:), 30, 'k', 'filled')
        plot([i-.4 i+.4], [mean(allpwin(i,:)) mean(allpwin(i,:))], 'r', 'LineW', 4)
    end
    BoxByRow(gca, allpwin, 1:length(T), .5)
    % Relabel x-ticks
    xticks(1:length(T)); % Set the x-tick positions to match the box plot positions
    xticklabels(T); % Set the custom x-tick labels
    xlabel('Tau')
    ylabel('% time at source')
    xlim([-.05 length(T)+.5])
    title(['n = ' num2str(ntrials) ' ,  avvars = ' num2str(avvars)])
    
    if print_figs
        % print_fig(1, gcf, [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt], 'eps', 1); 
        exportgraphics(gcf, [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.png'])
        print(gcf, '-depsc', '-painters', [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.eps'])
    end 
    
    
    %% heatmap by mode
if sum_figs
    add2txt = '';
    % Define grid resolution
    xEdges = -400:15:400; % X grid edges (0 to 10, step 1)
    yEdges = -1200:15:1700; % Y grid edges (0 to 10, step 1)
%     yEdges = -800:15:2300; % Y grid edges (0 to 10, step 1)
    
    whichT = 1:length(T);
    sigma = 1; % Standard deviation for Gaussian smoothing
    figure('Position', [100, 100, 600, 200*length(whichT)]); clf;
%     figure('Position', [100, 100, 600, 300]); clf;
    for t = 1:length(whichT); %1:length(P)
        allx = [];
        ally = [];
        alls = [];
        for i = 1:length(res(whichT(t),:))
            allx = [allx, res(whichT(t),i).x];
            ally = [ally, res(whichT(t),i).y];
            alls = [alls, res(whichT(t),i).s];
        end
    
        ax1 = subplot(length(whichT), 3, t*3-2);
        ax2 = subplot(length(whichT), 3, t*3-1);
        ax3 = subplot(length(whichT), 3, t*3);
        
        % Pass axis handles to the function
        FSCnavmodel_plotStateHeatmaps([ax1, ax2, ax3], allx, ally, alls, xEdges, yEdges, sigma);
    
        subplot(length(whichT), 3, t*3-2);
        ylabel(['tau=' num2str(T(whichT(t)))])
        if t == 1
            title(ax1, 'baseline');end
        if t == 1
            title(ax2, 'goal');end
        if t == 1
            title(ax3, 'search');end
            
    end
    sgtitle(['n = ' num2str(ntrials)])
    if print_figs
        % % print_fig(1, gcf, [subpath 'FSCmodel_stateheatmap_n' num2str(ntrials) add2txt], 'eps', 1); 
        print(gcf, '-depsc', '-painters', [subpath 'FSCmodel_stateheatmap_n' num2str(ntrials) add2txt '.eps'])
        exportgraphics(gcf, [subpath 'FSCmodel_stateheatmap_n' num2str(ntrials) add2txt '.png'])
    end
end
toc

function p = pRad(x, y, centerX, centerY, radius_mm, mmperpix)
    % Function to calculate the percentage of samples within a radius
    % from a specified (centerX, centerY) coordinate.

    % Calculate distances from the center point to all samples
    distances = sqrt((x - centerX).^2 + (y - centerY).^2);

    % Find samples within the radius
    radius = round(radius_mm/mmperpix);
    withinRadius = distances <= radius;

    % Calculate the percentage
    p = (sum(withinRadius) / numel(x));
end
