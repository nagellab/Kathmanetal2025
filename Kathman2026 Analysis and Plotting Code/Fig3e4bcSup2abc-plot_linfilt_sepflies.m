
clear all; 
close all

headpath = ['/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/data/all flies/'];
cd(headpath)

allflies_filenames_all

%%
subpath = '/Users/kathmn01/Dropbox (NYU Langone Health)/nagellab/Nick/exploratory figs/lin filt/';
subpath = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/2024 paper/sandbox/exploratory figs/lin filt/';
% use_flies = [21:27, 30:40]; sdm = 0.5;%[21:40, 91:93];%[21:26, 33, 37]; %[21:26,30, 32:35, 37:39]; %21:39;%[6,8:10];%1:17; %[6:10, 12:14, 16:17]; %[1,4,16];%[6:15,17]; %15,17];%, %%% 1,4,16 need debugging
% use_flies = 41:51; sdm = 0; %sdm not used in colornorm

use_flies = [27, 32, 36, 38, 111, 112];%[27, 30, 32, 36, 37, 38, 103, 109:112]; %[27, 32, 36, 38, 111]; %[27, 30, 32, 36, 37]; %[101:103, 109:112]; %[27, 30, 32, 36, 37, 38]; %plume trials (not including fly 37 bc too few trials)
add2txt = '_all6HACKflieswithactivity';
% use_flies = [47:50]; %plume trials
% add2txt = '_all4FC1flies';
% add2txt = '_hACkFC1overlay';
print_figs = 0;
plot_sepfly = 0; %also plots peak metrics scatters

bothlines{1} = [27, 32, 36, 38, 111, 112];
bothlines{2} = [6, 8, 9, 11, 12];
% bothlines{2} = [48, 49, 50, 56, 57];

fig_no = 1;
labls{1} = 'odor --> bump';
labls{2} = 'fvel --> bump';
labls{3} = 'abs(avel) --> bump';
widths = cell(1,2);
aucs = cell(1,2);
for l = 1:2
    use_flies = bothlines{l};
    j = 0;
    syncstruct_all = [];
    allDtrim = cell(1,3);
    figure(fig_no); %clf;
    % set(gcf, 'Position', [1 1 920 550]); %sep flies
    set(gcf, 'Position', [1 1 1420 1025]); 
%     set(gcf, 'Position', [1 1 1420 325]); 
    for f = 1:length(use_flies)
    %     filename{use_flies(f)}
        syncstruct = [];
        load(filename{use_flies(f)})
    
    
        pre = 30; post = 20;
        res = computebumpfilters(syncstruct,pre,post,0.0001);%, add2txt, print_figs, subpath);
        
        pretrim = 29; posttrim = 12; %sec
    %     pretrim = 23; posttrim = 8; %sec
        for d = 1:3
            Dtrim = [];
            Dtrim = res.D(pre*10-pretrim*10:pre*10+posttrim*10,d);
            allDtrim{d} = [allDtrim{d}, Dtrim];
            
    %         figure(1); clf; hold on
            if plot_sepfly
                subplot(length(use_flies),3,d + j*3); hold on
                t = [1:length(Dtrim)]/10-pretrim;
                plot(t,Dtrim, 'k', 'LineW', .5);


                %find filter metrics for comparis
                if d == 3 && l == 1
                    [pk, pki, w] = findpeaks(smoothdata(-Dtrim, 'movmean', 30), 'MinPeakProminence', .0000005, 'Annotate','extents');
                else
%                     [pk, pki, w] = findpeaks(smoothdata(Dtrim, 'movmean', 20), 'MinPeakProminence', .00000005, 'Annotate','extents'); %needed for fly57
                    [pk, pki, w] = findpeaks(smoothdata(Dtrim, 'movmean', 30), 'MinPeakProminence', .000005, 'Annotate','extents');
                end
                if ~isempty(pk)
                    [~, i] = max(pk);
                    if d == 1
                        widths{l} = [widths{l}, w(i)];
                    end

                    if d == 3
    %                     auc = trapz([1:length(Dtrim)]/10-pretrim, smoothdata(Dtrim, 'movmean', 30));
    %                     auc = trapz([1:length(Dtrim)]/10-pretrim, Dtrim);
%                         auc = sum(smoothdata(Dtrim, 'movmean', 30));
                        widthi = (pki(i)-round(w(i)/2)):(pki(i)+round(w(i)/2));
                        auc = trapz(t(widthi), Dtrim(widthi));
                        aucs{l} = [aucs{l}, auc];
                        title(['auc = ' num2str(auc)])
                    end
                end
%                 findpeaks(Dtrim, 'Threshold', 'Annotate','extents')
                plot([xlim], [0 0], '--k', 'LineW', 0.25)
                plot([0 0], [ylim], '--k', 'LineW', 0.25)
                if d == 1
                    ylabel(['Fly ' num2str(use_flies(f))])
                end
                if f ~= length(use_flies)
                    set(gca,'XTick',[]);
                end
            end
    
    % % % %         legend('odor', 'fvel', 'avel', '', '', 'Location', 'NorthWest')
    % % % %         title(['n = ' num2str(j)])
    % % % %         xlabel('Time(s)'); 
        end
        j = j + 1;
    end
    if ~plot_sepfly
        for d = 1:3
            subplot(1,3,d); hold on
            data = allDtrim{d}';
            t = [1:length(Dtrim)]/10-pretrim;
            m = mean(data);
            se = std(data)./sqrt(size(data,1));
            st = std(data);
        
            us = (m + se)/(max(st)-min(st)); %normalized by std
            ls = (m - se)/(max(st)-min(st));
            
            fill([t, t(end:-1:1)], [us, ls(end:-1:1)], 'k', 'FaceAlpha',0.1)
            plot(t, m/(max(st)-min(st)), 'k', 'LineW', 2) %normalized by std

            plot([xlim], [0 0], '--k', 'LineW', 0.5)
            plot([0 0], [ylim], '--k', 'LineW', 0.5)
            title(labls{d}, 'Interpreter', 'None')
        end
    end
    fig_no = fig_no + 1;
end

if print_figs
    if plot_sepfly
        print_fig(1, gcf, [subpath 'linfilt_split' add2txt], 'eps', 1); 
        exportgraphics(gcf, [subpath 'linfilt_split' add2txt '.png'])
    else
        print_fig(1, gcf, [subpath 'linfilt' add2txt], 'eps', 1); 
        exportgraphics(gcf, [subpath 'linfilt' add2txt '.png'])
    end
end
% if plot_sepfly
% 
%     widths{1} = [widths{1} 120]; %%%%%FIX THIS, THIS IS A ROUGH GUESS FOR 3RD HCK FLY, PEAK FINDER ISN'T PICKING UP%%%%
%     
%     figure(200);clf; hold on
%     for i = 1:2
%         scatter(ones(1,length(widths{i}))*i, widths{i}/10, 50, 'k', 'filled')
%         mean(widths{i}/10)
%         plot([i-.1 i+.1], [mean(widths{i}/10), mean(widths{i})/10], 'k', 'LineW', 4)
%     end
%     [~, p] = ttest2(widths{1}, widths{2});
%     xlim([0.5 2.5])
%     ylim([0 18])
%     ylabel('peak width(s)')
%     title(['1=hdck, 2=fc1; p=' num2str(p, '%.3f')])
%     sigstar([1 2], p)
% 
% 
%     figure(201);clf; hold on
%     for i = 1:2
%         scatter(ones(1,length(aucs{i}))*i, aucs{i}/10, 50, 'k', 'filled')
%         mean(aucs{i}/10)
%         plot([i-.1 i+.1], [mean(aucs{i}/10), mean(aucs{i})/10], 'k', 'LineW', 4)
%     end
%     [~, p] = ttest2(aucs{1}, aucs{2});
%     xlim([0.5 2.5])
% %     ylim([0 .1])
%     ylabel('AUC')
%     title(['1=hdck, 2=fc1; p=' num2str(p, '%.3f')])
%     sigstar([1 2], p)    
% end
