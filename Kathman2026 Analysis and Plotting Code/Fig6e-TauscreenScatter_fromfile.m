clear all
tic

plot_state_heatmaps = 0;
print_figs = 0;

tau_star = 7.8;
tau_star = 5.59;
y_star = 0.65;

add2txt = '';
subpath = '/Users/kathmn01/Desktop/';

% load('/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/FSC modeling/datafit_fscmodel_figs/boundary_tauscreen_n500_strtposx50y50-300_fixed_no-res.mat')
load('/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/FSC modeling/datafit_fscmodel_figs/boundaryair_tauscreen_n500_strtposx50y50-300_15mmgoal_15hzfs_0d5Othresh_5sAtau_no-res2.mat')
ntrials = length(pwin);

figure('Position', [100, 100, 200*size(allpwin,1), 600]); clf; hold on;

for i = 1:size(allpwin,1)
    jit = randn(1, length(allpwin(i,:))) * 0.1;
    scatter(gca, ones(1,length(allpwin(i,:))) * i + jit, allpwin(i,:), 30, 'k', 'filled')
    plot([i-.4 i+.4], [mean(allpwin(i,:)) mean(allpwin(i,:))], 'r', 'LineWidth', 4)
end

BoxByRow(gca, allpwin, 1:length(T), .5)

% Plot red star 
x_star = log2(tau_star / T(1)) + 1;
plot(x_star, y_star, 'r*', 'LineWidth', 2, 'MarkerSize', 12)

xticks(1:length(T))
xticklabels(T)
xlabel('Tau')
ylabel('% time at source')
xlim([0.5 length(T)+0.5])
title(['n = ' num2str(ntrials) ' ,  avvars = ' num2str(avvars)])

if print_figs
    exportgraphics(gcf, [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.png'])
    print(gcf, '-depsc', '-painters', [subpath 'FSCmodel_pWin_n' num2str(ntrials) add2txt '.eps'])
end

toc



% clear all
% close all
% 
% %% ---------------- USER SETTINGS ----------------
% file_to_load = '/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/nagellab/Nick/FSC modeling/datafit_fscmodel_figs/rigolli_tauscreen_n500_strtposx50y50-350_150mmgoal_15hzfs_0d5Othresh_5sAtau_no-res.mat';
% 
% print_figs = 0;
% subpath = '/Users/kathmn01/Desktop/';
% add2txt = '';
% 
% point_size = 30;
% jit_sd = 0.1;
% mean_line_color = 'r';
% point_color = 'k';
% 
% %% ---------------- LOAD SAVED DATA ----------------
% S = load(file_to_load);
% 
% % Pull required variables
% allpwin = S.allpwin;
% T = S.T;
% 
% % Optional metadata if present
% if isfield(S,'avvars')
%     avvars = S.avvars;
% else
%     avvars = [];
% end
% 
% if isfield(S,'start_pos')
%     start_pos = S.start_pos;
% else
%     start_pos = [];
% end
% 
% % Infer number of trials from saved matrix
% ntrials = size(allpwin,2);
% 
% %% ---------------- PLOT: SCATTER + MEAN + BOX ----------------
% figure('Position', [100, 100, 200*size(allpwin,1), 600]); clf; hold on;
% 
% for i = 1:size(allpwin,1)
%     jit = randn(1, size(allpwin,2)) * jit_sd;
%     scatter(gca, ones(1,size(allpwin,2))*i + jit, allpwin(i,:), ...
%         point_size, point_color, 'filled');
%     
%     plot([i-.4 i+.4], [mean(allpwin(i,:)) mean(allpwin(i,:))], ...
%         mean_line_color, 'LineWidth', 4);
% end
% 
% % Keep this if BoxByRow is already in your path and you want the same look
% BoxByRow(gca, allpwin, 1:length(T), .5);
% 
% % Axis labels/ticks
% xticks(1:length(T));
% xticklabels(T);
% xlabel('Tau');
% ylabel('% time at source');
% xlim([-.05 length(T)+.5]);
% 
% % Title
% if ~isempty(avvars)
%     title(['n = ' num2str(ntrials) ' ,  avvars = ' num2str(avvars)]);
% else
%     title(['n = ' num2str(ntrials)]);
% end
% 
% %% ---------------- OPTIONAL SAVE ----------------
% if print_figs
%     [~, base_name, ~] = fileparts(file_to_load);
% 
%     exportgraphics(gcf, [subpath base_name '_pWinFromSaved' add2txt '.png']);
%     print(gcf, '-depsc', '-painters', ...
%         [subpath base_name '_pWinFromSaved' add2txt '.eps']);
% end
