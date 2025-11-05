function out = run_nav_simulation(triallength, tau, plot_figs)
% RUN_NAV_SIMULATION  Simulates odor-guided navigation based on neural dynamics
%
% out = run_nav_simulation(triallength, tau, plot_figs)
%
% This function simulates a Drosophila-like agent navigating within an odor
% plume environment using a 3-state behavioral model (baseline, goal-directed,
% search). The odor concentration is sampled from an external plume dataset
% (Crimaldi 2017, 10 cm/s bounded plume).
%
% INPUTS:
%   triallength : number of samples to simulate (e.g. 1500 for ~100 s at fs=15 Hz)
%   tau         : mean decay time constant (in seconds) controlling transition probability
%   plot_figs   : logical (1 or 0), whether to display output figures
%
% OUTPUT:
%   out : structure with simulation results:
%         out.U  - DN activity matrix (5 x time)
%         out.v  - forward velocity (mm/s)
%         out.a  - angular velocity (rad/s)
%         out.x  - x position (plume pixel coordinates)
%         out.y  - y position (plume pixel coordinates)
%         out.odor - sampled odor concentration
%         out.C  - compressed odor signal
%         out.s  - behavioral state index (1=baseline, 2=goal, 3=search)
%
% Example:
%   out = run_nav_simulation(1800, 5, 1);
%   This runs a 2-minute simulation with tau=5s and displays figures.
%
% Note:
%   You must edit the variable `plume_path` below to point to the correct
%   .h5 plume dataset on your system.
%
% Authors: Nicholas Kathman, Katherine Nagel (2025)
% Manuscript: Neural dynamics for working memory and evidence integration during olfactory navigation in Drosophila

%%  1. SETUP AND PARAMETERS

% ---- Path to plume dataset (edit this path to your local copy) ----
plume_path = '/path/to/10302017_10cms_bounded_2.h5';  % <-- update this line

% ---- Simulation and model parameters ----
beta   = 0.025;      % compression term
thresh = 0.5;        % odor detection threshold

B   = 0.8;           % autocorrelation
XLR = 0;             % left-right weight (baseline)
XLR2 = -0.035;       % left-right weight (search)
XS  = -0.03;         % stop weight
XX  = 0;             % same-side DN interaction

W = [B XX XS XLR XLR; XX B XS XLR XLR; XS XS B XS XS; XLR XLR XS B XX; XLR XLR XS XX B];
Woff = [B XX XS XLR2 XLR2; XX B XS XLR2 XLR2; XS XS B XS XS; XLR2 XLR2 XS B XX; XLR2 XLR2 XS XX B];

%  Environment parameters (Crimaldi dataset)
fs = 15;            % Hz of plume video
xpix = 216; ypix = 406; xmm = 159.84; ymm = 300.44;
nsamp = 3600;       % samples per plume dataset
pxscale = xmm / xpix;

% Behavioral coefficients and initial conditions
avvars = [0.7 0.3 0.1 0.2];  % fitted from experimental data
start_pos = [50, 0, 50, 300]; % [xrange, xoffset, yrange, yoffset]

tauA = fs * 5;                % adaptive odor time constant (samples)
state = 'baseline';           % initial state
U = [0.1; 0.1; 0; 0.1; 0.1];  % initial DN activity
v(1) = 0; a(1) = 0; A(1) = 0; C(1) = 0;

% Initial position and headin
xrange = round(start_pos(1)*xpix/216);
xoffset = round(start_pos(2)*xpix/216);
yrange = round(start_pos(3)*ypix/406);
yoffset = round(start_pos(4)*ypix/406);
x(1) = randn(1)*xrange;
y(1) = rand(1)*yrange - yoffset;
theta(1) = 2*pi*rand(1);

% Random noise and transition probability
n = randn(5, triallength);
p = 1 / (tau * fs);   % convert decay time constant to probability

%% 2. MAIN SIMULATION LOOP
for i = 2:triallength
    tind = mod(i-1, nsamp) + 1;  % restart plume sampling periodically
    xind = round(x(i-1)) + 108;
    yind = -round(y(i-1));

    % Sample odor intensity from plume dataset
    if ismember(xind, 1:216) && ismember(yind, 1:406)
        odor(i) = max(0, h5read(plume_path, '/dataset2', [xind yind tind], [1 1 1]));
    else
        odor(i) = 0;
    end

    % Adaptive odor compression
    A(i) = A(i-1) + (odor(i) - A(i-1)) / tauA;
    C(i) = round(odor(i) ./ (odor(i) + beta + A(i)));

    % Behavioral state updates
    switch state
        case 'baseline'
            s(i) = 1;
            U(:,i) = act(W*U(:,i-1) + 0.1*[1;1;0;-1;-1]*sin(theta(i-1)+pi/2) + n(:,i));
            if C(i) > thresh
                state = 'goal';
            end

        case 'goal'
            s(i) = 2;
            U(:,i) = act(W*U(:,i-1) + [1;1;0;-1;-1]*sin(theta(i-1)-pi/2) + [0;0;-.1;0;0] + n(:,i));
            if C(i) < thresh && rand(1) < p
                state = 'search';
            end

        case 'search'
            s(i) = 3;
            U(:,i) = act(Woff*U(:,i-1) + n(:,i));
            if C(i) > thresh
                state = 'goal';
            elseif rand(1) < 0.01
                state = 'baseline';
            end
    end

    % Compute motion
    if U(3,i) > 0
        v(i) = 0; a(i) = 0;
    else
        v(i) = (avvars(1)*max(0,sum(U([1 5],i))) + avvars(2)*max(0,sum(U([2 4],i))));
        a(i) = (avvars(3)*diff(U([1 5],i)) + avvars(4)*diff(U([2 4],i)));
    end

    theta(i) = theta(i-1) + a(i)/fs;
    [dx, dy] = pol2cart(theta(i), v(i)/fs);
    x(i) = x(i-1) + dx;
    y(i) = y(i-1) + dy;
end

%% 3. PLOTTING
if plot_figs
    % Time series plots
    figure(1); set(gcf, 'Position', [78 34 367 751]);
    subplot(5,1,1); plot(odor); ylabel('odor');
    subplot(5,1,2); plot(C); ylabel('binary odor');
    subplot(5,1,3); hold off; plot(s==2,'r'); hold on; plot(s==3,'b'); ylabel('states');
    subplot(5,1,4); plot(v); ylabel('fvel');
    subplot(5,1,5); plot(a); ylabel('avel');

    % Trajectory plot
    figure(2); set(gcf, 'Position', [448 363 840 630]); hold on;
    plot(x, y, 'k');
    plot(x(s==1), y(s==1), 'k.');
    plot(x(s==2), y(s==2), 'c.');
    plot(x(s==3), y(s==3), 'g.');
    plot(x(1), y(1), 'xb', 'LineW', 6);
    plot([-107,108,108,-107,-107],[-405,-405,0,0,-405],'r-','LineWidth',3);
    plot([-107,0,108],[-405,0,-405],'r-','LineWidth',3);
    axis equal; xlabel('Crosswind (px)'); ylabel('Downwind (px)');
end

%% 4. OUTPUT STRUCTURE
out.U = U;
out.v = v;
out.a = a;
out.x = x;
out.y = y;
out.odor = odor;
out.C = C;
out.s = s;

end

%% HELPER FUNCTION
function y = act(x)
% Activation function preventing runaway excitation
    y = 20 ./ (1 + exp(-x/4)) - 10;
end
