
function [out, pxscale] = FSCnavmodel(triallength,environment,state, tau, start_pos, plot_figs) 
% out = FSCnavmodel(triallength, environment, plotting, ntrials) 
% navigation model with three states: 
%   1 — baseline
%   2 - goal directed
%   3 — search
%
% (B)aseline and (S)earch states set by Gattuso model
% (G)oal directed can be directly upwind or at an angle to the wind
% switch between states govered by statistics of h∆K bump


%% model parameters
% tauA = 150;
beta = 0.025;
thresh = 0.5;


B = 0.8;           % autocorrelative term for locomotion
XLR = 0;           % left-right weight for baseline locomotion
XLR2 = -0.035;     % left-right weight for local search
XS = -0.03;        % stop weight
XX = 0;            % interaction weight for DNs on same side

W = [B XX XS XLR XLR;       % network weights baseline period
    XX B XS XLR XLR;
     XS XS B XS XS;
    XLR XLR XS B XX;
    XLR XLR XS XX B];

Woff = [B XX XS XLR2 XLR2;       %network weights for offset period
    XX B XS XLR2 XLR2;
    XS XS B XS XS;
    XLR2 XLR2 XS B XX;
    XLR2 XLR2 XS XX B];


%% reminder: X is crosswind (short axis), Y is downwind (long axis)

%% environment parameters
if strcmp(environment, 'rigoliplume')

    avvars = [40 30 .9 .9]; %floris flight data
    fs = 2.78; % Hz of plume, 2.78, wind speed/mean horiz. speed: 25cm/s, at 53cm height
    xpix = 280;
    ypix = 1225;    
    xmm = 3000;
    ymm = 15000;
    nsamp = 2600;

else
    avvars = [.7 0.3 0.1 0.2]; % ball data
    fs = 15; % Hz of plume
    xpix = 216;
    ypix = 406;
    xmm = 159.84;
    ymm = 300.44;
    nsamp = 3600;
end
tauA = fs*5; %5s to match walking ball plume closed loop
pxscale = xmm/xpix; %x=1458mm/2489pix, y=799.3/1364; %mm/pixel ratio to convert pixels from the plume data to actual mm


%% initial values
state = 'baseline';
U = [0.1; 0.1 ;0 ;0.1; 0.1];
v(1) = 0;
a(1) = 0;

A(1) = 0;
C(1) = 0;

%Set start positions (using range and offset values scaled for plume pixel size to replicate positions used with orginal boundary plume data)
% start_pos = [50, 0, 50, 300] %used for walking trials
xrange = round(start_pos(1)*xpix/216); xoffset = round(start_pos(2)*xpix/216); 
yrange = round(start_pos(3)*ypix/406); yoffset = round(start_pos(4)*ypix/406);
x(1) = randn(1)*xrange; % random distribution of initial X positions centered around 0, with sigma=16
y(1) = rand(1)*yrange-yoffset; % random distribution of initial Y positions, between -250 and -300

theta(1) = 2*pi*rand(1); % initial heading

n = randn(5,triallength); %noise component for each neuron


%% compute goal transition probability
p = 1 / (tau * fs); % convert decay time constant to transition probability


%% run model
for i=2:triallength

    switch environment
        case {'boundaryplume'}

            tind = mod(i-1,3600)+1; % Restarts the count in case we want to run longer trials
            xind = round(x(i-1))+108; 
            yind = -round(y(i-1));

            if ismember(xind,[1:216]) && ismember(yind,[1:406]) %set odor value for x,y
                odor(i)=max(0,h5read('/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/Data/10302017_10cms_bounded_2.h5',...
                    '/dataset2', [xind yind tind], [1 1 1]));
            else
                odor(i) = 0;
            end

        case {'rigoliplume'}

            % check for odor
            tind = mod(i-1,nsamp)+1; % Restarts the count in case we want to run longer trials
            xind = round(x(i-1))+round(xpix/2); 
            yind = -round(y(i-1));
            if ismember(xind,[1:xpix]) && ismember(yind,[1:ypix]) %set odor value for x,y

                odor(i)=max(0,h5read('/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/Data/rigoli_9d5x1d9m_10cms_100hz_norm.h5',...
                    '/Data', [yind xind tind], [1 1 1]));
                odor(i) = odor(i)*1;
            else
                odor(i) = 0;
            end
            

        case {'pulse'}
            if (i>451)&(i<600)
                odor(i) = 1;
            else
                odor(i) = 0;
            end

        case {'none'}
            odor(i) = 0;
    end

    % adaptive threshold
    A(i) = A(i-1) + (odor(i) - A(i-1))/tauA; %adaptive odor
    C(i) = round(odor(i)./(odor(i)+beta+A(i))); %compressed odor

    % determine DN activity based on state
    switch state

        case {'baseline'}
            s(i) = 1;
            %U(:,i) = act(W*U(:,i-1) + n(:,i));
            U(:,i) = act(W*U(:,i-1) + 0.1*[1;1;0;-1;-1]*sin(theta(i-1)+pi/2) + n(:,i));

            if C(i)>thresh
                state='goal';
            end

        case {'goal'}
            s(i) = 2;
            U(:,i) = act(W*U(:,i-1) + [1;1;0;-1;-1]*sin(theta(i-1)-pi/2) + [0;0;-.1;0;0] + n(:,i));
            if C(i)<thresh
                if rand(1)<p %0.02 
                    state = 'search';
                end
            end  

        case {'search'}
            s(i) = 3;
            U(:,i) = act(Woff*U(:,i-1) + n(:,i));
            if C(i)>thresh
                state='goal';
            elseif rand(1)<0.01 %orig 0.02
                state = 'baseline';
            end

    end

    % compute forward and angular velocity (/sec)
    if U(3,i)>0
        v(i) = 0;
        a(i) = 0;
    else
        v(i) = (avvars(1)*max(0,sum(U([1 5],i))) + avvars(2)*max(0,sum(U([2 4],i)))); % should be in mm/sec
        a(i) = (avvars(3)*diff(U([1 5],i))+avvars(4)*diff(U([2 4],i)));% should be in rad/sec
    end

    theta(i) = theta(i-1)+a(i)/fs; % convert to rad/samp
    [dx, dy] = pol2cart(theta(i),v(i)/fs); %converting to x and y coordinates and make unit/sec
    x(i) = x(i-1)+dx;
    y(i) = y(i-1)+dy;

end

if plot_figs
    figure(10); set(gcf,'Position',[78    34   367   751]);
    subplot(5,1,1); plot(odor); ylabel('odor')
    subplot(5,1,2); plot(C); 
    ylabel('binary odor')

    subplot(5,1,3); hold off; plot(s==2,'r');  hold on; plot(s==3,'b'); ylabel('states')
    subplot(5,1,4); plot(v); ylabel('fvel')
    subplot(5,1,5); plot(a); ylabel('avel')
    
    figure(2); set(gcf,'Position',[448   363   560*1.5   420*1.5]); hold on;
    % figure(2); clf; set(gcf,'Position',[448   363   560   420]); hold on;
    % plume = h5read('/Users/kathmn01/NYU Langone Health Dropbox/Nicholas Kathman/Data/10302017_10cms_bounded_2.h5', '/dataset2', [1 1 1], [216, 406, 1]);
    % plume = plume'; plume = flipud(plume);
    % imagesc(gca, (1:216)-108, (1:406)-406, plume); cmap = colormap(gca,'gray'); colormap(flipud(cmap));
    plot(x,y,'k');
    plot(x(s==1),y(s==1),'k.');
    plot(x(s==2),y(s==2),'c.');
    plot(x(s==3),y(s==3),'g.');
    plot(x(1), y(1), 'xb', 'LineW', 6)
    % 
    if strcmp(environment, 'waterplume')
        plot([-136, 137, 137, -136, -136], [-512, -512, 0, 0, -512], 'r-', 'LineWidth', 3);
        plot([-136, 0, 137], [-512, 0, -512], 'r-', 'LineWidth', 3)
        circle(0,0,20*(512/405))

    else
        plot([-107, 108, 108, -107, -107], [-405, -405, 0, 0, -405], 'r-', 'LineWidth', 3);
        plot([-107, 0, 108], [-405, 0, -405], 'r-', 'LineWidth', 3)
        circle(0,0,round(15*(xmm/160)/pxscale)) %scales target to relative radius of boundary plume (rad=15mm, plume_xmm = 160)
    end
    % 
    axis equal
end

out.U = U;
out.v = v;
out.a = a;
out.x = x;
out.y = y;
out.odor = odor;
out.C = C;
out.s = s;

end

function y = act(x) %keeps the system from exploding
    y = 20./(1+exp(-x/4))-10;
end




