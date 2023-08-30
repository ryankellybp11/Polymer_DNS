% Author: Ryan Kelly 
% Last Edited: 7/19/22
% The purpose of this MATLAB script is to take input data from the DNS code
% called "mean_u_data.dat" and visualize how close it is to statistically
% stationary turbulence. 
%
% The code produces two figures: 
%   1) The x- and z-averaged velocity profiles are shown divided into 10
%   time intervals over which they are averaged. The flow can be considered
%   stationary when there is sufficiently little deviation between these
%   profiles.
%   2) The overall mean turbulent profile is compared to the laminar
%   profile. Due to the force balance of channel flow, the shear stress
%   (and therefore the velocity gradient at the wall) must be the same.
%
% The only inputs required are for the domain and Reynolds number of the
% flow

clear variables; close all; clc; format long g
set(0,'DefaultTextInterpreter','latex')


% ======================================================================= %
%                               Inputs                                    %
% ======================================================================= %

% nyp = 129; R_tau = 125; % Re = 125
nyp = 129; R_tau = 180; % Re = 180
% nyp = 193; R_tau = 400; % Re = 400

re = 100; % 1/nu
filename = 'Mean_u_data/mean_u_data_180.dat';

% ----------------------------------------------------------------------- %

YL = 2;
y = (1-cos(((1:nyp) - 1)*pi/(nyp-1)))*YL/2;
A = importfile(filename);
nsteps = length(A)/nyp;
B = zeros(nsteps,nyp);
U_m = zeros(nyp,1);
d = YL/2; % Channel half height
dpdx = -1/d*(2*R_tau/(YL*re))^2;

U0 = re*dpdx*((y.^2)/2 - d*y);
% plot(U0,y)

% Fill array from vector file input
for i = 1:nsteps
    for j = 1:nyp
        B(i,j) = A((i-1)*nyp + j);
    end
end

% Calculate mean for each time interval
dn = floor(nsteps/10);
for j = 1:nyp
    U_m(j) = mean(B(:,j));
    for n = 1:10
        Um(j,n) = mean(B((dn*(n-1)+1):(dn*n),j));
    end        
end

hold on
for i = 1:size(Um,2)
    plot(Um(:,i),y')
    % pause(0.5)
end
title('Development of Mean Velocity Profile')
xlabel('$U$ (cm/s)')
ylabel('$y$ (cm)')

% Compute urms solely for reference
urms = sqrt(1/length(U_m)*sum(U_m.^2));

%%
% Create Laminar/Turbulent profile comparison figure
figure(2)
plot(U_m,y','k','LineWidth',2)
hold on

plot(U0,y,'b','LineWidth',2)
xlabel('$U$ (cm/s)','FontSize',18)
ylabel('$y$ (cm)','FontSize',18)
title('Laminar and Turbulent Velocity Profiles','FontSize',18)
grid on

%% Autocorrelation of mean velocity profile
% clc; close all;
% C = B(99001:end,:);
% r = xcorr(C(:,1),'normalized');