clear all; close all;
%This script does the fit of the OU and Matern model for 200 drifter
%trajectories selected for their large relative change in Coriolis
%frequency, which is likely to make the time series of their velocities
%non-stationary.

%Sampling step, i.e. time between two consecutive observations. In hours
%here.
delta = 2; 
%Optimization parameters
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000,...
    'TolFun',1e-10,'TolX',1e-10,'Display','on');
% Include zero-frequency? 1=yes, 0=no;
ZEROF=0; 
%We make a progress bar
progressBar = waitbar(0, 'Loading the data');
%Loading the data
load blurreddrifters.mat
for drifter_id = 1:200
    %We go through all 200 drifter trajectories one by one
    %Update the progress bar percentage
    waitbar(drifter_id/200, progressBar, ...
        ['Processing drifter trajectories... ' ...
        num2str(drifter_id) '/200']);
    %Time series of velocities of that drifter
    X = blurreddrifters.cv{drifter_id}; 
    %Time series of latitudes
    drifter_lats = blurreddrifters.lat{drifter_id};
    %Coriolis frequency corresponding to the mean latitude. Used for the
    %stationary method.
    CF = coriolis_frequency(mean(drifter_lats));
    %Time series of Coriolis frequencies over the whole trajectory
    CF_t = coriolis_frequency(drifter_lats);
    %Define the minimum and maximum frequencies to use for the fit. These
    %need to be in radian per hour as the frequencies of the spectrum will
    %be in that unit.
    LPCNT = convertFrequency(-0.8, 'cycles per day', 'radians per hour'); 
    UPCNT = convertFrequency(0.8, 'cycles per day', 'radians per hour');
    % Only fit to one side (with the inertial oscillation)
    if CF > 0 
        LPCNT=0;
    else
        UPCNT=0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRELIMINARIES
    N = length(X); 
    omega = Fourier_frequencies(N, delta);
    MF = floor(N/2)+1; 
    %Indices of frequencies in the estimation range. Note that LPCNT and
    %omega should be in the same units.
    LB = frequenciesToIndices(LPCNT, omega);
    UB = frequenciesToIndices(UPCNT, omega);
    %Periodogram
    SZ = 1 / N / delta *(abs(fft(X))).^2;
    SZ = fftshift(SZ); 
    %% INITIAL PARAMETERS (same as in Sykulski et. al (2016) JRSSc for stationary drifters)
    %xb will contain initial values for the optimization algorithms.
    xb=zeros(1,6);
    valMF = SZ(MF);
    % set slope xb(5) to 1
    xb(5)=1;
    valMAX = max(SZ);
    %Estimate location of inertial peak
    IF = round(N * 0.5 * (1 + 0.5 * CF / omega(end))); 
    if CF > 0 
        [~,fmax] = max(SZ(IF:N)); 
        fmax=fmax+IF-1;
    else
        [~,fmax] = max(SZ(1:IF));
    end 
    NF=floor((abs(MF-fmax))/3); % Number of frequencies we now  analyze over
    xb(6)=quantile(SZ([MF-NF:MF-1 MF+1:MF+NF]).*(omega([MF-NF:MF-1 MF+1:MF+NF]))'.^2./(valMAX-SZ([MF-NF:MF-1 MF+1:MF+NF])),0.5); % solve simultaneous equations for xb(4) and xb(6), take median
    xb(4) = valMAX*xb(6); 
    xb(4) = abs(xb(4))^0.5; 
    xb(6) = abs(xb(6))^0.5; % square root of these to arrive at starting parameters
    valmax = SZ(fmax); 
    xb(2) = CF; % corresponding frequency peak in radians
    xb(3)=quantile(SZ([fmax-NF:fmax-1 fmax+1:fmax+NF]).*(omega([fmax-NF:fmax-1 fmax+1:fmax+NF])-xb(2))'.^2./(valmax-SZ([fmax-NF:fmax-1 fmax+1:fmax+NF])),0.5); % solve simultaneous equations for xb(1) and xb(3), take median
    xb(1)=valmax*xb(3); xb(1) = abs(xb(1))^0.5;
    xb(3) = abs(xb(3))^0.5; % square root of these to arrive at starting parameters
    %% Parameter estimation.
    %%The stationary likelihood of Sykulski et al. (2016) JRSSc
    %Minimum and maximum parameter ranges
    minbound = [0  1  pi*sqrt(3)/(N*xb(3)) 0 0.5/xb(5)...
        pi*sqrt(3)/(N*xb(6))]; 
    maxbound = [inf 1 inf inf 2.5 / xb(5) inf]; 
    [x1b, fval1, exitflag1]=fminsearchbnd(@(x) maternOUmodel(x,xb,SZ',N, ...
        LB,UB,MF,ZEROF, delta), ones(1,6),minbound,maxbound, options); 
    x1 = x1b .* xb;
    %%Our nonstationary likelihood
    %Same initialization for the optimization algorithm
    xb2 = xb;
    %Same constraints, except for the second parameter which is set to
    %zero. Indeed, for the non-stationary likelihood, the second parameter
    %represents the average shift to the coriolis frequency. Here we assume
    %no shift.
    minbound = [0  0  pi*sqrt(3)/(N*xb2(3)) 0 0.5/xb2(5) ...
        pi*sqrt(3)/(N*xb2(6))];
    maxbound = [inf 0 inf inf 2.5/xb(5) inf]; 
    [x2, fval2, exitflag2, time2, ker] = fitOUandMatern( SZ', CF_t, ...
        delta, LB, UB, MF, ZEROF, xb2, minbound, maxbound, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stat.ml(drifter_id) = fval1;
    nstat.ml(drifter_id) = fval2;
    stat.param(drifter_id, :) = x1;
    nstat.param(drifter_id, :) = x2;
end
%Update the progressbar and close it
waitbar(1, progressBar, 'Done... Preparing figures...');
pause(2);
close(progressBar);
%% FIGURE 3 IN PAPER
%Calculate the median of estimates of the damping timescale for both
%methods, in days (hence we multiply by 24)
median_stationary = median(1 ./ (stat.param(:, 3) * 24));
median_nonstationary = median(1 ./(nstat.param(:, 3) * 24));
%create figure
Arthur3 = figure('Position',[1   1   scrsz(3) / 1   scrsz(3) / 2.5]); 
%first subplot, contains the scatter of estimates of the damping time scale
subplot(1, 2, 1);
set(gcf, 'renderer', 'Painters');
scatter(1 ./ (stat.param(:, 3) * 24), 1 ./ (nstat.param(:,3) * 24))
axis square;
xlim([0 10]);
ylim([0 10]); 
xlabel(['1/\lambda (days, stationary fit, median = ' ...
    num2str(median_stationary) ')']);
ylabel(['1/\lambda (days, non-stationary fit, median = '...
    num2str(median_nonstationary) ' )']); 
title('(a)'); 
box on
hold on;
line([0 10],[0 10]);
%Second subplot, contains the histogram of the differences in the
%pseudo-likelihoods of both methods for each velocity time series.
subplot(1,2,2); 
histogram(stat.ml-nstat.ml, 50);
title('(b)');
ylabel('frequency'); 
xlabel('$l_S(\hat{\theta}_S)-l_{NS}(\hat{\theta}_{NS})$', 'Interpreter', 'Latex'); 
%Exporting for the paper
exportfig(Arthur3, 'Arthur3.eps', 'color', 'cmyk','Fontmode','fixed','FontSize', 12); 