clear all; close all;
%This script does the fit of the OU and Matern model to the GDP drifters.
%Sampling step (time between two consecutive observations)
delta = 2;
%Optimization algorithm parameters
options=optimset('GradObj','on','MaxFunEvals',100000,'MaxIter',10000, ...
    'TolFun',1e-10,'TolX',1e-10,'Display','on'); 
% Include zero-frequency? 1=yes, 0=no;
ZEROF=1; 
% Figure sub-headings
jj=1; 
FFF = ['a','b','c','d']; 
%Loading the data
load blurreddrifters.mat
load drifterulysses.mat
% The Figure size
scrsz = get(0,'ScreenSize'); 
Arthur2 = figure('Position',[1 1 scrsz(3)/2 scrsz(3)/2.6])
for drifter_id =  [85 123 149 201] 
    % these are the drifter IDs in the Figure example
    if drifter_id == 201
        %This drifter's velocity time series can be considered as
        %stationary.
        % Time series of velocities
        X = drifterulysses.cv(1:852);
        %Time series of latitudes used to compute Coriolis frequencies
        drifter_lats = drifterulysses.lat(1:852);
        %"average" Coriolis frequency, in radians/hour   
        CF = coriolis_frequency(mean(drifter_lats)); 
        % Coriolis frequencies, in radians/hour
        CF_t = coriolis_frequency((drifter_lats)); 
        %Minimal and maximal frequencies used for the fit. In cycles per
        %day here so we wil have to convert.
        LPCNT = convertFrequency(-1.5, 'cycles per day', 'radians per hour'); 
        UPCNT = convertFrequency(1.5, 'cycles per day', 'radians per hour'); 
    else
        X = blurreddrifters.cv{drifter_id};                 
        drifter_lats = blurreddrifters.lat{drifter_id};
        CF = coriolis_frequency(mean(drifter_lats));
        CF_t = coriolis_frequency(drifter_lats);
        LPCNT = convertFrequency(-0.8, 'cycles per day', 'radians per hour'); 
        UPCNT = convertFrequency(0.8, 'cycles per day', 'radians per hour'); 
    end
    % Only fit to one side (with the inertial oscillation)
    if CF>0 
        LPCNT=0;
    else
        UPCNT=0;
    end
    %Quick analysis of the sequence of Coriolis frequencies
    CFmin = min(CF_t);
    CFmax = max(CF_t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PRELIMINARIES
    N = length(X); 
    %Fourier frequencies. Unit of omega depends on the unit of delta. For
    %instance if delta is 2 hours, the Fourier frequencies returned are in
    %radians/h.
    omega = Fourier_frequencies(N, delta);
    MF = floor(N/2)+1; 
    %Indices of frequencies in the estimation range. First we convert to
    %radians per hour.
    LB = frequenciesToIndices(LPCNT, omega);
    UB = frequenciesToIndices(UPCNT, omega);
    SZ=1/N/delta*(abs(fft(X))).^2; 
    SZ=fftshift(SZ); % Power Spectrum
    %% INITIAL PARAMETERS (same as in Sykulski et. al (2016) JRSSc for stationary drifters)
    xb=zeros(1,6);
    valMF = SZ(MF); 
    xb(5)=1; 
    valMAX = max(SZ); 
    IF=round(N*0.5*(1+0.5*CF/pi*delta)); % set slope xb(5) to 1
    if CF > 0 % find location of inertial peak
        [dum1,fmax] = max(SZ(IF:N)); fmax=fmax+IF-1;
    else
        [dum1,fmax] = max(SZ(1:IF));
    end 
    NF=floor((abs(MF-fmax))/3); % Number of frequencies we now analyze over
    xb(6)=quantile(SZ([MF-NF:MF-1 MF+1:MF+NF]).*(omega([MF-NF:MF-1 MF+1:MF+NF]))'.^2./(valMAX-SZ([MF-NF:MF-1 MF+1:MF+NF])),0.5); % solve simultaneous equations for xb(4) and xb(6), take median
    xb(4)=valMAX*xb(6); 
    xb(4) = abs(xb(4))^0.5; 
    xb(6) = abs(xb(6))^0.5; % square root of these to arrive at starting parameters
    valmax = SZ(fmax);
    %Inertial frequency initialized to the mean coriolis frequency
    xb(2)=CF; 
    xb(3)=quantile(SZ([fmax-NF:fmax-1 fmax+1:fmax+NF]).*(omega([fmax-NF:fmax-1 fmax+1:fmax+NF])-xb(2))'.^2./(valmax-SZ([fmax-NF:fmax-1 fmax+1:fmax+NF])),0.5); % solve simultaneous equations for xb(1) and xb(3), take median
    xb(1)=valmax*xb(3); xb(1) = abs(xb(1))^0.5; 
    xb(3) = abs(xb(3))^0.5; % square root of these to arrive at starting parameters
    %% LIKELIHOOD1----------------------------------------------
    %%The stationary likelihood of Sykulski et al. (2016) JRSSc
    disp('Stationary likelihood estimation...');
    tic;
    %Parameter range for the constrained optimization
    minbound = [0  1  pi*sqrt(3)/(N*xb(3)) 0 0.5/xb(5) pi*sqrt(3)/(N*xb(6))]; 
    maxbound = [inf 1 inf inf 2.5/xb(5) inf];
    %Optimization
    [x1b,fval1,exitflag1]=fminsearchbnd(@(x) maternOUmodel(x,xb,SZ',N,LB,UB,MF,ZEROF, delta), ones(1,6),minbound,maxbound, options); 
    time1=toc;
    %%Our nonstationary likelihood
    disp('Nonstationary likelihood method estimation...');
    %Same initialization for the optimization algorithm
    xb2 = xb;
    %Same constraints, except for the second parameter which is set to
    %zero. Indeed, for the non-stationary likelihood, the second parameter
    %represents the average shift to the coriolis frequency. Here we assume
    %no shift.
    minbound = [0  0  pi*sqrt(3)/(N*xb2(3)) 0 0.5/xb2(5) pi*sqrt(3)/(N*xb2(6))];
    maxbound = [inf 0 inf inf 2.5/xb(5) inf]; 
    [x2, fval2, exitflag2, time2, ker] = fitOUandMatern( SZ', CF_t, delta, LB, UB, MF, ZEROF, xb2, minbound, maxbound, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% DISPLAY SOME OUTPUT
    x1 = x1b .* xb;
    disp('Likelihood 1 (stationary likelihood method):');
    disp(['Estimate: ', num2str(x1)]);
    disp(['CPU time: ', num2str(time1)]);
    disp('Likelihood 2 (nonstationary likelihood method):');
    disp(['Estimate: ', num2str(x2)]);
    disp(['CPU time: ', num2str(time2)]);
    %% FIGURE 2 IN PAPER
    %Convert frequencies to cycles per day for the figure
    omega = convertFrequency(omega, 'radians per hour', 'cycles per day');
    CF = convertFrequency(CF, 'radians per hour', 'cycles per day');
    CFmin = convertFrequency(CFmin, 'radians per hour', 'cycles per day');
    CFmax = convertFrequency(CFmax, 'radians per hour', 'cycles per day');
    %Colours
    colour3 = [0 0 0];
    colour2 = [0.5 0.5 0.5];
    colour1 = [0.5 0.5 0.5];
    subplot(2,2,jj);
    %Plot of the periodogram
    plot(omega,10*(log10(SZ)),'color', colour1); 
    xlim([-pi pi]);
    %stationary fit
    acv=maternacvs(x1(4:6),N,delta)+complexouacvs(x1(1:3),N,delta); % predicted autocovariance sequence
    ESF2=2*fft(acv.*(1-(0:N-1)/N))-acv(1); 
    ESF3=abs(real(fftshift(ESF2))); % fft of acvs with triangle kernel to find expected blurred spectrum
    %non-stationary fit
    ESF_=S_(ker, x2, delta);
    hold on; plot(omega,10*log10(ESF3),'color',colour2,'linestyle','--','linewidth',2); % stationary expected periodogram
    hold on; plot(omega,10*log10(ESF_),'color',colour3,'linestyle','--','linewidth',2); % nonstationary expected periodogram
    hold on; plot(omega(LB:UB),10*log10(ESF3(LB:UB)),'color', colour2,'linewidth',2); % stationary expected periodogram (modelled range)
    hold on; plot(omega(LB:UB),10*log10(ESF_(LB:UB)),'color', colour3,'linewidth',1.5); % nonstationary expected periodogram (modelled range)
    hold on; line([CF CF],[-10 60],'color','k'); % average Coriolis
    hold on; line([CFmin CFmin] , [-10 60],'color',colour2,'linestyle','--'); % minimum Coriolis
    hold on; line([CFmax CFmax] , [-10 60],'color',colour2,'linestyle','--'); % maximum Coriolis
    xlabel('\omega (cpd)', 'FontSize', 19); ylabel('PSD (10 log_{10} m^2s^{-2} cpd^{-1})', 'FontSize', 19);
    set(findall(gcf,'type','text'),'fontSize',19,'fontWeight','normal');
    if drifter_id == 201
        xlim([-1.5 1.5]);
    else
        xlim([-1 1]); 
    end
    ylim([-10 60]);
    title(['(' num2str(FFF(jj)) ')']);  jj=jj+1;
end
%%
exportfig(Arthur2, 'Arthur2.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 19); % For exporting fig into paper