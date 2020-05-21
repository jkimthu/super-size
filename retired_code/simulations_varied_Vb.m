%%  Simulations: varied Vb

%  Goal: determine whether initial birth size of simulations
%        affects steady-state cell size distribution

%  Conclusion: nope! Steady-state is robust to initial cell size


%  Strategy:
%
%  Part 0: initialize simulations (square wave mu function)
%  Part 1: loop through differet initial Vb and plot simulation results


%  Last edit: jen, 2020 Feb 6
%  commit: square wave, T = 15min, manually varied initial Vi



% Okie, go go let's go!

%% Part 0. initialize analysis
clc
clear

% 0. initialize periods of interest (seconds)
T = [900, 3600];

% 0. initialize timestep in simulation (seconds)
tstep = 15;
tstep_h = tstep/3600;

% 0. initialize parameters for instantaneous growth rate function 
mu_high = 1.9; % units = 1/h
mu_low = 0.7;  % units = 1/h
numstep = T/tstep; % seconds, conversion to h occurs in simulation
qT = numstep/4;    % seconds, conversion to h occurs in simulation

% 0. initialize colors per successive period
%colorSpectrum = {'Indigo','MediumSlateBlue','DodgerBlue','DeepSkyBlue','Teal','DarkGreen','MediumSeaGreen','GoldenRod','DarkOrange','Red'};
colorSpectrum = {'Indigo','DodgerBlue'};%,'DeepSkyBlue','Teal','MediumSeaGreen','GoldenRod','Gold'};
%colorSpectrum = {'SlateGray','DarkCyan','CadetBlue','DeepSkyBlue','DodgerBlue','Navy','Indigo',};


%% Part 1. simulate cells in fluctuations for each timescale 

%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations


% 0. start generation counter
gen = 1000;
gen_counter = 1;

% 0. initialize vectors for data storage
tau = nan(gen,length(T));
Vd = nan(gen,length(T));
Vb = nan(gen,length(T));

% 0. initiation size, timestep, tau
is_newCell = 1; % 0 = division did not just happen
                % 1 = birth timestep! 
                %     record birth size and reset current_tau
Vi = 2;         % first birth size, fixed across conditions
current_tau = -(17/60)*Vi + (64/60);
T_counter = 1; % time counter for growth rate
tsb = 0;       % time since birth for cell cycle
V = [];
t = [];

t_col = 1; % need to make this a loop for multiple timescales

% 1. define instantaneous growth rate function for condition
T_length = numstep(t_col);
T_quarter = qT(t_col);
mufun = [mu_high*ones(T_quarter,1); mu_low*ones(2*T_quarter,1); mu_high*ones(T_quarter,1)];

% 2. simulate!
counter_loop = 0; % each loop is a timestep
while gen_counter < gen+1
    
    counter_loop = counter_loop + 1;
    
    if is_newCell == 1
        % store current cell cycle Vi
        Vb(gen_counter,t_col) = Vi;
    end
    
    % i. calculate exponent (mu * timestep)
    expo = mufun(T_counter)*tstep_h; % mu(t)*t in hours
    
    % ii. calculate new volume after timestep
    if expo < 0
        Vt = Vi;
    else
        Vt = Vi * 2^(expo);
    end
    
    % iii. calculate time into cell cycle
    tsb = tsb + tstep_h;
    
    
    % iv. store timestep data for visualization
    V(counter_loop) = Vt;
    if counter_loop == 1
        t(counter_loop) = tstep_h;
    else
        t(counter_loop) = t(counter_loop-1) + tstep_h;
    end
    
    
    % v. determine whether cell divides
    if tsb >= current_tau 
        % a. division! store V_div and tau data
        Vd(gen_counter) = Vt;
        tau(gen_counter,t_col) = current_tau;
        
        % b. re-set parameters for next cell cycle
        gen_counter = gen_counter + 1;
        Vi = Vt/2;
        current_tau = -(17/60)*Vi + (64/60);
        is_newCell = 1;
        tsb = 0;
    else
        Vi = Vt;
        is_newCell = 0;
    end
    
    
    % vi. determine whether nutrient period restarts
    if T_counter == T_length
        T_counter = 1;
    else
        T_counter = T_counter + 1;
    end
    
end

figure(1)
hold on
plot(t,V)
axis([0 30 0.5 4])
xlabel('Time (h)')
ylabel('Volume')
title('visualization 15 min')
legend('1','2','3','4','0.5')

%% E. simulate cells in 60 min fluctuations, visualization
%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations


% 0. start generation counter
gen = 1000;
g60 = 1;

% 0. initialize vectors for data storage
tau_60 = nan(gen,1);
Vd_60 = nan(gen,1);
Vi_60 = nan(gen,1);

% 0. timestep for each timescale
t60 = 1/binsPerPeriod;


% 0. initiation size, timestep, tau
is_newCell = 1; % 0 = division did not just happen
Vi = 2; %
tau = -(17/60)*Vi + (64/60);
T_counter = 1; % time counter for growth rate
tsb = 0; % time counter for cell cycle
V = [];
t = [];

counter_loop = 0; % each loop is another timestep
while g60 < gen+1
    
    counter_loop = counter_loop + 1;
    
    if is_newCell == 1
        % store current cell cycle Vi
        Vi_60(g60) = Vi;
    end
    
    % calculate exponent
    expo = mufun60(T_counter)*t60; % mu(t)*t in hours
    
    % calculate new volume after timestep
    if expo < 0
        Vt = Vi;
    else
        Vt = Vi * 2^(expo);
    end
    
    % calculate time into cell cycle
    tsb = tsb + t60;
    
    % store timestep data
    V(counter_loop) = Vt;
    if counter_loop == 1
        t(counter_loop) = t60;
    else
        t(counter_loop) = t(counter_loop-1) + t60;
    end
    
    % determine whether cell divides
    
    if tsb >= tau 
        
        % division! store V_div and tau data
        Vd_60(g60) = Vt;
        tau_60(g60) = tau;
        
        % re-set parameters for next cell cycle
        g60 = g60 + 1;
        Vi = Vt/2;
        tau = -(17/60)*Vi + (64/60);
        is_newCell = 1;
        tsb = 0;
        
    else
        Vi = Vt;
        is_newCell = 0;
    end
    
    if T_counter == 30
        T_counter = 1;
    else
        T_counter = T_counter + 1;
    end
    
end

% figure(2)
% plot(t,V)
% %axis([0 40 0.5 3])
% axis([20 40 0.5 2])
% xlabel('Time (h)')
% ylabel('Volume')
% title('visualization 60 min')

clear ans counter_loop experimentFolder expType g15 g60
clear minTime new_cell period_frac specificGrowthRate
clear index growthRates_final muT t t_cc tau
clear timeInPeriodFraction V Vi Vt

%% F. simulation analysis: plot mean Vi, added V vs Vi, tau vs Vi

steadies = 50:1000;

figure(3)
histogram(Vi(steadies),'Facecolor',rgb(colorSpectrum{2}))
hold on
histogram(Vi_60(steadies),'Facecolor',rgb(colorSpectrum{1}))
ylabel('Cell count')
xlabel('Initial volume')
legend('15 min','60 min')


% added V vs initial volume
delta_15 = Vd - Vi;
delta_60 = Vd_60 - Vi_60;

figure(4)
plot(Vi(steadies),delta_15(steadies),'o','Color',rgb(colorSpectrum{2}))
hold on
plot(Vi_60(steadies),delta_60(steadies),'o','Color',rgb(colorSpectrum{1}))
ylabel('Added volume')
xlabel('Initial volume')
legend('15 min','60 min')


% added V binned by birth size
binsize = 0.02; % cubic microns
birth_bin15 = ceil(Vi/binsize);
birth_bin60 = ceil(Vi_60/binsize);

% sort added size by birth bin
adder15_binned = accumarray(birth_bin15(steadies),delta_15(steadies),[],@(x) {x});
adder15_mean = cellfun(@mean, adder15_binned);
adder15_std = cellfun(@std, adder15_binned);
adder15_count = cellfun(@length, adder15_binned);
adder15_sem = adder15_std./sqrt(adder15_count);

adder60_binned = accumarray(birth_bin60(steadies),delta_60(steadies),[],@(x) {x});
adder60_mean = cellfun(@mean, adder60_binned);
adder60_std = cellfun(@std, adder60_binned);
adder60_count = cellfun(@length, adder60_binned);
adder60_sem = adder60_std./sqrt(adder60_count);

figure(5)
errorbar(adder15_mean,adder15_sem,'o','Color',rgb(colorSpectrum{2}))
hold on
errorbar(adder60_mean,adder60_sem,'o','Color',rgb(colorSpectrum{1}))
ylabel('Added volume')
xlabel('Initial volume * 0.02 (bin size)')
legend('15 min','60 min')
axis([34 52 .6 1])

% tau vs birth size
figure(6)
plot(Vi(steadies),tau(steadies)*60,'o','Color',rgb(colorSpectrum{2}))
hold on
plot(Vi_60(steadies),tau_60(steadies)*60,'o','Color',rgb(colorSpectrum{1}))
ylabel('Interdivision time (min)')
xlabel('Initial volume')
legend('15 min','60 min')


