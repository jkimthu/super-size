%%  figure 6 - overlaid growth rate from successive periods


%  Goal: is timescale of growth rate fluctuation the cause of different
%        cell size?

%  Conclusion from Part 1: no, fluctuation timescale determines fluctuation
%                          in size but mean Vb and tau are not too different




%  Strategy:
%
%  Part 0: initialize simulations (square wave mu function)
%  Part 1: simulate square wave mu function for each timescale
%  Part 2: simulate mu function with fixed transition time in mu


%  Last edit: jen, 2020 Feb 6
%  commit: simplify growth rate signal to square wave
% 


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


%% Part 1. simulate cells in square fluctuations for each timescale 

%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations



% 0. initialize vectors for data storage
gen = 1000;
tau_mean = nan(gen,length(T));
Vd = nan(gen,length(T));
Vb = nan(gen,length(T));

for t_col = 1:length(T) % need to make this a loop for multiple timescales
    
    % 0. start generation counter
    gen_counter = 1;
    
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
    
    % 1. define instantaneous growth rate function for condition (SQUARE)
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
            Vd(gen_counter,t_col) = Vt;
            tau_mean(gen_counter,t_col) = current_tau;
            
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
    
    figure(t_col)
    hold on
    plot(t,V)
    axis([0 30 0.5 3.2])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('visualization:',num2str(T(t_col))))
    
end

mean(Vb)
mean(Vd)
mean(tau_mean)

%% Part 2. simulate cells in fixed transition fluctuations for each timescale 

%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations

% 0. initialize vectors for data storage
gn = 1000;
tau_fixed = nan(gn,length(T));
Vd_fixed = nan(gn,length(T));
Vb_fixed = nan(gn,length(T));

for t_col = 1:length(T) % need to make this a loop for multiple timescales
    
    % 0. start generation counter
    gn_counter = 1;
    
    % 0. initiation size, timestep, tau
    is_newCell = 1; % 0 = division did not just happen
    % 1 = birth timestep!
    %     record birth size and reset current_tau
    Vi_fixed = 2;         % first birth size, fixed across conditions
    cur_tau = -(17/60)*Vi_fixed + (64/60);
    T_counter = 1; % time counter for growth rate
    tsb = 0;       % time since birth for cell cycle
    V = [];
    t = [];
    
    % 1. define instantaneous growth rate function for condition (FIXED)
    fixed_trans_time = 2*60; % 2 min = 120 seconds
    fixed_tsteps = fixed_trans_time/tstep;
    T_length = numstep(t_col);
    T_quarter = qT(t_col);
    T_remainder = T_quarter - fixed_tsteps;
    mu_trans = linspace(mu_low,mu_high,fixed_tsteps+1);
    mufun = [mu_high*ones(T_quarter,1); mu_low*ones(2*T_quarter,1); mu_trans(2:end)'; mu_high*ones(T_remainder,1)];
    mu_functions{t_col} = mufun;
    
    
    % 2. simulate!
    counter_loop = 0; % each loop is a timestep
    while gn_counter < gn+1
        
        counter_loop = counter_loop + 1;
        
        if is_newCell == 1
            % store current cell cycle Vi
            Vb_fixed(gn_counter,t_col) = Vi_fixed;
        end
        
        % i. calculate exponent (mu * timestep)
        expo = mufun(T_counter)*tstep_h; % mu(t)*t in hours
        
        % ii. calculate new volume after timestep
        if expo < 0
            Vt_fixed = Vi_fixed;
        else
            Vt_fixed = Vi_fixed * 2^(expo);
        end
        
        % iii. calculate time into cell cycle
        tsb = tsb + tstep_h;
        
        
        % iv. store timestep data for visualization
        V(counter_loop) = Vt_fixed;
        if counter_loop == 1
            t(counter_loop) = tstep_h;
        else
            t(counter_loop) = t(counter_loop-1) + tstep_h;
        end
        
        
        % v. determine whether cell divides
        if tsb >= cur_tau
            % a. division! store V_div and tau data
            Vd_fixed(gn_counter,t_col) = Vt_fixed;
            tau_fixed(gn_counter,t_col) = cur_tau;
            
            % b. re-set parameters for next cell cycle
            gn_counter = gn_counter + 1;
            Vi_fixed = Vt_fixed/2;
            cur_tau = -(17/60)*Vi_fixed + (64/60);
            is_newCell = 1;
            tsb = 0;
        else
            Vi_fixed = Vt_fixed;
            is_newCell = 0;
        end
        
        
        % vi. determine whether nutrient period restarts
        if T_counter == T_length
            T_counter = 1;
        else
            T_counter = T_counter + 1;
        end
        
    end
    
    figure(t_col)
    hold on
    plot(t,V)
    axis([0 30 0.5 3.2])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('visualization:',num2str(T(t_col))))
    
end

figure(3)
histogram(Vb_fixed(50:1000,1))
hold on
histogram(Vb_fixed(50:1000,2))
legend('15','60')
xlabel('Birth volume')
ylabel('Cell count')

birthVol_mean = mean(Vb_fixed)
tau_mean = mean(tau_fixed)
lambda_mean = cellfun(@mean,mu_functions)

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
plot(Vi(steadies),tau_mean(steadies)*60,'o','Color',rgb(colorSpectrum{2}))
hold on
plot(Vi_60(steadies),tau_60(steadies)*60,'o','Color',rgb(colorSpectrum{1}))
ylabel('Interdivision time (min)')
xlabel('Initial volume')
legend('15 min','60 min')


