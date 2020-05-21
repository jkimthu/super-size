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


%  Last edit: jen, 2020 Feb 12
%  commit: extend to apply square waves in mu to all timescales
% 


% Okie, go go let's go!

%% Part 1. simulate T = 30 sec with square growth rate waves

clc
clear

% 0. initialize periods of interest (seconds)
T = 30;
t30 = 1;

% 0. initialize timestep in simulation (seconds)
tstep = 15;
tstep_h = tstep/3600;

% 0. initialize parameters for instantaneous growth rate function 
mu_high = 1.9; % units = 1/h
mu_low = 0.7;  % units = 1/h
numstep = T/tstep; % seconds, conversion to h occurs in simulation
qT = numstep/4;    % seconds, conversion to h occurs in simulation


% 0. initialize vectors for data storage
gen = 1000;
tau_mean = nan(gen,length(T));
Vd = nan(gen,length(T));
Vb = nan(gen,length(T));


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
T_length = numstep;
T_quarter = qT;
mufun = [mu_high; mu_low];
mu_functions{1} = mufun;


% 2. simulate!
counter_loop = 0; % each loop is a timestep
while gen_counter < gen+1
    
    counter_loop = counter_loop + 1;
    
    if is_newCell == 1
        % store current cell cycle Vi
        Vb(gen_counter) = Vi;
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
        tau_mean(gen_counter) = current_tau;
        
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

% 3. plot timescale data!
figure(1)
hold on
plot(t,V)
axis([0 30 0.5 3.2])
xlabel('Time (h)')
ylabel('Volume')
title(strcat('visualization:',num2str(T)))


figure(2)
histogram(Vb(50:end,1),'FaceColor',rgb('DodgerBlue'))
legend('30 sec')
xlabel('Birth volume')
ylabel('Cell count')
title(strcat('histogram:',num2str(T)))


% 4. store timescale data
mufun_fluc{t30} = mufun;
birthSize_fluc{t30} = Vb;
divSize_fluc{t30} = Vd;
tau_fluc{t30} = tau_mean;


%% Part 2A. simulate cells in fixed transition fluctuations for each timescale 

%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations

% 0. initialize periods of interest (seconds)
T = [300, 900, 3600];
t300 = 2; t900 = 3; t3600 = 4;

% 0. initialize timestep in simulation (seconds)
tstep = 15;
tstep_h = tstep/3600;

% 0. initialize parameters for instantaneous growth rate function 
mu_high = 1.9; % units = 1/h
mu_low = 0.7;  % units = 1/h
numstep = T/tstep; % seconds, conversion to h occurs in simulation
qT = numstep/4;    % seconds, conversion to h occurs in simulation

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
    T_length = numstep;
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
        if T_counter == T_length(t_col)
            T_counter = 1;
        else
            T_counter = T_counter + 1;
        end
        
    end
    
    figure(3)
    hold on
    plot(t,V)
    axis([0 30 0.5 3.2])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('visualization:',num2str(T(t_col))))
    
end
figure(3)
legend('5','15','60')
title('fixed transition')

figure(4)
histogram(Vb_fixed(50:end,1),'FaceColor',rgb('CadetBlue'))
hold on
histogram(Vb_fixed(50:end,2),'FaceColor',rgb('DodgerBlue'))
hold on
histogram(Vb_fixed(50:end,3),'FaceColor',rgb('Indigo'))
legend('5','15','60')
xlabel('Birth volume')
ylabel('Cell count')
title('fixed transition')


% 4. store timescale data
mufun_fluc{t300} = mu_functions{1};
birthSize_fluc{t300} = Vb_fixed(:,1);
divSize_fluc{t300} = Vd_fixed(:,1);
tau_fluc{t300} = tau_fixed(:,1);

mufun_fluc{t900} = mu_functions{2};
birthSize_fluc{t900} = Vb_fixed(:,2);
divSize_fluc{t900} = Vd_fixed(:,2);
tau_fluc{t900} = tau_fixed(:,2);

mufun_fluc{t3600} = mu_functions{3};
birthSize_fluc{t3600} = Vb_fixed(:,3);
divSize_fluc{t3600} = Vd_fixed(:,3);
tau_fluc{t3600} = tau_fixed(:,3);


%% Part 2B. simulate T = 300, 900, 3600 sec with square growth rate waves

%     start with same Vi
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations

clear
clc

% 0. initialize periods of interest (seconds)
T = [300, 900, 3600];
t300 = 2; t900 = 3; t3600 = 4;

% 0. initialize timestep in simulation (seconds)
tstep = 15;
tstep_h = tstep/3600;

% 0. initialize parameters for instantaneous growth rate function 
mu_high = 1.9; % units = 1/h
mu_low = 0.7;  % units = 1/h
numstep = T/tstep; % seconds, conversion to h occurs in simulation
qT = numstep/4;    % seconds, conversion to h occurs in simulation

% 0. initialize vectors for data storage
gn = 1000;
tau_sq = nan(gn,length(T));
Vd_sq = nan(gn,length(T));
Vb_sq = nan(gn,length(T));

for t_col = 1:length(T) % need to make this a loop for multiple timescales
    
    % 0. start generation counter
    gn_counter = 1;
    
    % 0. initiation size, timestep, tau
    is_newCell = 1; % 0 = division did not just happen
    % 1 = birth timestep!
    %     record birth size and reset current_tau
    Vi_sq = 2;         % first birth size, fixed across conditions
    cur_tau = -(17/60)*Vi_sq + (64/60);
    T_counter = 1; % time counter for growth rate
    tsb = 0;       % time since birth for cell cycle
    V = [];
    t = [];
    
    % 1. define instantaneous growth rate function for condition (FIXED)
    %fixed_trans_time = 2*60; % 2 min = 120 seconds
    %fixed_tsteps = fixed_trans_time/tstep;
    T_length = numstep;
    T_quarter = qT(t_col);
    mufun = [mu_high*ones(T_quarter,1); mu_low*ones(2*T_quarter,1); mu_high*ones(T_quarter,1)];
    mu_functions{t_col} = mufun;
    
    
    % 2. simulate!
    counter_loop = 0; % each loop is a timestep
    while gn_counter < gn+1
        
        counter_loop = counter_loop + 1;
        
        if is_newCell == 1
            % store current cell cycle Vi
            Vb_sq(gn_counter,t_col) = Vi_sq;
        end
        
        % i. calculate exponent (mu * timestep)
        expo = mufun(T_counter)*tstep_h; % mu(t)*t in hours
        
        % ii. calculate new volume after timestep
        if expo < 0
            Vt_sq = Vi_sq;
        else
            Vt_sq = Vi_sq * 2^(expo);
        end
        
        % iii. calculate time into cell cycle
        tsb = tsb + tstep_h;
        
        
        % iv. store timestep data for visualization
        V(counter_loop) = Vt_sq;
        if counter_loop == 1
            t(counter_loop) = tstep_h;
        else
            t(counter_loop) = t(counter_loop-1) + tstep_h;
        end
        
        
        % v. determine whether cell divides
        if tsb >= cur_tau
            % a. division! store V_div and tau data
            Vd_sq(gn_counter,t_col) = Vt_sq;
            tau_sq(gn_counter,t_col) = cur_tau;
            
            % b. re-set parameters for next cell cycle
            gn_counter = gn_counter + 1;
            Vi_sq = Vt_sq/2;
            cur_tau = -(17/60)*Vi_sq + (64/60);
            is_newCell = 1;
            tsb = 0;
        else
            Vi_sq = Vt_sq;
            is_newCell = 0;
        end
        
        
        % vi. determine whether nutrient period restarts
        if T_counter == T_length(t_col)
            T_counter = 1;
        else
            T_counter = T_counter + 1;
        end
        
    end
    
    figure(5)
    hold on
    plot(t,V)
    axis([0 30 0.5 3.2])
    xlabel('Time (h)')
    ylabel('Volume')
    title(strcat('visualization:',num2str(T(t_col))))
    
end
figure(5)
title('square')
legend('5','15','60')

figure(6)
histogram(Vb_sq(50:end,1),'FaceColor',rgb('CadetBlue'))
hold on
histogram(Vb_sq(50:end,2),'FaceColor',rgb('DodgerBlue'))
hold on
histogram(Vb_sq(50:end,3),'FaceColor',rgb('Indigo'))
legend('5','15','60')
xlabel('Birth volume')
ylabel('Cell count')
title('square')

% mufun_fluc{t300} = mu_functions{1};
% birthSize_fluc{t300} = Vb_sq(:,1);
% divSize_fluc{t300} = Vd_sq(:,1);
% tau_fluc{t300} = tau_sq(:,1);
% 
% mufun_fluc{t900} = mu_functions{2};
% birthSize_fluc{t900} = Vb_sq(:,2);
% divSize_fluc{t900} = Vd_sq(:,2);
% tau_fluc{t900} = tau_sq(:,2);
% 
% mufun_fluc{t3600} = mu_functions{3};
% birthSize_fluc{t3600} = Vb_sq(:,3);
% divSize_fluc{t3600} = Vd_sq(:,3);
% tau_fluc{t3600} = tau_sq(:,3);


%% Part 2C. simulate timescales with asymetric waves (varied duty cycle)

%     start with same Vi (2 cubic microns)
%     use same tau(Vi)
%     use different mu(t)
%     simulate 1000 generations

clear
clc

% 0. initialize different signals
%      i) asymetric, most high
%     ii) asymetric, somewhat high
%    iii) symetric, square wave
%     iv) asymetric, somewhat low
%      v) asymetric, most low
mu_bases{1} = [1,1,1,1,1,0];
mu_bases{2} = [1,1,1,1,0,0];
mu_bases{3} = [1,1,1,0,0,0];
mu_bases{4} = [1,1,0,0,0,0];
mu_bases{5} = [1,0,0,0,0,0];


% 0. initialize periods of interest (seconds)
%T = [30, 300, 900, 3600];
T = [3600, 900, 300, 30]; % for easier plotting


% 0. initialize timestep in simulation (seconds)
tstep = 5;
tstep_h = tstep/3600;


% 0. initialize parameters for instantaneous growth rate function 
mu_high = 1.9; % units = 1/h
mu_low = 0.7;  % units = 1/h
numstep = T/tstep; % seconds, conversion to h occurs in simulation


% 0. initialize vectors for data storage
generations = 1000;
tau_compiled = cell(length(mu_bases),length(T));
Vd_compiled = cell(length(mu_bases),length(T));
Vb_compiled = cell(length(mu_bases),length(T));
delta_compiled = cell(length(mu_bases),length(T));
mufun_compiled = cell(length(mu_bases),length(T));


for tscale = 1:length(T) % need to make this a loop for multiple timescales
    
    
    % 0. for each nutrient timescale, determine number of timesteps
    T_length = numstep(tscale);
    
    % 1. for each growth rate signal...
    for musig = 1:length(mu_bases)
        
        % 1. start generation counter
        gen_counter = 1;
        
        
        % 1. initiation size, timestep, tau
        is_newCell = 1; % 0 = division did not just happen
                        % 1 = birth timestep!
                        %     record birth size and reset current_tau
        Vi = 2;         % first birth size, fixed across conditions
        gen_tau = -(17/60)*Vi + (64/60);
        T_counter = 1; % time counter for growth rate
        tsb = 0;       % time since birth for cell cycle
        V = [];
        t = [];
        
        
        % 1. initialize vectors for current simulation data
        tau = nan(generations,1);
        Vdiv = nan(generations,1);
        Vbirth = nan(generations,1);
        
        % 1. define instantaneous growth rate function
        current_signal = mu_bases{musig};
        segment_length = T_length/length(current_signal);
        
        mufun = [];
        for ff = 1:length(current_signal)
            
            isHigh = current_signal(ff);
            if isHigh == 1
                mufun = [mufun; mu_high*ones(segment_length,1)];
            else
                mufun = [mufun; mu_low*ones(segment_length,1)];
            end
            
        end
        mufun_compiled{musig,tscale} = mufun;
        clear ff isHigh
        
        
        % 2. simulate!
        counter_loop = 0; % each loop is a timestep
        while gen_counter < generations+1
            
            counter_loop = counter_loop + 1;
            
            if is_newCell == 1
                % store current cell cycle Vi
                Vbirth(gen_counter,1) = Vi;
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
            if tsb >= gen_tau
                
                % a. division! store V_div and tau data
                Vdiv(gen_counter,1) = Vt;
                tau(gen_counter,1) = gen_tau;
                
                % b. re-set parameters for next cell cycle
                gen_counter = gen_counter + 1;
                Vi = Vt/2;
                gen_tau = -(17/60)*Vi + (64/60);
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
        
        % 3. store data for current signal of current timescale
        tau_compiled{musig,tscale} = tau;
        Vb_compiled{musig,tscale} = Vbirth;
        Vd_compiled{musig,tscale} = Vdiv;
        delta_compiled{musig,tscale} = Vdiv - Vbirth;
        clear tau Vbirth Vdiv
        clear V Vi Vt tsb t T_counter gen_tau gen_counter expo is_newCell counter_loop
        
    end
    
end

% 4. analyze simulated data from all signals and timescales
Vb_means = cellfun(@mean,Vb_compiled);
fun_means = cellfun(@mean,mufun_compiled);
delta_means = cellfun(@mean,delta_compiled);

figure(1)
for hh = 1:length(T)
    histogram(Vb_compiled{5,hh})
    hold on
end
legend('3600','900','300','30')
xlabel('birth volume')

figure(2)
for hh = 1:length(T)
    histogram(tau_compiled{5,hh})
    hold on
end
legend('3600','900','300','30')
xlabel(' interdivision time')

