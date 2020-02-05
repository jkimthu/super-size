%% Figure 1: size vs lambda (scatter)



%  Goal: plot population-level division size as a function of mean growth rate
%        across the cell cycle
%
%        steady-state populations have been demonstrated to follow a strong
%        positive correlation between the natural log of cell size and
%        growth rate.
%
%        do populations under fluctuations follow the same trend as
%        steady-state populations?
%
%        input: image processed data from each experimental replicate
%        output: plot, modeled after Figure 1C of Taheri et al., (2014)



%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. calculate difference from mean line
%  Part 3. bootstrap distance to assess variability of difference from mean line


%  Last edit: Jen Nguyen, 2019 Feb 5
%  Commit: calculate difference and error from steady mean line  

%  OK let's go!


%% Part 0. initialize data, saved from figure1A_division.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')
lamb = 1; % column in compiled_data (meta) that is lambda


% 0. initialize plotting parameters
palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
environment_order = {'low',30,300,900,3600,'ave','high'};
shape = 'o';


%% Part 1. sort data by nutrient condition, keeping replicates apart
%          this code is the same as figure1A_scatter.m (2020 Feb 4)

% 1. re-organize data by nutrient condition, keep replicates separated
t0_5 = 1:3; t5 = 4:6; t15 = 7:10; t60 = 11:13; % row number in data structure 
fluc = 1; low = 2; ave = 3; high = 4;          % row number in data structure 

counter = zeros(1,length(environment_order));
organized_data = [];
for exp = 1:length(compiled_data)
    
    for exp_cond = 1:4
        
        % assign to correct environment column based on condition
        if exp_cond == 1 % fluc
            
            if ismember(exp,t0_5) == 1
                envr = 2;
            elseif ismember(exp,t5) == 1
                envr = 3;
            elseif ismember(exp,t15) == 1
                envr = 4;
            elseif ismember(exp,t60) == 1
                envr = 5;
            end
            
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
            
        else % steady
            
            if exp_cond == 2 % low
                envr = 1;
            elseif exp_cond == 3
                envr = 6;
            elseif exp_cond == 4
                envr = 7;
            end
            
            counter(1,envr) = counter(1,envr) + 1;
            data = compiled_data{exp,1}{exp_cond,1};
            organized_data{counter(1,envr),envr} = data;
            
        end
        
    end
    
end
clear envr exp_cond data counter


%% Part 2. calculate difference from mean line

% 0. initialize which environments pertain to which columns in organized_data
steady = [1,6,7];
metric = 1; % 1 = volume at birth
lamb = 1;

% 1. compile all cells from each condition
for col = 1:length(environment_order)
    
    cond_data = organized_data(:,col);
    
    current_lambda = [];
    current_size = [];
    for rr = 1:length(cond_data)
        
        rr_data = cond_data{rr,1};
        if isempty(rr_data) == 1
            continue
        end
        lambda = rr_data.meta(:,lamb); % note: mu is all instananeous vals in each cell cycle
        sizes = rr_data.cc(:,metric);
        
        current_lambda = [current_lambda; lambda];
        current_size = [current_size; sizes];
    end
    compiled_lambda{1,col} = current_lambda;
    compiled_size{1,col} = current_size;
end
clear current_lambda current_size lambda sizes rr_data rr col



% 2. determine best fit line from steady conditions
compiled_lambda_mean = cellfun(@mean,compiled_lambda);
compiled_size_mean = cellfun(@mean,compiled_size);

steady_lambda = compiled_lambda_mean(steady);
steady_size = compiled_size_mean(steady);
fit = polyfit(steady_lambda,log(steady_size),1);
clear steady_lambda steady_size


% 3. plot best fit line
fit_slope = fit(1);
fit_B = fit(2);

fit_x = linspace(0,4,100);
fit_y = fit_slope*fit_x + fit_B;

figure(1)
plot(fit_x,fit_y,'Color',rgb('SlateGray'))
ylabel('ln birth volume')
xlabel('lambda')



% 4. calculate perpendicular distance from line for each cell
compiled_distance = [];
for cond = 1:length(environment_order) % loop through each condition
    
    current_x = compiled_lambda{1,cond}; % x = individual lambdas
    current_y = compiled_size{1,cond};   % y = indiviudal birth volumes (not log scale)
    current_y_log = log(current_y);      % convert birth volume to log scale
    % test(cond) = length(find(current_y_log < 0));
    % only 59 negative values in low, none in other conditions
    
    distance = NaN(length(current_x),1);
    sign = NaN(length(current_x),1);
    for cycle = 1:length(current_x)    % loop through each cell
        
        cc_x = current_x(cycle);
        cc_y = current_y_log(cycle);
        
%         figure(1) % plot (x,y) of current cell 
%         hold on
%         plot(cc_x,cc_y,'o','Color',rgb('DarkCyan'))
        
        
        % i. find perpendicular line to fit
        cc_m = (1/fit_slope)*(-1);
        cc_B = cc_y - (cc_m * cc_x);
        cc_line_x = linspace(0,4,100);
        cc_line_y = cc_m*cc_line_x + cc_B;
        
%         figure(1) % plots perpendicular line that goes through cell
%         hold on
%         plot(cc_line_x,cc_line_y,'Color',rgb('DarkCyan'))
%         axis([0 4 0 4])
%         grid on
%         pbaspect([1 1 1])
        
        
        % ii. find intercept between two lines
        intercept_x = (cc_B - fit_B)/(fit_slope - cc_m);
        intercept_y = (cc_m * intercept_x) + cc_B;
        
%         figure(1) % plot intercept between two lines (fit and cell)
%         hold on
%         plot(intercept_x,intercept_y,'o','Color',rgb('SlateGray'))
      
        
        % iii. calculate distance between cell and fit
        distance(cycle) = sqrt( (intercept_x-cc_x)^2 + (intercept_y-cc_y)^2 );
        
        
        % iv. record sign
        if intercept_x < cc_x
            sign(cycle) = -1;
        else
            sign(cycle) = 1;
        end
        
    end
    clear intercept_x intercept_y cc_B cc_m
    clear cc_y cc_x current_x current_y current_y_log
    
    dd = distance.*sign;
    compiled_distance{1,cond} = dd;
    
end
clear distance sign dd cycle


% 5. calculate mean distance and standard error for each condition
d_mean = cellfun(@mean,compiled_distance);
d_std = cellfun(@std,compiled_distance);
d_count = cellfun(@length,compiled_distance);
d_sem = d_std./sqrt(d_count);


% 6. plot mean distance from line and error
figure(2)
bar(1:7,d_mean)              
hold on
%er = errorbar(1:7,d_mean,d_sem,d_sem);    
er = errorbar(1:7,d_mean,d_std,d_std);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
fname = {'low','30 s','5 min','15 min','60 min','ave','high'}; 
set(gca, 'XTick', 1:length(fname),'XTickLabel',fname);
ylabel('distance')

%legend('low','30','300','900','3600','ave','high')
clear d_* er fname cond


%% Part 3. bootstrap distance to assess variability of mean (x,y)

% for a given fluctuating condition, sample 1000 populations with
% replacement from cell cycle data (fixing pairs) for each true population size

palette = {'Indigo','DarkTurquoise','SteelBlue','DeepSkyBlue','DodgerBlue','GoldenRod','FireBrick'};
counter_figs = 0;
%fluctuating = [2,3,4,5];

for condition = 1:length(compiled_distance)
    
    % 1. isolate measured distances
    curr_distance = compiled_distance{:,condition};
    
    % 2. simulate 1000 populations of various size, sampling with replacement
    popsize = [100, 1000, 10000];
    npops = 1000;
    for pp = 1:length(popsize)
        
        sim_means = NaN(1000,1);
        for ii = 1:npops
            
            sim_dist = datasample(curr_distance,popsize(pp));
            sim_means(ii) = mean(sim_dist);
            
        end
        compiled_simulation_means{pp,condition} = sim_means;

        figure(condition)
        subplot(1,length(popsize),pp)
        histogram(sim_means,'FaceColor',rgb(palette(condition)))
        hold on
        title(strcat('mean:',num2str(mean(sim_means)),'; std:',num2str(std(sim_means))))
    end
end


%legend('low','30','300','900','3600','ave','high')
