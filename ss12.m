%% ss12: delta_i vs Vb_i compiled from all replicates


%  Goal: compile replicate data from fluctuating and steady conditions and
%        plot delta_i vs Vb_i as a scatter


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keep replicates apart
%  Part 2. compile & plot replicate data from each fluctuating condition
%  Part 3. compile & plot replicate data from each steady condition


%  Last edit: Jen Nguyen, 2020 May 21
%  Commit: add part 3 to plot tau_i vs delta_i for steady conditions

%  OK let's go!

%% Part 0. initialize data, saved from figure1A_division.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')


% 0. initialize plotting parameters
palette = {'DarkTurquoise','Chocolate'};
environment_order = {'low',30,300,900,3600,'ave','high'};
palette_steady = {'DarkSlateBlue','DarkGoldenrod','DarkRed','SlateBlue','Goldenrod','Crimson'};
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
clear t0_5 t5 t15 t60


%% Part 2. compile replicate data from each fluctuating condition

% 0. initialize parameters for which calculate stats in organized_data
vol_birth = 1;     % volume at birth = col in cc (compiled in figure1A_division.m)
vol_div = 2;     % volume at division = col in cc
tau = 2;    % interdivision time = col 2 in meta


% 0. for each condition of interest
conditionsOI = [900; 3600];
counter = 0;
numcells = [];
delta_i = [];
Vb_i = [];
for cond = 1:length(conditionsOI)
    
    % 1. isolate replicate data from each conditions
    currCond = conditionsOI(cond);
    if currCond == 900
        col = 4;
    elseif currCond == 3600
        col = 5;
    end
    currData = organized_data(:,col);
    noData = cellfun(@isempty,currData);
    currData = currData(noData == 0);
    
    
    % 2. compile data from each replicate 
    for rep = 1:length(currData)
        
        % i. tick for figure counter
        counter = counter + 1;
        
        % ii. gather parameters
        repData = currData{rep,1};
        Vb = repData.cc(:,vol_birth);
        Vdiv = repData.cc(:,vol_div);
        
        % iii. keep only data within 95% of mean Vb
        Vb_mean = mean(Vb);
        sigma = std(Vb);
        
        Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
        Vdiv_trim1 = Vdiv(Vb < Vb_mean + 2*sigma);
        
        Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        Vdiv_trim2 = Vdiv_trim1(Vb_trim1 > Vb_mean - 2*sigma);
        clear Vb_trim1 Vdiv_trim1 sigma
        clear repData Vdiv Vb sigma
        
        
        % iv. concatenate data
        delt = Vdiv_trim2 - Vb_trim2;
        delta_i = [delta_i; delt];
        Vb_i = [Vb_i; Vb_trim2];
        
    end
    clear Vb_trim2 Vdiv_trim2 delt dt_trim2
    
    % iv. plot
    figure(cond)
    scatter(Vb_i,delta_i,'MarkerEdgeColor',rgb(palette{cond}))
    xlabel('Vb_i (min)')
    ylabel('delta_i')
    numcells = length(delta_i);
    title(strcat('compiled-t',num2str(currCond/60),'-n',num2str(numcells)))
    
    % overlay fit
    fit = polyfit(Vb_i,delta_i,1);
    x = linspace(min(Vb_i),max(Vb_i),10);
    y = fit(1).*x + fit(2);
    hold on
    if currCond == 900
        text(min(Vb_i)-0.5,min(Vb_i),strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette{cond}))
    elseif currCond == 3600
        text(max(Vb_i)-1,max(Vb_i)+5,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette{cond}))
    end
    plot(x,y,'Color',rgb(palette{cond}))
        
end
clear currCond currData col rep cond
clear conditionsOI counter
clear Vb_mean x y tau fluc high low ave


%% Part 3. compile replicate data from each steady condition


% 0. initialize parameters for which calculate stats in organized_data
vol_birth = 1;     % volume at birth = col in cc (compiled in figure1A_division.m)
vol_div = 2;     % volume at division = col in cc
tau = 2;    % interdivision time = col 2 in meta


% 0. for each fluctuating condition of interest
conditionsOI = [900; 3600];
numcells = [];

for cond_i = 1:length(conditionsOI)
    
    
    % 1. loop through steady conditions and compile replicate data
    currCond_i = conditionsOI(cond_i);
    if currCond_i == 900
        rows_steady = 7:10;
    elseif currCond_i == 3600
        rows_steady = 11:13;
    end
    
    conditions_steady = {'low','ave','high'};
    color_counter = 0;
    
    for cond_ii = 1:length(conditions_steady)
        
        color_counter = color_counter + 1;
        delta_i = [];
        Vb_i = [];
       
        
        % 1. isolate replicate data from each steady condition
        currCond = conditions_steady{cond_ii};
        
        if strcmp(currCond,'low') == 1
            col = 1;
        elseif strcmp(currCond,'ave') == 1
            col = 6;
        elseif strcmp(currCond,'high') == 1
            col = 7;
        end
        condData = organized_data(:,col);
        clear col
        
        
        % 2. use only steady data from fluctuating experiments
        currData = condData(rows_steady);
        noData = cellfun(@isempty,currData);
        currData = currData(noData == 0);
        clear noData
        
        
        % 3. compile data from each replicate condition
        for rep = 1:length(currData)
           
            % i. gather parameters
            repData = currData{rep,1};
            Vb = repData.cc(:,vol_birth);
            Vdiv = repData.cc(:,vol_div);
            
            % ii. keep only data within 95% of mean Vb
            Vb_mean = mean(Vb);
            sigma = std(Vb);
            
            Vb_trim1 = Vb(Vb < Vb_mean + 2*sigma);
            Vdiv_trim1 = Vdiv(Vb < Vb_mean + 2*sigma);
            
            Vb_trim2 = Vb_trim1(Vb_trim1 > Vb_mean - 2*sigma);
            Vdiv_trim2 = Vdiv_trim1(Vb_trim1 > Vb_mean - 2*sigma);
            clear Vb_trim1 Vdiv_trim1 sigma
            clear repData Vdiv Vb sigma
            
            % iii. concatenate data
            delt = Vdiv_trim2 - Vb_trim2;
            delta_i = [delta_i; delt];
            Vb_i = [Vb_i; Vb_trim2];
    
        end
        clear Vb_trim2 Vdiv_trim2 delt
        
        
        % 5. plot
        figure(cond_i)
        scatter(Vb_i,delta_i,'MarkerEdgeColor',rgb(palette_steady{color_counter}))
        numcells = [numcells; length(delta_i)];
        
        % overlay fit
        fit = polyfit(Vb_i,delta_i,1);
        x = linspace(1,8,10);
        %x = linspace(min(Vb_i),max(Vb_i),10);
        y = fit(1).*x + fit(2);
        hold on
        if cond_ii == 1
            text(1.5,25,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette_steady{color_counter}))
        elseif cond_ii == 2
            text(1.5,22,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette_steady{color_counter}))
        elseif cond_ii == 3
            text(1.5,19,strcat('y=',num2str(fit(1)),'x+',num2str(fit(2))),'Color',rgb(palette_steady{color_counter}))
        end
        plot(x,y,'Color',rgb(palette_steady{color_counter}))
        
    end
    
end

figure(1) % cond_i = 1
axis([1 8 -5 30])
legend(num2str(numcells(1:3)))
title('T = 15 min')
xlabel('Vb_i')
ylabel('delta_i (min)')

figure(2) % cond_i = 2
axis([1 8 -5 30])
legend(num2str(numcells(4:6)))
title('T = 60 min')
xlabel('Vb_i')
ylabel('delta_i (min)')
    