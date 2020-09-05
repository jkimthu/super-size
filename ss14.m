%% ss14: single cell tau_i vs Vb_i in steady environments


%  Goal: determine whether single cells display slope between tau_i and Vb_i.
%        plot single-cell tau_i and Vb_i over time in steady environments:
%         1. steady low
%         2. steady ave
%         3. steady high


%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keeping replicates apart
%  Part 2. loop through steady conditions and plot single-cell tau_i and Vb_i
%          of lineages with most generations
%  Part 3. distributions of single-lineage slopes per condition



%  Last edit: Jen Nguyen, 2020 June 3
%  Commit: plot tau_i and Vb_i for single-cell trajectories


%  OK let's go!

%% Part 0. initialize analysis

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('storedMetaData.mat')
load('A1_div_ccSize.mat')


% 0. initialize plotting parameters
environment_order = {'low',30,300,900,3600,'ave','high'};
%palette_fluc = {'DarkOrchid','DodgerBlue','Chocolate','LimeGreen'};
palette_steady = {
                  'DarkSlateBlue', 'Chocolate',     'DarkRed';
                  'SlateBlue',     'Peru',          'LightCoral';
                  'DarkMagenta',   'DarkGoldenrod', 'IndianRed';
                  'DarkOrchid',    'Goldenrod',     'FireBrick';
                  'Indigo',        'Sienna',        'Crimson';
                  'Amethyst',      'Tan',           'Coral';
                  'Violet',        'Gold',          'Tomato';
                  'MediumPurple', 'DarkKhaki',     'Maroon';
                  'Orchid',        'Khaki',         'Red';
                  'Plum',          'PaleGoldenRod', 'Salmon';
                  };
%shape = 'o';


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
clear compiled_data low ave high fluc exp


%% Part 2. loop through steady conditions and plot single-cell tau_i and Vb_i

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss14')

% for each steady condition
steadies = find(cellfun(@isstr,environment_order)); 
mingens{1} = [6; 7; 7; 6; 9; 8; 2; NaN; 7; 6; 6; 6; 7]; % low
mingens{2} = [9; 8; 15; 11; 14; 3; 4; 6; 6; 5; 11; 10; 10]; % ave
mingens{3} = [7; 8; 19; 9; 18; NaN; 5; NaN; 17; 15; 18; 9; 17]; % high
for sc = 1:3
    
    % 1. isolate data from current steady condition
    col = steadies(sc);
    palette = palette_steady(:,sc);
    current_condData = organized_data(:,col);
    
    % 2. from each replicate, determine single cell lineages with most generations
    %    note all data from A1_div_ccSize.mat is already trimmed to post 3h
    numreps = length(current_condData);
    pop_fits = nan(numreps,2);
    
    % determine shortest length of replicates to plot 
%     for rr = 1:numreps
%         if isempty(current_condData{rr,1}) == 1
%             continue
%         else
%             rep_meta = current_condData{rr,1}.meta;
%             rep_lineages = rep_meta(:,3);   % each row = cc; number = lineage (track num)
%             clear rep_meta rep_cc
%             
%             % ii. determine lineages with greatest number of generations
%             [unique_lineages,~,idx] = unique(rep_lineages);
%             numgens = accumarray(idx(:),1);
%             figure(rr)
%             histogram(numgens)     % during dev, plot histogram of lineage lengths
%             clear idx              % manually select number of generations for each rep
%         end
%     end
    
    for rr = 1:numreps
        
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            % i. isolate rep data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_lineages = rep_meta(:,3);   % each row = cc; number = lineage (track num)
            %rep_birthtimes = rep_meta(:,4); % each row = cc; number = time of birth (h)
            %rep_divtimes = rep_meta(:,5);   % each row = cc; number = time of division (h)
            
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % ii. determine lineages with greatest number of generations
            [unique_lineages,~,idx] = unique(rep_lineages);
            numgens = accumarray(idx(:),1);
            clear idx
            
            % iii. plot population data of current replicate
            figure(rr)
            scatter(rep_birthvols,rep_taus,100,rgb('Silver'))
            
            % overlay population level fit and store slope and intercept data
            fit_p = polyfit(rep_birthvols,rep_taus,1);
            x = linspace(1,8,10);
            y = fit_p(1).*x + fit_p(2);
            hold on
            plot(x,y,'Color',rgb('Silver'))
            
            pop_fits(rr,1) = fit_p(1);
            pop_fits(rr,2) = fit_p(2);
            clear x y fit_p
            
            
            % iv. plot tau_i vs Vb_i from each of longest lineages
            pl_gennum = mingens{sc}(rr);
            pl_lineages = unique_lineages(numgens >= pl_gennum,1);
            if length(pl_lineages) > 10
                pl_lineages = pl_lineages(1:10,1);
            end
            clear pl_gennum
            
            for ln = 1:length(pl_lineages)
                
                curr_linnum = pl_lineages(ln);
                curr_rows = find(rep_lineages == curr_linnum);
                curr_color = rgb(palette{ln});
                
                tau_i = rep_taus(curr_rows);
                Vb_i = rep_birthvols(curr_rows);
                
                figure(rr)
                hold on
                scatter(Vb_i,tau_i,100,'filled','MarkerFaceColor',curr_color)
                
                % lineage fit
                curr_fit = polyfit(Vb_i,tau_i,1);
                xl = linspace(1,8,10);
                yl = curr_fit(1).*xl + curr_fit(2);
                hold on
                plot(xl,yl,'Color',curr_color)
                
                ln_fits(ln,:) = curr_fit;
                clear xl yl curr_fit curr_color tau_i Vb_i 
                
            end
            clear rep_birthvols rep_lineages rep_taus
            
            % 3. save and close plot (prepare to loop through next condition)
            figure(rr)
            saveas(gcf,strcat('ss14-steady-',environment_order{col},'-rep-',num2str(rr)),'epsc')
            close(gcf)
            
            % 4. store slope data from fits
            steady_lineage_fits{rr,sc} = ln_fits;
            clear ln ln_fits curr_linnum curr_rows pl_lineages

        end
    end
    
    steady_population_fits{1,sc} = pop_fits;
    clear pop_fits rr
    
end
clear col sc exp


%% Part 3.








