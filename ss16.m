%% ss16: single cell tau_i vs Vb_i in steady environments


%  Goal: determine whether distributions of slope (tau vs Vb) vary between
%        single cells and independent replicates.

%        Plot distributions of slope from each steady condition (low, ave, high)



%  Strategy: 
%
%  Part 0. initialize analysis
%  Part 1. sort data by nutrient condition, keeping replicates apart
%  Part 2. single cell distributions, loop through steady conditions and:
%            - calculate single-cell slopes (tau_i divided by Vb_i)
%            - plot distribution of slopes, draw line at mean
%  Part 3. replicate distributions
%            - plot distribution of slopes measured from each replicate,
%              draw a line at mean



%  Last edit: Jen Nguyen, 2020 July 28
%  Commit: first commit, check if distributions of slope differ between
%          single cell and replicate average


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


%% Part 2. distributions of single-cell tau_i and Vb_i

%            - calculate single-cell slopes (tau_i divided by Vb_i) from all
%              lineages that have at least 3 cell cycles
%            - plot distribution of slopes, draw line at mean

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss16')

% for each steady condition
steadies = find(cellfun(@isstr,environment_order));
slope_i_means = nan(13,3);
slope_i_stds = nan(13,3);
slope_i_all = cell(13,3);
Vb_i = cell(13,3);
tau_i = cell(13,3);
slopes_ave = nan(13,3);
yints_ave = nan(13,3);
for sc = 1:3 % sc = steady condition
    
    % 1. isolate data from current steady condition
    col = steadies(sc);
    palette = palette_steady(1,sc);
    current_condData = organized_data(:,col);
    
    % 2. determine replicates to loop through
    numreps = length(current_condData);
    for rr = 1:numreps
        
        if isempty(current_condData{rr,1}) == 1
            continue
        else
            % 3. isolate replicate data
            rep_meta = current_condData{rr,1}.meta;
            rep_taus = rep_meta(:,2);       % each row = cc; number = tau (min)
            rep_cc = current_condData{rr,1}.cc;
            rep_birthvols = rep_cc(:,1);    % each row = cc; number = birth volume (cubic um)
            clear rep_meta rep_cc
            
            
            % 4. calculate single-cell slopes
            rep_slopes = rep_taus./rep_birthvols;
            
            
            % 5. remove items over 3 std from mean
            mm = mean(rep_slopes);
            thresh = 3 * std(rep_slopes);
            upper = mm + thresh;
            lower = mm - thresh;
            cut_up = rep_slopes > upper;
            cut_down = rep_slopes < lower;
            toCut = cut_up + cut_down;
            
            slopes_final = rep_slopes(toCut == 0);
            taus_final = rep_taus(toCut == 0);
            Vb_final = rep_birthvols(toCut == 0);
            clear mm thresh upper lower cut_up cut_down toCut
            clear rep_taus rep_birthvols rep_slopes
            
            
            % 6. plot distribution of single-cell slopes replicate
            figure(rr)
            histogram(slopes_final,100,'FaceColor',rgb(palette));
            xlabel('single-cell slope')
            ylabel('# cells')
            title(strcat('replicate',num2str(rr)))
            hold on
            
            
            % 7. store all means and standard deviations
            mf = mean(slopes_final);
            sf = std(slopes_final);
            slope_i_means(rr,sc) = mf;
            slope_i_stds(rr,sc) = sf;
            slope_i_all{rr,sc} = slopes_final;
            Vb_i{rr,sc} = Vb_final;
            tau_i{rr,sc} = taus_final;
            clear mf sf
            
            
            % 8. plot scatter of slope_i
            figure(rr+20)
            scatter(Vb_final,taus_final,100,rgb(palette))
            
            % lineage fit
            fit = polyfit(Vb_final,taus_final,1);
            xx = linspace(0,9,10);
            yy = fit(1).*xx + fit(2);
            hold on
            plot(xx,yy,'Color',rgb(palette),'LineWidth',2)
            axis([0 9 0 90])
            ylabel('tau_i')
            xlabel('Vb_i')
            title(strcat('replicate',num2str(rr)))
            
            slopes_ave(rr,sc) = fit(1);
            yints_ave(rr,sc) = fit(2);
            clear xx yy fit Vb_final taus_final
            
        end
                
    end

end
clear sc rr

slopes_i.vals = slope_i_all;
slopes_i.means = slope_i_means;
slopes_i.stds = slope_i_stds;

save('ss16_output.mat','slopes_ave','yints_ave','slopes_i','tau_i','Vb_i')


for rep = 1:numreps
    
    % 8. save and close plot (prepare to loop through next condition)
    figure(rep)
    saveas(gcf,strcat('ss16-i-rep-',num2str(rep)),'epsc')
    close(gcf)
    
    figure(rep+20)
    saveas(gcf,strcat('ss16-ave-rep-',num2str(rep)),'epsc')
    close(gcf)
    
end

clear col sc slope_i_all slope_i_means slope_i_stds


%% Part 3. replicate distributions
%          run Part 0 and 1, then load items from Part 2.

cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss16')
load('ss16_output.mat')

for ss = 1:3 % sc = steady condition
    
    % 1. isolate data from current steady condition
    palette = palette_steady(1,ss);
    slopes_cond = slopes_ave(:,ss);

    % 2. plot distribution of slopes per condition
    figure(1)
    histogram(slopes_cond,'FaceColor',rgb(palette))
    hold on
    
end
clear ss

figure(1)
axis([-22 0 0 8])
xlabel('average slope')
ylabel('# replicates')






