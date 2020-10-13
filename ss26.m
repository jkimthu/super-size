%% ss26: using residuals compare linear vs hyperbolic fit single condition data

%  Goal: which is a better fit for a given nutrient condition?


%  Strategy:
%
%       For a given experimental condition:
%        1.  Plot scatter
%        2.  Fit line and fit hyperbola
%        3.  Calculate residuals R to determine which gives less error


%  Last edit: jen, 2020 Oct 12
%  Commit: first commit, residuals for linear and hyperbolic fit to single nutrient conditions


%  OK let's go!

%% Part 0. initialize cell cycle data, saved from ss17.m

clear
clc

% 0. initialize complete meta data
cd('/Users/jen/super-size/')
load('ss17.mat')

clear fluc_lambda_i


%% Part 1. restructure steady tau_i, slope_i and Vb_i into columns by condition
%  adapted from Part 2 of ss24.m


% 1. initialize data vectors for slope_i and Vb_i
numcells_steady = [];
steady_tau_i = [];
steady_slope_i = [];
steady_Vb_i = [];


% 2. loop through timescale groups to concatenate data into columns by
%    nutrient conditiion
for ii = 1:length(compiled_steady)
    
    num_ii = cellfun(@length,compiled_steady{1,ii}.tau_i);
    numcells_steady = [numcells_steady; num_ii];
    
    tau_ii = compiled_steady{1,ii}.tau_i;
    steady_tau_i = [steady_tau_i; tau_ii];
    
    vb_ii = compiled_steady{1,ii}.Vb_i;
    steady_Vb_i = [steady_Vb_i; vb_ii];
    
    slope_ii = compiled_steady{1,ii}.slope_i;
    steady_slope_i = [steady_slope_i; slope_ii];
    
end
numcells_steady(numcells_steady == 0) = NaN;
clear ii num_ii vb_ii tau_ii slope_ii


%% Part 2. restructure data for replicate based plotting

%  re-organize data so that fluc and steady data are in the same matrix,
%  such that:

%  column order denotes nutrient condition: fluc, low, ave, high
%  row order denotes replicate of fluc timescale: 1-4 (30 s), 5-8 (5 min), 9-12 (15 min), 13-16(60 min)

%  strategy:
%  1. restructure fluc data from 4x4 to 16x1
%  2. concatenate 16x3 steady data such that fluc data is 1st column

%  OK let's go!

%  1. restructure fluc data from 4x4 to 16x1

fluc_Vb_i_restruct = [];
fluc_tau_i_restruct = [];
fluc_slope_i_restruct = [];

for col = 1:4
    fluc_Vb_i_restruct = [fluc_Vb_i_restruct; fluc_Vb_i(:,col)];
    fluc_tau_i_restruct = [fluc_tau_i_restruct; fluc_tau_i(:,col)];
    fluc_slope_i_restruct = [fluc_slope_i_restruct; fluc_slopes_i(:,col)];
end


%  2. concatenate 16x3 steady data such that fluc data is 1st column
Vb_i_data = [fluc_Vb_i_restruct, steady_Vb_i];
tau_i_data = [fluc_tau_i_restruct, steady_tau_i];
slope_i_data = [fluc_slope_i_restruct, steady_slope_i];


clear col fluc_Vb_i fluc_tau_i fluc_slope_i steady_Vb_i steady_tau_i steady_slope_i
clear fluc_Vb_i_restruct fluc_tau_i_restruct fluc_slope_i_restruct


%% Part 3. calculate residuals R for linear and hyperbolic fits

% 0. initialize num of conditions
palette = {'DodgerBlue','Indigo','GoldenRod','FireBrick'};
ss = size(tau_i_data);
rows = ss(1); % replicate experiments
cols = ss(2); % nutrient conditions (fluc, low, ave, high)
clear ss


% 1. loop through each experiment
residuals_lin_slope = nan(rows,cols);
residuals_lin_tau = nan(rows,cols);
residuals_hypr_slope = nan(rows,cols);
residuals_hypr_tau = nan(rows,cols);

fit_lin_slope = cell(rows,cols);
fit_lin_tau = cell(rows,cols);
fit_hypr_slope = cell(rows,cols);
fit_hypr_tau = cell(rows,cols);
    
counter = 0;
for ee = 1:rows %
    
    ee_taus = tau_i_data(ee,:);
    ee_vbs = Vb_i_data(ee,:);
    ee_slopes = slope_i_data(ee,:);

    
    % 2. for each condition...
    for cc = 1:length(ee_vbs)
        
        % 3. isolate condition tau_i, slope_i & Vb_i
        cc_taus = ee_taus{cc};
        if isempty(cc_taus) == 1
            continue
        end
        counter = counter + 1;
        cc_vbs = ee_vbs{cc};
        cc_slopes = ee_slopes{cc};
        cc_color = rgb(palette{cc});
        
        
        % 4. plot tau_i vs. Vb_i and slope_i vs. Vb_i
%         figure(1)
%         scatter(cc_vbs,cc_taus,'MarkerFaceColor',cc_color,'MarkerEdgeColor',cc_color)
%         hold on
%         
%         figure(2)
%         scatter(cc_vbs,cc_slopes,'MarkerFaceColor',cc_color,'MarkerEdgeColor',cc_color)
%         hold on
        
        
        % 5. fit hyperbola to aggregate single-cell data
        %    hyperbola = @(b,x) b(1) + b(2)./(x);
        %    residual norm cost func = @(b) norm(y - hyprb(b,x));
        %    initial parameters, B0 = [1,1]
        %    estimated parameters, B = [b(1), b(2)];
        y_tau = cc_taus;
        y_slope = cc_slopes;
        x_vb = cc_vbs;
        B0 = [1; 1];
        
        % tau = y parameter
        hyprb = @(b,x_vb) b(1) + b(2)./(x_vb);
        rncf = @(b) norm(y_tau - hyprb(b,x_vb));
        B_tau = fminsearch(rncf, B0);
        
        % slope = y parameter
        hyprb = @(c,x_vb) c(1) + c(2)./(x_vb);
        rncf = @(c) norm(y_slope - hyprb(c,x_vb));
        B_slope = fminsearch(rncf, B0);
        
        
        % 6. plot hyperbola over scattered data
%         x = linspace(1,9,20);
%         y = B_tau(1) + B_tau(2)./x;
%         figure(1)
%         plot(x,y,'Color',cc_color,'LineWidth',3)
%         
%         y = B_slope(1) + B_slope(2)./x;
%         figure(2)
%         plot(x,y,'Color',cc_color,'LineWidth',3)
        
        
        %7. fit and plot line over scattered data
        linear_slope = polyfit(x_vb,y_slope,1);
        linear_tau = polyfit(x_vb,y_tau,1);
        
%         y = linear_tau(2) + linear_tau(1).*x;
%         figure(1)
%         plot(x,y,'Color',cc_color,'LineWidth',1)
%         
%         y = linear_slope(2) + linear_slope(1).*x;
%         figure(2)
%         plot(x,y,'Color',cc_color,'LineWidth',1)
        
        
        % 8. store fits
        fit_lin_slope{ee,cc} = linear_slope;
        fit_lin_tau{ee,cc} = linear_tau;
        fit_hypr_slope{ee,cc} = B_slope;
        fit_hypr_tau{ee,cc} = B_tau;
        
        
        % 9. calculate residuals for each condition
        %    Residual = sum( (y_i,model - y_i,observed)^2 )
        
        % i. calculate y_i,model
        y_model_linear_slope = (linear_slope(1).*x_vb) + linear_slope(2);
        y_model_linear_tau = (linear_tau(1).*x_vb) + linear_tau(2);
        y_model_hypr_slope = B_slope(1) + B_slope(2)./x_vb;
        y_model_hypr_tau = B_tau(1) + B_tau(2)./x_vb;
        

        % ii. calculate difference between y_i,model and y_i,observed
        diff_linear_slope = y_model_linear_slope - y_slope;
        diff_linear_tau = y_model_linear_tau - y_tau;
        diff_hypr_slope = y_model_hypr_slope - y_slope;
        diff_hypr_tau = y_model_hypr_tau - y_tau;
        clear y_model*
        
        % iii. square and sum difference
        R_linear_slope = sum(diff_linear_slope.^2);
        R_linear_tau = sum(diff_linear_tau.^2);
        R_hypr_slope = sum(diff_hypr_slope.^2);
        R_hypr_tau = sum(diff_hypr_tau.^2);
        clear diff_*
        
        % 10. store residuals for each condition
        residuals_lin_slope(ee,cc) = R_linear_slope;
        residuals_lin_tau(ee,cc) = R_linear_tau;
        residuals_hypr_slope(ee,cc) = R_hypr_slope;
        residuals_hypr_tau(ee,cc) = R_hypr_tau;
        clear R_linear* R_hypr*

    end
    clear cc_slopes cc_taus cc_vbs cc_color cc x x_vb y y_slope y_tau
    clear rncf linear_slope linear_tau hyprb B0 B_slope B_tau
    clear ee_*

    
%     cd('/Users/jen/Documents/StockerLab/Writing/manuscript 2/superSize_figs/ss26/')
%     
%     figure(1)
%     xlabel('Vb_i (cubic um)')
%     ylabel('tau_i (min)')
%     title('hyperbolic fit')
%     saveas(gcf,strcat('ss26-condition-fit-ytau-rep',num2str(counter)),'epsc')
%     close(gcf)
%     
%     figure(2)
%     xlabel('Vb_i (cubic um)')
%     ylabel('slope_i (min/cubic um)')
%     title('hyperbolic fit')
%     saveas(gcf,strcat('ss26-condition-fit-yslope-rep',num2str(counter)),'epsc')
%     close(gcf)

    
end
clear ee counter


%% Part 4. visualize residuals for all experiments in bar plot


% 0. initialize matrices for data reformatting
r_lin_slope = nan(rows,7);
r_lin_tau = nan(rows,7);
r_hy_slope = nan(rows,7);
r_hy_tau = nan(rows,7);


% 1. organize data for easy plotting
fluc = 1; low = 2; ave = 3; high = 4;
t30 = 1:3; t300 = 5:7; t900 = 9:12; t3600 = 13:15;


% i. steady data        
r_lin_slope(:,1:3) = residuals_lin_slope(:,low:high);
r_lin_tau(:,1:3) = residuals_lin_tau(:,low:high);
r_hy_slope(:,1:3) = residuals_hypr_slope(:,low:high);
r_hy_tau(:,1:3) = residuals_hypr_tau(:,low:high);

% ii. fluc data
r_lin_slope(1:3,4) = residuals_lin_slope(t30,fluc);
r_lin_tau(1:3,4) = residuals_lin_tau(t30,fluc);
r_hy_slope(1:3,4) = residuals_hypr_slope(t30,fluc);
r_hy_tau(1:3,4) = residuals_hypr_slope(t30,fluc);

r_lin_slope(1:3,5) = residuals_lin_slope(t300,fluc);
r_lin_tau(1:3,5) = residuals_lin_tau(t300,fluc);              
r_hy_slope(1:3,5) = residuals_hypr_slope(t300,fluc);
r_hy_tau(1:3,5) = residuals_hypr_tau(t300,fluc);  

r_lin_slope(1:4,6) = residuals_lin_slope(t900,fluc);
r_lin_tau(1:4,6) = residuals_lin_tau(t900,fluc);              
r_hy_slope(1:4,6) = residuals_hypr_slope(t900,fluc);
r_hy_tau(1:4,6) = residuals_hypr_tau(t900,fluc); 

r_lin_slope(1:3,7) = residuals_lin_slope(t3600,fluc);
r_lin_tau(1:3,7) = residuals_lin_tau(t3600,fluc);              
r_hy_slope(1:3,7) = residuals_hypr_slope(t3600,fluc);
r_hy_tau(1:3,7) = residuals_hypr_tau(t3600,fluc); 
               

% 2. re-organize further: slopes and taus
%    linear, then hyperbolic

% i. slopes
slopes_compiled = [r_lin_slope(:,1), r_hy_slope(:,1), r_lin_slope(:,2), r_hy_slope(:,2), r_lin_slope(:,3), r_hy_slope(:,3), r_lin_slope(:,4), r_hy_slope(:,4), r_lin_slope(:,5), r_hy_slope(:,5), r_lin_slope(:,6), r_hy_slope(:,6), r_lin_slope(:,7), r_hy_slope(:,7)];
taus_compiled = [r_lin_tau(:,1), r_hy_tau(:,1), r_lin_tau(:,2), r_hy_tau(:,2), r_lin_tau(:,3), r_hy_tau(:,3), r_lin_tau(:,4), r_hy_tau(:,4), r_lin_tau(:,5), r_hy_tau(:,5), r_lin_tau(:,6), r_hy_tau(:,6), r_lin_tau(:,7), r_hy_tau(:,7)];


% 3. box plot comparing linear and hyperbolic fit per condition
figure(1)
boxplot(slopes_compiled)
ylabel('R')
xlabel('Condition (lin vs hypr)')
title('residuals for y=slope')

figure(2)
boxplot(taus_compiled)
ylabel('R')
xlabel('Condition (lin vs hypr)')
title('residuals for y=tau')

     
