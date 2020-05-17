# super-size
Scripts plotting the effect of nutrient fluctuations on cell size




## Figures:

1A. mean birth size vs mean growth rate of each nutrient condition
    size metrics: birth volume of cell cycle, birth volume of next cell cycle
                  birth length of cell cycle, birth length of next cell cycle
                  birth width of cell cycle, birth width of next cell cycle
                  birth SA:V ratio of cell cycle, birth SA:V ratio of next cell cycle



## Supplementary Figures:





## Data processing and organization:

  All figures use input data compiled into a data structure by the following scripts:

  i. size_metrics (to be organized from figure1A_division.m)
     This script loops through data files from each experimental replicate to compile a matrix of cell size measurements, organized by cell cycle.

     currently, figure1A_division.m saves a file " A1_ccSize.mat "
     this includes a data structure "compiled_data" which includes 4 matrices:
        1) cc: 8 parameters (col)
           each row, values of each parameter for current cell cycle (cc)
        		1. Vbirth
        		2. Vdiv
        		3. Lbirth
        		4. Ldiv
        		5. Wbirth
        		6. Wdiv
        		7. SA2Vbirth
        		8. SA2Vdiv 

        2) meta_data: 5 parameters (col)
           each row, values of each parameter for current cell cycle (cc)
        		1. lambdas
        		2. tau_all/60
        		3. tracks_all
        		4. Tbirth_all
        		5. Tdiv_all
                
        3) plus1_data: 4 parameters (col)
           each row, values of each parameter at birth of next cell cycle (cc + 1)
        		1. plus1_vol
        		2. plus1_len
        		3. plus1_width
        		4. plus1_sa2v
                      
        4) mus_insta: each row = cell array of mu values for one individual e. coli



## Figures to figure stuff out

saved as ss...(#).m
1. 
2. simulate growth at different fluctuation timescales with measured mu(t)
3. probability of division over time in single shift environments
4. individual volume trajectories starting with cell experiencing a single nutrient upshift
5. individual volume trajectories starting with cell experiencing a single nutrient downshift
6. 
7. individual interdivision time vs individual growth rate (lambda) or birth volume
8. tau_i vs Vb_i for each replicate of T = 15 and 60 min, same plot with fit slopes
9. predicting V_mean from linear fit of tau_i vs Vb_i
10. tau_i vs Vb_i with data compiled from all replicates for each condition T = 15 and 60 min



## Code for on-going analyses

### ss1
1. histograms of Vb_i for each experimental replicate
2. histograms of tau for each experimental replicate

### ss2 
1. data vectors of mu(t) from T = 15 and 60 min.
- mu(t) is the mean of all mu_i calculated from all cells from all replicate conditions.
2. MATLAB figure of mu(t) from T = 15 and 60 min plotted over the nutrient period (high:low:low:high)

### ss3
1. probability of division in single shift experiments

### ss4
1. volume over time of single cell cycles that experience a single nutrient upshift

### ss5
1. single lineage volume tracks starting from cell cycle experiencing a single upshift

### ss6
1. probability of division as a function of nutrient score

### ss7
1. seven plots (4 for T = 15 min and 3 for T = 60 min) of tau_i vs Vb_i.
- n = number of individial cells (division to division event) that contribute to the scatter
- subplot 1 is tau_i vs lambda_i
- subplot 2 is tau_i vs Vb_i
- linear fit to scatter is displayed for A and B

### ss8
1. scatter of tau_i and birth size_i of each replicate overlaid
Q: are linear fit slopes different between conditions?
A: yes, T=15 min slopes are steeper! (slighter slopes suggest larger sizes)

### ss9
1. does tau_i(Vb_i) predict Vb_mean?
Q: can we predict V_mean from mu_mean and linear fit of tau_i and Vb_i? 
A: pretty ok, but no clear difference between timescales. suggests a transition?

### ss10
1. tau_i vs Vb_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)
scatter of tau_i and Vb_i 

### ss11
1. tau_i vs delta_i compiled from all replicates
two plots (one for each condition, T = 15 and 60 min)
scatter of tau_i and delta_i 
delta_i = Vdiv_i - Vb_i
