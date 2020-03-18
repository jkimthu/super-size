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

  i. size_metrics
     This script loops through data files from each experimental replicate to compile a matrix of cell size measurements, organized by cell cycle.


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


## Figures and data for Jonasz

### ss2 
1. data vectors of mu(t) from T = 15 and 60 min.
- mu(t) is the mean of all mu_i calculated from all cells from all replicate conditions.
2. MATLAB figure of mu(t) from T = 15 and 60 min plotted over the nutrient period (high:low:low:high)

### ss7
1. seven plots (4 for T = 15 min and 3 for T = 60 min) of tau_i vs Vb_i.
- n = number of individial cells (division to division event) that contribute to the scatter
- subplot 1 is tau_i vs lambda_i
- subplot 2 is tau_i vs Vb_i
- linear fit to scatter is displayed for A and B

