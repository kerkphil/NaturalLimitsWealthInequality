%-------------------------------------------------------------------------%
% Perform balanced growth simulation experiments for simple Piketty model
%-------------------------------------------------------------------------%
% This program calls the following .m and .mat file(s)
%   * kappa_ss.m: solves for the steady-state kappa_bar and k1tilbar
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Opening Commands
%-------------------------------------------------------------------------%
clear ;
%cd('/Users/rwe2/Documents/BYU Economics/Paper Ideas/simple-piketty/MATLAB') ;

starttime = clock ; % start timer

%-------------------------------------------------------------------------%
% Set parameters initial values
%-------------------------------------------------------------------------%
% A        = scalar, production function scale parameter
% alpha    = scalar, capital share of income, production function parameter
% delta_an = scalar, annual rate of capital depreciation
% delta    = scalar, generational rate of capital depreciation
% theta    = scalar, percent of population that is type 1
% g_an     = scalar, annual growth rate of labor augmenting technological
%            change
% g        = scalar, generational growth rate of labor augmenting
%            technological change
% s1       = scalar, savings rate of type 1 individuals
% s2       = scalar, savings rate of type 2 individuals
% gam1     = scalar, level parameter for utility of bequests for type 1
% gam2     = scalar, level parameter for utility of bequests for type 2
% l1       = scalar, labor endowment of type 1 households
% l2       = scalar, labor endowment of type 2 households
% L        = scalar, aggregate labor
% Delta    = scalar, interest rate wedge: Delta = r1t - r2t
%-------------------------------------------------------------------------%
A = 1 ;
alpha = 0.35 ;
delta_an = 0.08 ;
delta = 1 - (1 - delta_an)^30 ;
theta = 0.20 ;
g_an = 0.02 ;
g = (1 + g_an)^30 - 1 ;
s1 = 0.236 ;
s2 = 0.072 ;
gam1 = s1/(1-s1) ;
gam2 = s2/(1-s2) ;
l1 = 1 ;
l2 = 1 ;
L = theta*l1 + (1-theta)*l2 ;
Delta = 0 ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of (s1,s2)
%-------------------------------------------------------------------------%
% smat      = sims x 2 matrix, (s1,s2) pairs to test for kappa_bar and
%             omega_bar
% omkapmat1 = sims x 2 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to (s1,s2) combinations
% sims1     = integer, number of combinations of (s1,s2) for which to
%             compute kappa_bar and omega_bar
%-------------------------------------------------------------------------%
smat = [0.99 0.01; 0.75 0.1; 0.75 0.2; 0.5 0.3; 0.5 0.49; 0.236 0.072] ;
omkapmat1 = zeros(size(smat)) ;
sims1 = size(smat) ;
sims1 = sims1(1) ;

for sim = 1:sims1
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar by simulating a time series
    %---------------------------------------------------------------------%
    % k1init   = scalar, initial value of type 1 household wealth
    % k2init   = scalar, initial value of type 2 household wealth
    % kapinit  = scalar, initial value of kappa
    % k1       = scalar, current period's value of type 1 capital
    % k2       = scalar, current period's value of type 2 capital
    % kap      = scalar, current period's value of kappa_t
    % per      = integer, iteration number of simulation
    % kdist    = scalar, distance measure of kappa' from kappa
    % mindist  = scalar, tolerance level of distance governing convergence
    % maxper   = integer, maximum number of periods to simulate
    % K        = scalar, current period value of aggregate capital K_t
    % Y        = scalar, current period value of aggregate output Y_t
    % w        = scalar, current period value of real wage w_t
    % r        = scalar, current period value of real rental rate r_t
    % r1       = scalar, current period value of real return for type 1 hh
    % r2       = scalar, current period value of real return for type 2 hh
    % k1pr     = scalar, next period's value of type 1 capital
    % k2pr     = scalar, next period's value of type 2 capital
    % kappr    = scalar, next period's value of kappa'
    %---------------------------------------------------------------------%
    s1 = smat(sim,1) ;
    s2 = smat(sim,2) ;
    k1init = 100 ;
    k2init = 5 ;
    kapinit = k2init/k1init ;
    k1 = k1init ;
    k2 = k2init ;
    kap = kapinit ;
    per = 1 ;
    kdist = 10 ;
    mindist = 10^(-10) ;
    maxper = 50 ;

    while kdist > mindist && per < maxper
        per = per + 1 ;
        K = theta*k1 + (1 - theta)*k2 ;
        Y = A*(K^alpha)*((exp(g*per)*L)^(1 - alpha)) ;
        w = (1 - alpha)*(Y/L) ;
        r = alpha*(Y/K) - delta ;
        r1 = r + (Delta*(1-theta)*k2)/K ;
        r2 = r - (Delta*theta*k1)/K ;
        k1pr = s1*(w*l1 + (1 + r1)*k1) ;
        k2pr = s2*(w*l2 + (1 + r2)*k2) ;
        kappr = k2pr/k1pr ;
        kdist = abs(kappr - kap) ;
        k1 = k1pr ;
        k2 = k2pr ;
        kap = kappr ;
    end
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar using nonlinear system of
    % stationary steady-state equations
    %---------------------------------------------------------------------%
    % kapkinit = 1 x 2 vector, initial values for kappa_bar and k1tilbar
    % params   = 1 x 10 vector, parameters to be passed in to fsolve
    % options  = structural array, options for fsolve
    % kapk     = 1 x 2 vector, solution for kappa_bar and k1tilbar
    % kappa    = scalar, solution for kappa_bar
    % k1tilbar = scalar, solution for k1tilbar
    %---------------------------------------------------------------------%
    kapkinit = [kapinit k1init] ;
    params  = [A,alpha,delta,theta,Delta,g,l1,l2,s1,s2] ;
    options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                   'TolFun',1e-15) ;
    [kapk,~] = fsolve(@kappa_ss,kapkinit,options,params) ;
    kappa = kapk(1) ;
    k1tilbar = kapk(2) ;
    
    %---------------------------------------------------------------------%
    % If two methods give the same answer then record them
    %---------------------------------------------------------------------%
    if abs(kappa - kap) < 10^(-5)
        omkapmat1(sim,1) = theta*k1tilbar/...
                          (theta*k1tilbar + (1 - theta)*kappa*k1tilbar) ;
        omkapmat1(sim,2) = kappa ;
    else
        display(kap) ;
        display(kappa) ;
        display(kap - kappa) ;
        error('two methods for solving kappa_bar gave different answers') ;
    end
end
    
display([smat smat(:,1)./smat(:,2) omkapmat1]) ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of l1 given l2 = 1
%-------------------------------------------------------------------------%
% lmat      = sims x 1 vector, l1 values to test for kappa_bar and
%             omega_bar
% omkapmat2 = sims x 2 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to l1 values
% sims2     = integer, number of combinations of values of l1 for which to
%             compute kappa_bar and omega_bar
%-------------------------------------------------------------------------%
lmat = [2 5 10]' ;
s1 = 0.105 ;
s2 = 0.105 ;
omkapmat2 = zeros(size(lmat)) ;
sims2 = size(lmat) ;
sims2 = sims2(1) ;

for sim = 1:sims2
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar by simulating a time series
    %---------------------------------------------------------------------%
    % k1init   = scalar, initial value of type 1 household wealth
    % k2init   = scalar, initial value of type 2 household wealth
    % kapinit  = scalar, initial value of kappa
    % k1       = scalar, current period's value of type 1 capital
    % k2       = scalar, current period's value of type 2 capital
    % kap      = scalar, current period's value of kappa_t
    % per      = integer, iteration number of simulation
    % kdist    = scalar, distance measure of kappa' from kappa
    % mindist  = scalar, tolerance level of distance governing convergence
    % maxper   = integer, maximum number of periods to simulate
    % K        = scalar, current period value of aggregate capital K_t
    % Y        = scalar, current period value of aggregate output Y_t
    % w        = scalar, current period value of real wage w_t
    % r        = scalar, current period value of real rental rate r_t
    % r1       = scalar, current period value of real return for type 1 hh
    % r2       = scalar, current period value of real return for type 2 hh
    % k1pr     = scalar, next period's value of type 1 capital
    % k2pr     = scalar, next period's value of type 2 capital
    % kappr    = scalar, next period's value of kappa'
    %---------------------------------------------------------------------%
    l1 = lmat(sim) ;
    k1init = 100 ;
    k2init = 5 ;
    kapinit = k2init/k1init ;
    k1 = k1init ;
    k2 = k2init ;
    kap = kapinit ;
    per = 1 ;
    kdist = 10 ;
    mindist = 10^(-10) ;
    maxper = 50 ;

    while kdist > mindist && per < maxper
        per = per + 1 ;
        K = theta*k1 + (1 - theta)*k2 ;
        Y = A*(K^alpha)*((exp(g*per)*L)^(1 - alpha)) ;
        w = (1 - alpha)*(Y/L) ;
        r = alpha*(Y/K) - delta ;
        r1 = r + (Delta*(1-theta)*k2)/K ;
        r2 = r - (Delta*theta*k1)/K ;
        k1pr = s1*(w*l1 + (1 + r1)*k1) ;
        k2pr = s2*(w*l2 + (1 + r2)*k2) ;
        kappr = k2pr/k1pr ;
        kdist = abs(kappr - kap) ;
        k1 = k1pr ;
        k2 = k2pr ;
        kap = kappr ;
    end
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar using nonlinear system of
    % stationary steady-state equations
    %---------------------------------------------------------------------%
    % kapkinit = 1 x 2 vector, initial values for kappa_bar and k1tilbar
    % params   = 1 x 10 vector, parameters to be passed in to fsolve
    % options  = structural array, options for fsolve
    % kapk     = 1 x 2 vector, solution for kappa_bar and k1tilbar
    % kappa    = scalar, solution for kappa_bar
    % k1tilbar = scalar, solution for k1tilbar
    %---------------------------------------------------------------------%
    kapkinit = [kapinit k1init] ;
    params  = [A,alpha,delta,theta,Delta,g,l1,l2,s1,s2] ;
    options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                   'TolFun',1e-15) ;
    [kapk,~] = fsolve(@kappa_ss,kapkinit,options,params) ;
    kappa = kapk(1) ;
    k1tilbar = kapk(2) ;
    
    %---------------------------------------------------------------------%
    % If two methods give the same answer then record them
    %---------------------------------------------------------------------%
    if abs(kappa - kap) < 10^(-5)
        omkapmat2(sim,1) = theta*k1tilbar/...
                          (theta*k1tilbar + (1 - theta)*kappa*k1tilbar) ;
        omkapmat2(sim,2) = kappa ;
    else
        display(kap) ;
        display(kappa) ;
        display(kap - kappa) ;
        error('two methods for solving kappa_bar gave different answers') ;
    end
end
    
display([lmat omkapmat2]) ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of Delta
%-------------------------------------------------------------------------%
% Dmat      = sims x 1 vector, Delta values to test for kappa_bar and
%             omega_bar
% omkapmat3 = sims x 2 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to Delta values
% sims3     = integer, number of combinations of values of Delta for which
%             to compute kappa_bar and omega_bar
%-------------------------------------------------------------------------%
Dmat = [0.01 0.04 0.10 0.15 0.20]' ;
l1 = 1 ;
omkapmat3 = zeros(size(Dmat)) ;
sims3 = size(Dmat) ;
sims3 = sims3(1) ;

for sim = 1:sims3
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar by simulating a time series
    %---------------------------------------------------------------------%
    % k1init   = scalar, initial value of type 1 household wealth
    % k2init   = scalar, initial value of type 2 household wealth
    % kapinit  = scalar, initial value of kappa
    % k1       = scalar, current period's value of type 1 capital
    % k2       = scalar, current period's value of type 2 capital
    % kap      = scalar, current period's value of kappa_t
    % per      = integer, iteration number of simulation
    % kdist    = scalar, distance measure of kappa' from kappa
    % mindist  = scalar, tolerance level of distance governing convergence
    % maxper   = integer, maximum number of periods to simulate
    % K        = scalar, current period value of aggregate capital K_t
    % Y        = scalar, current period value of aggregate output Y_t
    % w        = scalar, current period value of real wage w_t
    % r        = scalar, current period value of real rental rate r_t
    % r1       = scalar, current period value of real return for type 1 hh
    % r2       = scalar, current period value of real return for type 2 hh
    % k1pr     = scalar, next period's value of type 1 capital
    % k2pr     = scalar, next period's value of type 2 capital
    % kappr    = scalar, next period's value of kappa'
    %---------------------------------------------------------------------%
    Delta = Dmat(sim) ;
    k1init = 100 ;
    k2init = 5 ;
    kapinit = k2init/k1init ;
    k1 = k1init ;
    k2 = k2init ;
    kap = kapinit ;
    per = 1 ;
    kdist = 10 ;
    mindist = 10^(-10) ;
    maxper = 50 ;

    while kdist > mindist && per < maxper
        per = per + 1 ;
        K = theta*k1 + (1 - theta)*k2 ;
        Y = A*(K^alpha)*((exp(g*per)*L)^(1 - alpha)) ;
        w = (1 - alpha)*(Y/L) ;
        r = alpha*(Y/K) - delta ;
        r1 = r + (Delta*(1-theta)*k2)/K ;
        r2 = r - (Delta*theta*k1)/K ;
        k1pr = s1*(w*l1 + (1 + r1)*k1) ;
        k2pr = s2*(w*l2 + (1 + r2)*k2) ;
        kappr = k2pr/k1pr ;
        kdist = abs(kappr - kap) ;
        k1 = k1pr ;
        k2 = k2pr ;
        kap = kappr ;
    end
    
    %---------------------------------------------------------------------%
    % Solve for kappa_bar and omega_bar using nonlinear system of
    % stationary steady-state equations
    %---------------------------------------------------------------------%
    % kapkinit = 1 x 2 vector, initial values for kappa_bar and k1tilbar
    % params   = 1 x 10 vector, parameters to be passed in to fsolve
    % options  = structural array, options for fsolve
    % kapk     = 1 x 2 vector, solution for kappa_bar and k1tilbar
    % kappa    = scalar, solution for kappa_bar
    % k1tilbar = scalar, solution for k1tilbar
    %---------------------------------------------------------------------%
    kapkinit = [kapinit k1init] ;
    params  = [A,alpha,delta,theta,Delta,g,l1,l2,s1,s2] ;
    options = optimset('Display','off','MaxFunEvals',100000,'MaxIter',1000,...
                   'TolFun',1e-15) ;
    [kapk,~] = fsolve(@kappa_ss,kapkinit,options,params) ;
    kappa = kapk(1) ;
    k1tilbar = kapk(2) ;
    
    %---------------------------------------------------------------------%
    % If two methods give the same answer then record them
    %---------------------------------------------------------------------%
    if abs(kappa - kap) < 10^(-5)
        omkapmat3(sim,1) = theta*k1tilbar/...
                          (theta*k1tilbar + (1 - theta)*kappa*k1tilbar) ;
        omkapmat3(sim,2) = kappa ;
    else
        display(kap) ;
        display(kappa) ;
        display(kap - kappa) ;
        error('two methods for solving kappa_bar gave different answers') ;
    end
end
    
display([Dmat omkapmat3]) ;

%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save BalGrSims.mat ;
