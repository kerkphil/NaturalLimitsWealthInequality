%-------------------------------------------------------------------------%
% Perform unbalanced growth simulation experiments for simple Piketty model
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
% etavec   = etasize x 1 vector, values of eta to test, production function
%            parameter proportional to the elasticity of substitution
%            between capital and labor
% etasize  = integer, number of values of eta to check
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
etavec = [-1 0 2/3]' ;
etasize = size(etavec) ;
etasize = etasize(1) ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of l1 given l2 = 1
%-------------------------------------------------------------------------%
% lvec      = simsL x 1 vector, l1 values to test for kappa_bar and
%             omega_bar
% omkapmatL = simsL x 6 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to l1 values and varying eta (-1, 0, 2/3)
% simsL     = integer, number of combinations of values of l1 for which to
%             compute kappa_bar and omega_bar
% eta       = scalar, production function parameter proportional to the
%             elasticity of substitution between capital and labor
%-------------------------------------------------------------------------%
lvec = [2 5 10]' ;
s1 = 0.105 ;
s2 = 0.105 ;
Delta = 0 ;
simsL = size(lvec) ;
simsL = simsL(1) ;
omkapmatL = zeros(simsL,6) ;

k1init = 1 ;
k2init = 1 ;
kapinit = k2init/k1init ;
Kinit = theta*k1init + (1 - theta)*k2init ;

for et = 1:etasize
    eta = etavec(et) ;
    for sim = 1:simsL

        %-----------------------------------------------------------------%
        % Solve for kappa_bar and omega_bar by simulating a time series
        %-----------------------------------------------------------------%
        % k1init   = scalar, initial value of type 1 household wealth
        % k2init   = scalar, initial value of type 2 household wealth
        % kapinit  = scalar, initial value of kappa
        % k1       = scalar, current period's value of type 1 capital
        % k2       = scalar, current period's value of type 2 capital
        % kap      = scalar, current period's value of kappa_t
        % per      = integer, iteration number of simulation
        % kdist    = scalar, distance measure of kappa' from kappa
        % mindist  = scalar, tolerance level of distance governing
        %            convergence
        % maxper   = integer, maximum number of periods to simulate
        % K        = scalar, current period value of aggregate capital K_t
        % Gamma    = scalar, current period value of Gamma_t
        % Y        = scalar, current period value of aggregate output Y_t
        % w        = scalar, current period value of real wage w_t
        % r        = scalar, current period value of real rental rate r_t
        % r1       = scalar, current period value of real return for type 1
        %            households
        % r2       = scalar, current period value of real return for type 2
        %            households
        % k1pr     = scalar, next period's value of type 1 capital
        % k2pr     = scalar, next period's value of type 2 capital
        % kappr    = scalar, next period's value of kappa'
        %-----------------------------------------------------------------%
        l1 = lvec(sim) ;
        k1 = k1init ;
        k2 = k2init ;
        kap = kapinit ;
        per = 1 ;
        kdist = 10 ;
        mindist = 10^(-10) ;
        maxper = 200 ;

        while kdist > mindist && per < maxper
            per = per + 1 ;
            K = theta*k1 + (1 - theta)*k2 ;
            Gamma = (alpha*(K^eta))/(alpha*(K^eta) + (1 - alpha)*(L^eta)) ;
            if eta == 0
                Y = exp(g*per)*A*(K^alpha)*(L^(1 - alpha)) ;
            else
                Y = exp(g*per)*A*(alpha*(K^eta) + ...
                    (1 - alpha)*(L^eta))^(1/eta) ;
            end
            w = (1 - Gamma)*(Y/L) ;
            r = Gamma*(Y/K) - delta ;
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
        
        omkapmatL(sim,(2*et)-1) = theta*k1pr/(theta*k1pr + ...
                                  (1 - theta)*k2pr) ;
        omkapmatL(sim,2*et) = kap ;
    end
end
    
display([lvec omkapmatL]) ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of Delta = r1 - r2
%-------------------------------------------------------------------------%
% Dvec      = simsD x 1 vector, Delta values to test for kappa_bar and
%             omega_bar
% simsD     = integer, number of combinations of values of Delta for which
%             to compute kappa_bar and omega_bar
% omkapmatD = simsD x 6 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to Delta values and varying eta (-1, 0, 2/3)
% eta       = scalar, production function parameter proportional to the
%             elasticity of substitution between capital and labor
%-------------------------------------------------------------------------%
Dvec = [0.01 0.04 0.10 0.15 0.20]' ;
l1 = 1 ;
s1 = 0.105 ;
s2 = 0.105 ;
simsD = size(Dvec) ;
simsD = simsD(1) ;
omkapmatD = zeros(simsD,6) ;

for et = 1:etasize
    eta = etavec(et) ;
    for sim = 1:simsD

        %-----------------------------------------------------------------%
        % Solve for kappa_bar and omega_bar by simulating a time series
        %-----------------------------------------------------------------%
        Delta = Dvec(sim) ;
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
            Gamma = (alpha*(K^eta))/(alpha*(K^eta) + (1 - alpha)*(L^eta)) ;
            if eta == 0
                Y = exp(g*per)*A*(K^alpha)*(L^(1 - alpha)) ;
            else
                Y = exp(g*per)*A*(alpha*(K^eta) + ...
                    (1 - alpha)*(L^eta))^(1/eta) ;
            end
            w = (1 - Gamma)*(Y/L) ;
            r = Gamma*(Y/K) - delta ;
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
        
        omkapmatD(sim,(2*et)-1) = theta*k1pr/(theta*k1pr + ...
                                  (1 - theta)*k2pr) ;
        omkapmatD(sim,2*et) = kap ;
    end
end
    
display([Dvec omkapmatD]) ;

%-------------------------------------------------------------------------%
% Do long-run solution of the model for different values of s1 and s2
%-------------------------------------------------------------------------%
% Smat      = simsS x 2 matrix, (s1,s2) pairs to test for kappa_bar and
%             omega_bar
% simsS     = integer, number of combinations of (s1,s2) pairs for which to
%             compute kappa_bar and omega_bar
% omkapmatS = simsS x 6 matrix, (omega_bar,kappa_bar) combinations
%             corresponding to (s1,s2) pairs and varying eta (-1, 0, 2/3)
% eta       = scalar, production function parameter proportional to the
%             elasticity of substitution between capital and labor
%-------------------------------------------------------------------------%
smat = [0.99 0.01; 0.75 0.1; 0.75 0.2; 0.5 0.3; 0.5 0.49; 0.236 0.072] ;
l1 = 1 ;
Delta = 0  ;
simsS = size(smat) ;
simsS = simsS(1) ;
maxper = 200 ;
omkapmatS = zeros(simsS,6) ;
kapmat = zeros(maxper,simsS) ;

for et = 1:etasize
    eta = etavec(et) ;
    for sim = 1:simsS

        %-----------------------------------------------------------------%
        % Solve for kappa_bar and omega_bar by simulating a time series
        %-----------------------------------------------------------------%
        s1 = smat(sim,1) ;
        s2 = smat(sim,2) ;
        k1 = k1init ;
        k2 = k2init ;
        kap = kapinit ;
        per = 1 ;
        kdist = 10 ;
        mindist = 10^(-10) ;
        omvec = zeros(maxper,1) ;
        k1vec = zeros(maxper,1) ;
        k2vec = zeros(maxper,1) ;
        k1vecg = zeros(maxper,1) ;
        k2vecg = zeros(maxper,1) ;
        kapmat(per,sim) = kap ;
        omvec(per) = (theta*k1)/(theta*k1 + (1 - theta)*k2) ;
        k1vec(per) = k1 ;
        k2vec(per) = k2 ;

        while kdist > mindist && per < maxper && kap > 0
            per = per + 1 ;
            K = theta*k1 + (1 - theta)*k2 ;
            Gamma = (alpha*(K^eta))/(alpha*(K^eta) + (1 - alpha)*(L^eta)) ;
            if eta == 0
                Y = exp(g*per)*A*(K^alpha)*(L^(1 - alpha)) ;
            else
                Y = exp(g*per)*A*(alpha*(K^eta) + ...
                    (1 - alpha)*(L^eta))^(1/eta) ;
            end
            w = (1 - Gamma)*(Y/L) ;
            r = Gamma*(Y/K) - delta ;
            r1 = r + (Delta*(1-theta)*k2)/K ;
            r2 = r - (Delta*theta*k1)/K ;
            k1pr = s1*(w*l1 + (1 + r1)*k1) ;
            k2pr = s2*(w*l2 + (1 + r2)*k2) ;
            kappr = k2pr/k1pr ;
            kdist = abs(kappr - kap) ;
            k1vecg(per-1) = (k1pr - k1)/k1 ;
            k2vecg(per-1) = (k2pr - k2)/k2 ;
            k1 = k1pr ;
            k2 = k2pr ;
            kap = kappr ;
            kapmat(per,sim) = kap ;
            omvec(per) = (theta*k1)/(theta*k1 + (1 - theta)*k2) ;
        end
        
        omkapmatS(sim,(2*et)-1) = theta*k1pr/(theta*k1pr + ...
                                  (1 - theta)*k2pr) ;
        omkapmatS(sim,2*et) = kap ;
    end
end
    
display([smat smat(:,1)./smat(:,2) omkapmatS]) ;

figure(1)
plot(1:1:12,kapmat(1:12,1),1:1:16,kapmat(1:16,2),1:1:21,kapmat(1:21,3),1:1:44,kapmat(1:44,4),1:1:44,kapmat(1:44,5),1:1:25,kapmat(1:25,6))
legend('s_1/s_2=99.00','s_1/s_2=7.50','s_1/s_2=3.75','s_1/s_2=1.67','s_1/s_2=1.02','s_1/s_2=3.28')
xlabel('Period'); ylabel('\kappa_t')



%-------------------------------------------------------------------------%
runtime = etime(clock,starttime) ; % end timer
save UnbalGrSims.mat ;
