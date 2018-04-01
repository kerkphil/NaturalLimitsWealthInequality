function diff = kappa_ss(kapk,params)
%--------------------------------------------------------------------------
% kap      = scalar, value of kappa_bar
% k1tilbar = scalar, value of k1tilbar
% A        = scalar, production function scale parameter
% alpha    = scalar, capital share of income, production function parameter
% delta    = scalar, generational rate of capital depreciation
% theta    = scalar, percent of population that is type 1
% Delta    = scalar, interest rate wedge: Delta = r1t - r2t
% g        = scalar, generational growth rate of labor augmenting
%            technological change
% l1       = scalar, labor endowment of type 1 households
% l2       = scalar, labor endowment of type 2 households
% s1       = scalar, savings rate of type 1 individuals
% s2       = scalar, savings rate of type 2 individuals
%--------------------------------------------------------------------------
kap = kapk(1) ;
k1tilbar = kapk(2) ;

A = params(1) ;
alpha = params(2) ;
delta = params(3) ;
theta = params(4) ;
Delta = params(5) ;
g = params(6) ;
l1 = params(7) ;
l2 = params(8) ;
s1 = params(9) ;
s2 = params(10) ;

%--------------------------------------------------------------------------
% Calculate the zero (difference) function
%--------------------------------------------------------------------------
% L       = scalar, aggregate labor
% Ktilbar = scalar, stationary steady-state aggregate capital stock
% Ytilbar = scalar, stationary steady-state aggregate output
% wtilbar = scalar, stationary steady-state real wage
% rbar    = scalar, steady-state real rental rate
% r1bar   = scalar, steady-state real return for type 1 households
% r2bar   = scalar, steady-state real return for type 2 households
% diff1   = scalar, zero version of equation 1 for kappa_bar
% diff2   = scalar, zero version of equation 2 for k1tilbar
% diff    = 2 x 1 vector, two values of zero version of equations 1 and 2
%--------------------------------------------------------------------------
L = theta*l1 + (1 - theta)*l2 ;
Ktilbar = theta*k1tilbar + (1 - theta)*k1tilbar*kap ;
Ytilbar = A*(Ktilbar^alpha)*(L^(1 - alpha)) ;
wtilbar = (1 - alpha)*(Ytilbar/L) ;
rbar = alpha*(Ytilbar/Ktilbar) - delta ;
r1bar = rbar + (Delta*(1 - theta)*k1tilbar*kap)/Ktilbar ;
r2bar = rbar - (Delta*theta*k1tilbar)/Ktilbar ;

diff1 = kap - (s2*wtilbar*l2)/(s1*(wtilbar*l1 + (1 + r1bar)*k1tilbar) - ...
        s2*(1 + r2bar)*k1tilbar) ;
diff2 = k1tilbar - (s1*wtilbar*l1)/(exp(g) - s1*(1 + r1bar)) ;
diff = [diff1 diff2] ;