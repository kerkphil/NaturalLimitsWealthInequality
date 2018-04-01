function diff = kappaunbal_ss(kapk,params)
%--------------------------------------------------------------------------
% kapk     = 1 x 2 vector, values for kappa and k1bartil
% params   = 1 x 11 vector of parameters passed in to function
% kap      = scalar, value of kappa_bar
% k1bartil = scalar, value of k1bartil
% A        = scalar, production function scale parameter
% alpha    = scalar, capital share of income, production function parameter
% delta    = scalar, generational rate of capital depreciation
% eta      = scalar, production function parameter proportional to the
%            elasticity of substitution between capital and labor
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
k1bartil = kapk(2) ;

A = params(1) ;
alpha = params(2) ;
delta = params(3) ;
eta = params(4) ;
theta = params(5) ;
Delta = params(6) ;
g = params(7) ;
l1 = params(8) ;
l2 = params(9) ;
s1 = params(10) ;
s2 = params(11) ;

%--------------------------------------------------------------------------
% Calculate the zero (difference) function
%--------------------------------------------------------------------------
% L       = scalar, aggregate labor
% Kbartil = scalar, stationary steady-state aggregate capital stock
% Gamma   = scalar, long-run relative product of capital
% Ybartil = scalar, stationary steady-state aggregate output
% wbartil = scalar, stationary steady-state real wage
% rbar    = scalar, steady-state real rental rate
% r1bar   = scalar, steady-state real return for type 1 households
% r2bar   = scalar, steady-state real return for type 2 households
% diff1   = scalar, zero version of equation 1 for kappa_bar
% diff2   = scalar, zero version of equation 2 for k1bartil
% diff    = 2 x 1 vector, two values of zero version of equations 1 and 2
%--------------------------------------------------------------------------
L = theta*l1 + (1 - theta)*l2 ;
Kbartil = theta*k1bartil + (1 - theta)*k1bartil*kap ;
if eta < 0
    Gamma = 0 ;
    Ybartil = A*((1-alpha)^(1/eta))*L ;
    wbartil = (1-Gamma)*Ybartil/L ;
    rbar = Gamma*(Ybartil/Kbartil) - delta ;
    r1bar = rbar + (Delta*(1 - theta)*k1bartil*kap)/Kbartil ;
    r2bar = rbar - (Delta*theta*k1bartil)/Kbartil ;
elseif eta == 0
    Gamma = alpha ;
    Ybartil = A*(Kbartil^alpha)*(L^(1 - alpha)) ;
elseif eta > 0
    Gamma = 1 ;
end

diff1 = kap - (s2*wbartil*l2)/(s1*(wbartil*l1 + (1 + r1bar)*k1bartil) - ...
        s2*(1 + r2bar)*k1bartil) ;
diff2 = k1bartil - (s1*wbartil*l1)/(exp(g) - s1*(1 + r1bar)) ;
diff = [diff1 diff2] ;