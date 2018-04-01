clear

% set simulation parameters
nyr = 25;        %number of period per year
nobs = 100;      %number of periods per simulation
ndyn = 1000;     %number of dynasties per simulation
nmc = 1000;      %number of monte carlos
startper = .01;  %initial population percent of type 1 dynasties.
startscale = 1;  %type 1 wealth relative to type 2 wealth 
startW2 = 0.01;  %initial wealth for type 2
startW1 = startW2*startscale;    %initial wealth for type 1

% set model parameters
E1 = 1;         %effective labor for type 1 agent
E2 = 1;         %effective labor for type 2 agent
alf = .35;      %capital quasi-share in GDP
eta = .1;       %maps to elasticity of substitution
A = 1;          %scaling factor in production function
del_ann = .08;  %annual depreciation rate
Del_ann = .0;   %annual wedge between type returns
g_ann = .01;    %annual growth rate of technology
rho1 = .9;      %type 1 Markov transition probability
%rho2 = 1;      %type 2 Markov transition probability
gam1 = 3;       %utility weight on type 1 bequests
gam2 = 1/3;     %utility weight on type 2 bequests

% set transition prob for type 2 to maintain a fixed percent
SSper = startper;
rho2 = (1-SSper*(2-rho1))/(1-SSper);

% calclate per period parameter values
del = 1-(1-del_ann)^nyr;
Del = (1+Del_ann)^nyr-1;
g = (1+g_ann)^nyr-1;

% create Markov transition matrix
Mark = [rho1 1-rho1; 1-rho2 rho2];

% set up arrys to store monte carlo moving averages
omega_avg = zeros(nobs,1);
Tperc_avg = zeros(nobs,1);
gini_avg = zeros(nobs,1);
Wmat_avg = zeros(nobs,ndyn);

% being monte carlos
for n=1:nmc
    if mod(n,100)==0;
        n
    end
    % setup within simulation matrices
    Wmat = zeros(nobs,ndyn);   %history of dynasty wealths
    Tmat = zeros(nobs,ndyn);   %history of dynasty types
    omega = zeros(nobs,1);     %history of relative wealth
    Tperc = zeros(nobs,1);     %history of number of type 1 dynasties
    gini = zeros(nobs,1);      %history of Gini coefficients

    % set initial types
    cutoff = round(startper*ndyn);
    Tmat(1,1:cutoff) = ones(1,cutoff);

    % begin simulation
    % find wealths by type
    Wmat(1,:) = startW1*Tmat(1,:) + startW2*(ones(1,ndyn)-Tmat(1,:));
    Wtot = sum(Wmat(1,:));
    W1tot = Wmat(1,:)*Tmat(1,:)';
    W2tot = Wtot - W1tot;
    Tperc(1) = sum(Tmat(1,:))/ndyn;
    omega(1) = W1tot/Wtot;
    gini(1) = ginicalc(Wmat(1,:));

    for t = 1:nobs-1
        % computer type indicators
        ind1 = Tmat(t,:);
        ind2 = ones(1,ndyn)-Tmat(t,:);
        % compute GE wages and returns
        K = Wtot;
        L = sum(ind1)*E1 + sum(ind2)*E2;
        if eta == 0
            Y = exp(g*(t-1))*A*K^alf*L^(1-alf);
        else
            Y = exp(g*(t-1))*A*(alf*K^eta+(1-alf)*L^eta)^(1/eta);
        end
        r = alf*Y*K^(eta-1)/(alf*K^eta+(1-alf)*L^eta);
        w = (1-alf)*Y*L^(eta-1)/(alf*K^eta+(1-alf)*L^eta);
        r1 = r - del + Del*W2tot/K;
        r2 = r - del - Del*W1tot/K;
        %[t r1 r2 r w Y L K]
        % get disposable wealth for all dynasties
        disp = w*(ind1*E1 + ind2*E2) + ...
            Wmat(t,:).*(ones(1,ndyn) + ind1*r1 + ind2*r2);
        % mean disposable wealth
        dbar = mean(disp);
        % utility weight on bequests for all dynasties
        gam = ind1*gam1 + ind2*gam2;
        % new wealth, coming entirely from bequests
        Wmat(t+1,:) = (gam./(1+gam)).*disp;
        % draw a vector of random draws and calculate dynasty transitions
        dynrand = rand(1,ndyn);
        indic1 = logical(dynrand<rho1);  %1 stays 1
        indic2 = logical(dynrand>rho2);  %0 stays 1
        Tmat(t+1,:) = Tmat(t,:).*indic1 ...
                    + (ones(1,ndyn)-Tmat(t,:)).*(indic2);
        % find wealths by type next period
        Wtot = sum(Wmat(t+1,:));
        W1tot = Wmat(t+1,:)*Tmat(t+1,:)';
        W2tot = Wtot - W1tot;
        Tperc(t+1) = sum(Tmat(t+1,:))/ndyn;
        omega(t+1) = W1tot/Wtot;
        gini(t+1) = ginicalc(Wmat(t+1,:));
    end
    Wmatperc = Wmat./repmat(sum(Wmat,2),[1,ndyn]);
    omega_avg = ((n-1)/n)*omega_avg + (1/n)*omega;
    Tperc_avg = ((n-1)/n)*Tperc_avg + (1/n)*Tperc;
    gini_avg = ((n-1)/n)*gini_avg + (1/n)*gini;
    Wmat_avg = ((n-1)/n)*Wmat_avg + (1/n)*Wmat./Wmatperc;
end
twoplots(omega_avg,gini_avg)
%threeplots(omega_avg,Tperc_avg,gini_avg)