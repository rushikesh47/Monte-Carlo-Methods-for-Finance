%%
% $$\textbf{Problem\hspace{.05cm} 1:}$$
%
% $$ \frac{dS_t}{S_t} = (r-q)dt + \sigma dW_t $$
%
% 
%
% $$ Asian\hspace{.1cm} option \hspace{.1cm} payoff = \left( \frac{1}{m} \sum_{i=1}^{m} S_{t_i} - K  \right)^{+} $$
%
%
%
% $$ Auxillary\hspace{.1cm}payoff\hspace{.1cm}using\hspace{.1cm}geometric\hspace{.1cm}mean = \left(  \left( \prod_{i=1}^{m} S_{t_i} \right)^{\frac{1}{m}} - K \right)^{+} $$
%
% $$ Antithetic: \tilde{W} = - W $$
%
%%





%% Declaring variables


clear all
S0 = 100;
sigma = .2;
r = .05;
q = .02;
T = 1;
m = 4;
K = 90:5:120;
D = exp(-r*T);
t = linspace(T/m,T,m);

%% MonteCarlo

N = 100;  %% Monte Carlo sample size
L = 100;  %% number of trials

C_mc = zeros(length(K),L);
C_cv = zeros(length(K),L);
C_ant  = zeros(length(K),L);

for ctr = 1:L
    %%%% Monte Carlo Asian Call 
    W = cumsum(randn(N,m)*sqrt(T/m),2);
    S = S0*exp((r-q-.5*sigma^2)*t + sigma*W);
    M = mean(S,2);
    %M = exp(mean(log(S),2));
    [Mv,Kv] = meshgrid(M,K);
    psi_asian = D*max(Mv-Kv,0);
    
    C_mc(:,ctr) = mean(psi_asian,2);
end

%% Control Variate

for ctr = 1:L
    W = cumsum(randn(N,m)*sqrt(T/m),2);   
    S = S0*exp((r-q-.5*sigma^2)*t + sigma*W);
    M = mean(S,2);
    M_geo = geomean(S,2);
    %M = exp(mean(log(S),2));
    [Mv,Kv] = meshgrid(M,K);
    [Mv_geo,Kv_geo] = meshgrid(M_geo,K);
    psi_asian = D*max(Mv-Kv,0);
    
    for i = 1:m
        M1(i) = exp((((r-q)-.5*sigma^2)*((T*i)/m^2))+((sigma^2*T*i*i)/(2*m^3)));
        M2(i) = exp((2*((r-q)-.5*sigma^2)*((T*i)/(m^2)))+(((2*sigma^2*T)/m^3)*(i*i)));
    end

    M1_new = prod(M1);
    F=S0*M1_new;
    M2_new = prod(M2);
    sig_new = sqrt(log(M2_new/(M1_new^2))/T);
    
    psi_geo = D*max(Mv_geo-Kv_geo,0);
    Epsi_geo = blkprice(F,K,r-q,T,sig_new);
    b = (psi_asian - mean(psi_asian,2))*(psi_geo - Epsi_geo')'/N;
    b = diag(b)./diag((psi_geo-Epsi_geo')*(psi_geo-Epsi_geo')'/N);
    psi_cv = psi_asian - diag(b)*(psi_geo - Epsi_geo');

    C_cv(:,ctr) = mean(psi_cv,2);
  
end

%% Antithetic

for ctr = 1:L
    %%%% Monte Carlo Asian Call 
    W = cumsum(randn(N,m)*sqrt(T/m),2);
    S = S0*exp((r-q-.5*sigma^2)*t + sigma*W);
    M_ant_1 = mean(S,2);
    [Mv_ant_1,Kv_ant_1] = meshgrid(M_ant_1,K);
    psi_asian_ant_1 = .5*D*max(Mv_ant_1-Kv_ant_1,0);
    
    S1 = S0*exp((r-q-.5*sigma^2)*t - sigma*W);
    M_ant_2 = mean(S1,2);
    [Mv_ant_2,Kv_ant_2] = meshgrid(M_ant_2,K);
    psi_asian_ant_2 = .5*D*max(Mv_ant_2-Kv_ant_2,0);


    C_ant(:,ctr) = mean((psi_asian_ant_1+psi_asian_ant_2),2);

end

err = [std(C_mc,[],2),std(C_cv,[],2),std(C_ant,[],2)];
avg = [mean(C_mc,2),mean(C_cv,2), mean(C_ant,2)];

hold on;
plot(K,avg);
legend('Monte Carlo','Control Variate', 'Antithetic');
hold off;


table = table(K', err(:,1), err(:,2), err(:,3),'VariableNames',{'strike', 'MC', 'Ctrl Var', 'Antithetic' });
disp(table);


%% Reduction in Variance 
% There is a significant correlation between our arithmetic and geometric 
% means, which results in a considerable reduction in variance. 
% Specifically, the reduction is equal to (1-rho^2) in control variates. In
% contrast, using antithetic variables reduces variance by 1/N, 
% resulting in a variance reduction of 0.01. While using the arithmetic and
% geometric means leads to a larger variance reduction due to high 
% correlation, the random variables used in both sections are different, so
% we may not observe the exact theoretical variance reductions.

