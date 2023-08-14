%% Asian Option approxiamtion
% Simulate $S_T$,
%
% $$S_T = S_0 \exp\left(\left(\mu-\frac{\sigma^2}{2}\right) T + \sigma W_T\right)\ .$$
%
% with $W_T\sim N(0,T)$. 


%% Initializing parameters
clear all
close all
T = 1; 
n = 5000; 
S0 = 100; 
r = .05; 
sigma = .2; 
K = 90:120; 
m = 4;
C_mc = MonteCarloAsianCallOption(T,m,S0,n,K,r,sigma);
C_lna = logNormalApproximation(S0,r,T,m,sigma,K);
C_bach = bachelier(S0, r, T, m, K, sigma,n);
T1 = table(C_mc,C_lna',C_bach','VariableNames',{'Monte Carlo','Log Normal','Bachelier'});
disp(T1);
%% Plots for T = 1 and m=4
fig = figure(1);
hold on;
title('Asian Call Option Price vs Strike Price (T=1, m=4)');
xlabel('Strike Price');
ylabel('Call option price');
plot(K,C_mc);
plot(K,C_lna);
plot(K, C_bach);
legend('Monte Carlo','Log Normal Approx.','Bachelier');
hold off;
%% All simulations for T=5 and m=20
T=5;
m=20;

C_mc = MonteCarloAsianCallOption(T,m,S0,n,K,r,sigma);
C_lna = logNormalApproximation(S0,r,T,m,sigma,K);
C_bach = bachelier(S0, r, T, m, K, sigma,n);
T1 = table(C_mc,C_lna',C_bach','VariableNames',{'Monte Carlo','Log Normal','Bachelier'});
disp(T1);
%% Plots for T = 5 and m=20
fig = figure(2);
hold on;
title('Asian Call Option Price vs Strike Price (T=5, m=20)');
xlabel('Strike Price');
ylabel('Call option price');
plot(K,C_mc);
plot(K,C_lna);
plot(K, C_bach);
legend('Monte Carlo','Log Normal Approx.','Bachelier');
hold off;



%% Monte Carlo Simulation function

function C_mc = MonteCarloAsianCallOption(T,m,S0,n,K,r,sigma)
% Generate Samples

ti = linspace(T/m, T, m);
dt = ti(2) - ti(1);
W = cumsum(randn(m,n),1)*sqrt(dt);
S = S0*exp((r-0.5*sigma^2)*repmat(ti,n,1)'+sigma*W);
S_bar = mean(S,1);

% Compute the payoffs and 
[Sv,Kv] = meshgrid(S_bar,K);
D =exp(-r*T); 
C_mc = D*mean(max(Sv(:,:)-Kv(:,:),0),2);


end


%% Log normal approximation

function C_lna = logNormalApproximation(S0,r,T,m, sigma,K)

ti = linspace(T/m, T, m);
Fi = S0*exp(r*ti);
Gij = zeros(m*m,1);

k = 1;
for i = 1:m
    for j = 1:m
        Gij(k) = Fi(i)*Fi(j)*exp((sigma^2)*min(ti(i),ti(j)));
        k=k+1;
    end
end

M1 = sum(Fi)/m;
M2 = sum(Gij)/(m^2);

S0_hat = exp(-1*r*T)*M1;
sigma_sq_hat = log(M2/(M1^2))/T;

% Black Scholes calculation

C_lna = blsprice(S0_hat,K,r,T,sqrt(sigma_sq_hat));

end


%% Bachelier


function C_bach = bachelier(S0, r, T, m, K, sigma,n)

ti = linspace(T/m, T, m);
dt = ti(2) - ti(1);
W = cumsum(randn(m,n),1)*sqrt(dt);
S = S0*exp((r-0.5*sigma^2)*repmat(ti,n,1)'+sigma*W);
S_bar = mean(S,1);

[mean_hat, variance] = normfit(S_bar);

F = mean_hat;

Z = (F-K)/variance;

C_bach = exp(-r*T)*((F-K).*normcdf(Z) + variance.*normpdf(Z));

end




