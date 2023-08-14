%% Heston: dX=k*(theta-X)*dt+sigma*sqrt(X)*dW
% dS = r*S*dt + sqrt(X)*(rho*dW+sqrt(1-rho^2)*dB)
clear all;
close all;
X0 = .2;
S0 = 100;
K =90:120;
k = 3;
theta = .2;
sigma = sqrt(2*theta*k);
sigma_mul = [0.35,0.75,1];
T = 3/12;
r = 0.05;
rho = [-0.2,0,0.2];
dt = 1/365;
n = ceil(T/dt);
% n=ceil(252/4);
M = 20000;
implied_volatility = zeros(length(rho)*length(sigma_mul),length(K));
implied_volatility_gt = zeros(length(rho)*length(sigma_mul),length(K));
log_moneyness = log(K*exp(-r*T)/S0);
row = 1;
for i = 1:length(rho)
    for j = 1:length(sigma_mul)
        
        implied_volatility(row,:) = heston_model(n,k,theta,M,S0,X0,K,r,dt,T,sigma*sigma_mul(i),rho(j));
        figure(1);
        subplot(length(sigma_mul),length(rho),row)
        plot(log_moneyness, implied_volatility(row,:))
        s = "Sigma Multiplier is (" + num2str(sigma_mul(i));
        s1 = "), Rho is (" + num2str(rho(j)) + ")";
        title(s+s1)

        row = row+1;

    end
end

row = 1;
for i = 1:length(rho)
    for j = 1:length(sigma_mul)
        
        implied_volatility_gt(row,:) = heston_model_gt(n,k,theta,M,S0,X0,K,r,dt,T,sigma*sigma_mul(i),rho(j));
        figure(2);
        subplot(length(sigma_mul),length(rho),row)
        plot(log_moneyness, implied_volatility_gt(row,:))
        s = "Sigma Multiplier is (" + num2str(sigma_mul(i));
        s1 = "), Rho is (" + num2str(rho(j)) + ")";
        title(s+s1)

        row = row+1;

    end
end


%% Heston

function implied_volatility = heston_model(n,k,theta,M,S0,X0,K,r,dt,T,sigma,rho)
X = zeros(M,n);
X(:,1) = X0;
Y = zeros(M,n);
Y(:,1) = log(S0);

for t= 2:n
    dW = randn(M,1)*sqrt(dt);
    dB = randn(M,1)*sqrt(dt);
    X(:,t) = (1-k*dt)*X(:,t-1) + k*theta*dt + sigma*sqrt(max(X(:,t-1),0)).*dB;
    Y(:,t) = Y(:,t-1) + (r - .5*X(:,t-1))*dt + sqrt(max(X(:,t-1),0)).*(rho*dB+sqrt(1-rho^2)*dW);
end

S = exp(Y(:,n));
[Sv,Kv] = meshgrid(S,K);
D =exp(-r*T); 
C_mc = D*max(Sv(:,:)-Kv(:,:),0);
C_mc_avg = mean(C_mc,2);
implied_volatility = blsimpv(S0,K',r,T,C_mc_avg);
end


function implied_volatility = heston_model_gt(n,k,theta,M,S0,X0,K,r,dt,T,sigma,rho);
Settle = datetime(2023,2,27);
Maturity = datemnth(Settle, 3);
c = optByHestonNI(r,S0,Settle, Maturity,'call',K,X0,theta,k,sigma,rho,'DividendYield',0);
implied_volatility = blsimpv(S0,K,r,T,c');
end

%% Effect of vol-of-vol parameter
% Higher vol-of-vol parameter leads to steeper slope and greater curvature
% in the smile, meaning that out-of-the money options have higher implied
% volatilities relative to at-the-money options.
%% Effect of rho 
% As rho increases, the implied volatality for out of money call option
% increases. Specifically, for a negative rho, the slope for the smile is
% negative and it keeps increasing as rho increases. For zero correlation,
% we get the closest to a smile for implied volatality.
%% Leverage effect
% A negative rho creates a volatality leverage effect as when the stock
% price decreases, the volatality increases and vice versa. This can be
% observed as when the strike price increases and the option goes more out
% of the money, the volatality decreases as seen by the downward slope.