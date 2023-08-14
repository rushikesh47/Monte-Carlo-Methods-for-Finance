%%
% 
% $$\textbf{Problem\hspace{.05cm} 2:}$$
%
% $$\textbf{Heston\hspace{.1cm}stochastic\hspace{.1cm}volatility\hspace{.1cm}model:}$$
% 
% $$ \frac{dS_t}{S_t} = rdt + \sqrt{X_t} \left( \rho dB_t + \sqrt{1- \rho^2} dW_t   \right) $$
% 
% $$ dX_t = \kappa(\theta - X_t)dt + \sigma\sqrt{X_t}dB_t $$ 
%
% $$\textbf{Stochastic\hspace{.1cm}volatility\hspace{.1cm}model\hspace{.1cm}call\hspace{.1cm}price\hspace{.1cm}using\hspace{.1cm}Romano\hspace{.1cm}and\hspace{.1cm}Touzi\hspace{.1cm}formula:}$$
% 
% $$ C^{stoch\hspace{0.1cm} - \hspace{0.1cm} vol}(s,x,T) = \mathrm{E} \left[ C^{bs} (se^Z, \sigma(X), T) | S_0 = s, X_0 = x \right] $$
%
% $$C^{bs}(s,x,T)\hspace{.1cm} denotes\hspace{.1cm} the \hspace{.1cm} Black-Scholes\hspace{.1cm}  call \hspace{.1cm} price$$
%
% $$ Z = \rho \int_{0}^{T}\sqrt{X_t}dB_t\ - \frac{\rho^2}{2} \int_{0}^{T}X_tdt\ $$ 
%
% $$ \sigma^2(X) = \frac{1-\rho^2}{T} \int_{0}^{T}X_tdt\    $$
%
%
% $$\textbf{Using\hspace{.05cm} the \hspace{.05cm}Heston\hspace{.05cm}
% Model:}$$
%
% $$ \sigma\sqrt{X_t}dB_t =  dX_t - \kappa(\theta - X_t)dt  $$
%
% $$ \sqrt{X_t}dB_t = \frac{dX_t - \kappa(\theta - X_t)dt}{\sigma} $$
%
% $$Substituting\hspace{.1cm}the\hspace{.1cm}above\hspace{.1cm}value\hspace{.1cm}in\hspace{.1cm}equation\hspace{.1cm}of\hspace{.1cm}Z$$
% 
% $$ Z = \frac{\rho}{\sigma} \int_{0}^{T}dX_t\ -  \kappa(\theta - X_t)dt - \frac{\rho^2}{2} \int_{0}^{T}X_tdt\   $$
%
% $$ = \frac{\rho}{\sigma} \left(  \int_{0}^{T}dX_t\ - \int_{0}^{T} \kappa \theta dt\  + \int_{0}^{T} \kappa X_t dt\  \right) - \frac{\rho^2}{2} \int_{0}^{T} X_t dt\ $$
%
% $$ Z = \frac{\rho}{\sigma}   \left(  X_t - X_0 - \int_{0}^{T} \kappa \theta dt\  + \int_{0}^{T} \kappa X_t dt\  \right) -   \frac{\rho^2}{2} \int_{0}^{T} X_t dt\    $$
%
% $$Using\hspace{.1cm}the\hspace{.1cm}above\hspace{.1cm}equation\hspace{.1cm}to\hspace{.1cm}compute\hspace{.1cm}Z\hspace{.1cm}values$$
%
%
%
%%


%% Part 1: Standard MC


clear all;
close all;

S0 = 100;
r = 0.05;
T = 3/12;
K = 90:5:120;
X0 = 0.2;
k = 3;
theta = 0.2;
sigma = sqrt(2*theta*k);
dt = 1/365;
n = floor(T/dt);
rho_seq = [0, -0.3, -0.7, -0.9];
% Number of Monte Carlo samples
M = 20000;
% Number of Estimators/Trials
L = 100;


for rho_iter = 1:length(rho_seq) 

rho = rho_seq(rho_iter);
imp_vol_mc_est = zeros(7,100);
call_price_mc_est = zeros(7,100);

X = zeros(L,M,n);
Y = zeros(L,M,n);

for i = 1:L

X(i,:,1) = X0;

Y(i,:,1) = log(S0);

for t= 2:n
    dW = randn(M,1)*sqrt(dt);
    dB = randn(M,1)*sqrt(dt);
    X(i,:,t) = (1-k*dt)*X(i,:,t-1) + k*theta*dt + sigma*sqrt(max(X(i,:,t-1),0)).*dB';
    Y(i,:,t) = Y(i,:,t-1) + (r - .5*X(i,:,t-1))*dt + sqrt(max(X(i,:,t-1),0)).*(rho.*dB'+sqrt(1-rho^2)*dW');
end

S = exp(Y(i,:,n));
[Sv,Kv] = meshgrid(S,K);
D =exp(-r*T); 
C_mc = D*max(Sv(:,:)-Kv(:,:),0);
C_mc_avg = mean(C_mc,2);
implied_volatility_mc = blsimpv(S0,K',r,T,C_mc_avg);
call_price_mc = C_mc_avg;

imp_vol_mc_est(:,i) = implied_volatility_mc;
call_price_mc_est(:,i) = call_price_mc;
end



%% Part 2: Romano and Touzi
int_X = sum(X,3)*dt;
sigma_square = ((1-rho^2)/T).*int_X;
Z = (rho*(X(:,:,end)-X0 - k*theta*T+k*int_X))/sigma -(0.5*rho^2)*int_X;
C_bs = zeros(L,M,length(K));
S_new = S0*exp(Z);
sigma_x = sqrt(sigma_square);
for i = 1:length(K)

C_bs(:,:,i) = blsprice(S_new,K(i),r,T,sigma_x);

end

call_price_cond_mc_est = squeeze(mean(C_bs,2));
imp_vol_cond_mc_est = zeros(L,length(K));
for i = 1:length(K)
imp_vol_cond_mc_est(:,i) = blsimpv(S0,K(i),r,T,call_price_cond_mc_est(:,i));
end
call_price_cond_mc_est = call_price_cond_mc_est';
imp_vol_cond_mc_est = imp_vol_cond_mc_est';
%% Part 3 Heston price


Settle = datetime(2023,9,27);
Maturity = datemnth(Settle, 3);
c_heston = optByHestonNI(r,S0,Settle, Maturity,'call',K,X0,theta,k,sigma,rho,'DividendYield',0);
implied_volatility_heston = blsimpv(S0,K,r,T,c_heston');
implied_volatility_heston = implied_volatility_heston';
avg_mc = mean(call_price_mc_est,2);
avg_cond_mc = mean(call_price_cond_mc_est,2);
avg_mc_iv = mean(imp_vol_mc_est,2);
avg_cond_mc_iv = mean(imp_vol_cond_mc_est,2);
std_err_mc_call = std(call_price_mc_est,0,2);
std_err_cond_mc_call = std(call_price_cond_mc_est,0,2);
std_err_mc_imp_vol = std(imp_vol_mc_est,0,2);
std_err_cond_mc_imp_vol = std(imp_vol_cond_mc_est,0,2);
ratio_std_err_call = std_err_mc_call./std_err_cond_mc_call;
ratio_std_err_vol = std_err_mc_imp_vol./std_err_cond_mc_imp_vol;
s = "Rho is (" + num2str(rho) + ")";

disp(s);

T1 = table(K', c_heston, avg_mc, avg_cond_mc, std_err_mc_call, std_err_cond_mc_call, ratio_std_err_call,'VariableNames',{'strike', 'Heston Price', 'Avg. MC', 'Avg. C-MC', 'std. err. MC', 'std. err. C-MC', 'Ratio MC/C-MC' });
disp(T1);


T2 = table(K', implied_volatility_heston, avg_mc_iv, avg_cond_mc_iv, std_err_mc_imp_vol, std_err_cond_mc_imp_vol, ratio_std_err_vol,'VariableNames',{'strike', 'Heston Imp. Vol', 'Avg. MC', 'Avg. C-MC', 'std. err. MC', 'std. err. C-MC', 'Ratio MC/C-MC' });
disp(T2);
end

%% Answer for variance reduction
% Fo rho = 0, in case of standard monte carlo procedure a single random
% process is incorporated in the underlying process leading to highest
% underlying volatility compared to when rho between 0 and 1. When Rho is between 0 and
% 1, the division of two random processes results in a lower underlying
% volatility. 
% 1. For SMC, intial stock price is constant as rho changes. Therefore, increase in
% modulus of rho descreases the option price and implied volaitlity for out 
% of money options conversely for in the money options. 
% 2. For CMC, as rho increases the underlying volatility always decrease. 
% Also, as rho changes initial stock price changes in bls price. Therefore, increasing the standard
% deviation in the calculations of option prices and implied volatilities.
% 3. If initial stock price was constant, there would be significant
% reductions in the standard deviation. Now, since initial stock price is
% not constant the reduced standard deviation is increased due to this
% variance. 
% 4. Therfore as modulus of rho increases the standard error of
% calculations of option prices and implied volatilities decrease for SMC
% and CMC for out of money options and inversely for in the money options.
% But the ratio of standard errors decreases.
























