%% Q3
%
% $$ \frac{dS_t^l}{S_t^l} = (r-\delta_l)dt + \sigma_ldW_t^l $$
%  
% $$ for \hspace{0.1cm} l = 1,2,3... $$
% r = 0.05
% $\delta_l = 0.02$
% $\sigma_l = 0.3$
% 
% $\rho_{ll'} = 0.2$ 
% 
% $$ dW_t^ldW_t^{l'} = \rho_{ll'}dt$$
% 
% $$ \Bigl( \sum_{l=1}^{3} \lambda_T^l S_T^l - K  \Bigr)^{+} $$
% 
% $$ S_0^1 = S_0^2 = S_0^3 = K = 100 \hspace{0.1cm} T = 1$$
% 
% $$ \lambda_t^l = \frac{S_t^l}{S_t^1 + S_t^2 + S_t^3}  $$
% $$  N = 50,000 $$
% 
% $$ \Bigl( \sum_{l=1}^{3} \lambda_T^l S_T^l - K  \Bigr)^{+} $$
% 
%
%
%
%
%
%
%
%
%

%% least squares Monte Carlo
clear all
close all
r =0.05;
delta_l = 0.02;
sigma_l = 0.3;
rho_ll_dash = 0.2;
corr_l = [1 0.2 0.2; 0.2 1 0.2; 0.2 0.2 1];
mu = [0 0 0];
T=1;
N= 50000;
d=10;
%% Part 1: Monte Carlo European Basket option
K=100;
S0_l = K;
n=1;
dt=T/n;
dW_l = mvnrnd(mu, corr_l,N).*sqrt(dt);

ST_l = S0_l*exp((r-delta_l-.5*sigma_l^2)*dt+sigma_l.*dW_l);
lambda_l = ST_l./sum(ST_l,2);
payoff = max(sum(lambda_l.*ST_l,2)-K,0);
payoff_mc = mean(payoff)*exp(-r*T);


%% Part 2: Least square monte-carlo Bermudan option
dt = 1/12;
n = 13;
Sl = zeros(3,N,n);
Sl(:,:,1) = K;
lambda_l = zeros(3,N,n);
lambda_l(:,:,1) = Sl(:,:,1)./sum(Sl(:,:,1),1);
for i = 2:n
    dW_l = mvnrnd(mu, corr_l,N).*sqrt(dt);
    Sl(:,:,i)= Sl(:,:,i-1).*exp((r-delta_l-.5*sigma_l^2)*dt+sigma_l.*dW_l');
    lambda_l(:,:,i) = Sl(:,:,i)./sum(Sl(:,:,i),1);
end



%% LSMC price for different powers

LSMC_Prc = LSMC_Bermudan(Sl,K,exp(-r*dt), lambda_l);


%% LSMC ad-hoc

LSMC_Prc_adhoc = LSMC_Bermudan_adhoc(Sl,K,exp(-r*dt), lambda_l);

T1 = table(payoff_mc,LSMC_Prc_adhoc,LSMC_Prc,'VariableNames',{'Payoff Monte Carlo','LSMC Ad-hoc','LSMC multi-index'});
disp(T1);

%% LSMC Function
function LSMC_Prc = LSMC_Bermudan(Sl,K,D, lambda_l)
d=10;
p1 = repmat((0:d)',1,d+1,d+1);
p2 = permute(p1,[3,1,2]);
p3 = permute(p1,[2,3,1]);
p_all = p1+p2+p3;
index = p_all<=d;
multi_index = [p1(index),p2(index),p3(index)];

[~, n]=size(Sl(1,1,:));
aggregated_Sl = squeeze(sum(Sl.*lambda_l,1));
payoff = max(aggregated_Sl-K,0);

V = zeros(size(aggregated_Sl));
V(:,end) = payoff(:,end);

for t = n:-1:2
    V(:,t-1) = V(:,t)*D;
    ind = payoff(:,t-1)>0;
    psi = power(Sl(1,ind,t-1),multi_index(:,1)).*power(Sl(2,ind,t-1),multi_index(:,2)).*power(Sl(3,ind,t-1),multi_index(:,3));
    [Q, ~] = qr(psi',0);
    b = Q'*V(ind,t)*D;
%     b = pinv(psi)'*V(ind,t)*D;
    V(ind,t-1) = max(Q*b,payoff(ind,t-1));
end

LSMC_Prc = mean(V(:,1));
end


%% LSMC Function adhoc
function LSMC_Prc_adhoc = LSMC_Bermudan_adhoc(Sl,K,D, lambda_l)
d=10;

[~, n]=size(Sl(1,1,:));
aggregated_Sl = squeeze(sum(Sl.*lambda_l,1));
payoff = max(aggregated_Sl-K,0);

V = zeros(size(aggregated_Sl));
V(:,end) = payoff(:,end);

for t = n:-1:2
    V(:,t-1) = V(:,t)*D;
    ind = payoff(:,t-1)>0;
    psi = power(aggregated_Sl(ind,t-1),0:d);
    [Q, ~] = qr(psi,0);
    b = Q'*V(ind,t)*D;
%     b = pinv(psi)*V(ind,t)*D;
    V(ind,t-1) = max(Q*b,payoff(ind,t-1));
end

LSMC_Prc_adhoc = mean(V(:,1));
end


%% LSMC multi-index vs European
% The LSMC price applies to a Bermudan basket option that allows for early 
% exercise, resulting in a higher option price compared to a European basket 
% option that can only be exercised at maturity which is calculated using 
% the Monte Carlo simulation.
%% LSMC multi-index vs LSMC ad hoc
% In the ad hoc case, the coefficients of the polynomials are reduced as 
% the stock values are multiplied by lambdas, which are less than 1. As a 
% result, the Q*b term in the ad hoc case shrinks slightly in comparison to
% the LSMC multi-index method, where coefficients are calculated directly 
% based on the stock values raised to power combinations. Resulting in a
% slightly lower price compared to the LSMC multi index price.