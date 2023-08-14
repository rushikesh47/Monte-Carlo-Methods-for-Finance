%% Jump diffusion process
clear all
close all
S0 = 100;
r = 0.05;
T = 1/12;
K = 90:120;
b = 0.1;
sigma = 0.3;

% Need range for these variables
lambda = [2,5];
a = [-0.1,0.1];


dt = 1/(4*365);
n = 20000;
dtn = ceil(T/dt);


vega = 0;

implied_volatility = zeros(4,length(K));
log_moneyness = log(K*exp(-r*T)/S0);
k = 1;
for i = 1:length(lambda)
    for j = 1:length(a)
        
        implied_volatility(k,:) = jump_diffusion_imp_vol(S0, dtn, n, r, sigma, T, dt, K, b, lambda(i),a(j));
        figure(10);
        subplot(length(lambda),length(a),k)
        plot(log_moneyness, implied_volatility(k,:))
        s = "Lambda is (" + num2str(lambda(i));
        s1 = "), a is (" + num2str(a(j)) + ")";
        title(s+s1)

        k = k+1;

    end
end



function implied_volatility = jump_diffusion_imp_vol(S0, dtn,n,r,sigma,T,dt,K, b, lambda, a)
X = log(S0)*ones(dtn,n);
J = ones(dtn,n);
f = 1;
for t = 2:dtn
    Z = normrnd(0,1,1,n);
%     Nt = pois_process(lambda,dt);
%     Nt = poissrnd(lambda*dt, 1,n);

%     logY=randn(1,Nt);
%     logY = randn(1);
%     logY=a+b*logY;
%     M=sqrt(Nt)*logY;
    Nt=poissrnd(lambda*dt,n,1);
    M = Nt*a + b*sqrt(Nt).*normrnd(0,1,n,1);
    vega = -1*(exp(a+(0.5*(b^2)))-1)*lambda;
%     if isempty(logY)
%         vega = 0;
%     else
%         vega = -1*(exp(a+(0.5*b^2))-1)*lambda;
%     end

    X(t,:)=X(t-1,:)+(r-.5*sigma^2+vega)*dt+sigma*sqrt(dt).*Z+M';
    J(t,:)=J(t-1,:).*exp(M');
    j_sum = sum(sum(J));

end

% X_bar = mean(X,1);
S = exp(X(dtn,:));
[Sv,Kv] = meshgrid(S,K);
D =exp(-r*T); 
% C_jd = D*mean(max(Sv(:,:)-Kv(:,:),0),2);
% figure(1)
% plot(K, C_jd);
po = D*max(Sv(:,:)-Kv(:,:),0);
po = mean(po, 2);
imp_vol = blsimpv(S0,K',r,T,po);
% figure(randi(100));
% plot(X(:,n));
% imp_vol = calcBSImpVol(1,po,S0,Kv',T,r,0);
% disp(imp_vol);
% imp_vol_bar = mean(imp_vol,2);
% log_moneyness = log(K.*exp(r*T)./S0);
% figure(2);
% plot(log_moneyness,imp_vol_bar);
% implied_volatility = imp_vol_bar;
implied_volatility = imp_vol;
end

%% generate Poisson process
function Nt=pois_process(lambda,t)

p=exp(-lambda*t);
F=p;
n=0;

U=rand;

while U>F
    n=n+1;
    p=p*lambda*t/n;
    F=F+p; 
end

Nt=n;
end

%% Effect of lambda
% The jump intensity parameter, λ, determines the frequency at which jumps
% occur in the stock price process. A higher λ value means that jumps occur
% more frequently, and this can lead to a higher implied volatility for
% out-of-the-money options. This is because the occurrence of jumps can
% increase the probability of extreme price movements, which in turn
% increases the implied volatility of out-of-the-money options.
%% Effect of jump size parameter 'a'
% The jump size parameter, a, determines the magnitude of the jumps in the
% stock price process. A higher a value means that jumps are larger, and 
% this can also lead to a higher implied volatility for out-of-the-money
% options. This is because larger jumps have a greater impact on the
% distribution of the stock price, and can lead to a higher implied
% volatility of out-of-the-money options. 
