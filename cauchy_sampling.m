%% Question 1:
% Generate 1000 daily \returns" Xi for i = 1; 2; : : : ; 1000 from each of the two distributions,
% the Cauchy and the logistic. Choose the parameters so that the median is zero and P(jXij <
% :06) = :95. Graph the total return over an n day period versus n. Is there a qualitative
% difference in the two graphs? Repeat with a graph of the daily return averaged over days
% 1; 2; : : : ; n.
%
% A Cauchy random variable has cdf
%
% $$ F(x) = \frac{1}{2} +
% \frac{1}{\pi}\tan^{-1}\left(\frac{x-x_0}{\gamma}\right) $$
%
% where $x_0$ is a location and $\gamma>0$ is a scale.
% 
% Logistic random variable has cdf
%
% $$ F(x) = \frac{1}{1+e^{-\left(\frac{x-\mu}{s}\right)}} $$
% 
% where  $\mu$  is a location and $s>0$ is a scale 


%% Answer
clear;
% Setting number of samples and scale parameters

n = [1:1000];
gamma = 0.0047;
s = -0.016377505;
U = rand(1000,1);
U_cauchy = U;
U_log = U;
% Generating Cauchy samples
% U_cauchy = rand(1000,1);
X_cauchy = tan((U_cauchy-0.5).*pi).*gamma;



% Generating Logistic samples
% U_log = rand(1000,1);
X_logistic = s*log((1./U_log)-1);



% Calculating the cumulative sum of the samples
X_cauchy_cumsum = cumsum(X_cauchy);
X_logistic_cumsum = cumsum(X_logistic);


fig1 = figure(1);
title('Cumulative returns vs n');
xlabel('n');
ylabel('Returns');
hold on
plot( n, X_logistic_cumsum);
plot( n, X_cauchy_cumsum);
legend('Logistic','Cauchy');
hold off

% Calculating Average of cumulative sum
X_avg_cum_log = zeros(1,1000);
X_avg_cum_cau = zeros(1,1000);

for i=1:1000
    X_avg_cum_log(i) = X_logistic_cumsum(i)/n(i);
    X_avg_cum_cau(i) = X_cauchy_cumsum(i)/n(i);
end

% Plotting the cumulative sum and average cumulative sum for both distributions
fig2 = figure(2);
title('Average Cumulative returns vs n');
xlabel('n');
ylabel('Returns');
hold on
plot(n, X_avg_cum_log);
plot(n, X_avg_cum_cau);
legend('Logistic','Cauchy');
hold off




















