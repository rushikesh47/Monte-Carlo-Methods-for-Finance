%% Q2
%
% $$ dZ_t = \frac{1}{\epsilon} (1 + \alpha - Z_t)dt + \sqrt{\frac{2Z_t}{\epsilon}} dW_t $$
%
% $$ \psi_k (z) = \frac{z^{-\alpha} e^{z}}{k!} \frac{d^k}{dz^k} (e^{-z}z^{k + \alpha})  $$
%
% $$ for \hspace{0.1cm} k = 0,1,2... $$
%
% $$ \mu (z) = \frac{i}{\Gamma (\alpha + 1)} z^{\alpha} e^{-z} $$
%
% $$ p_t(z|z_0) = \sum_{k=0}^{\inf} c_k e^{\frac{-kt}{\epsilon}} \psi_k (z_0) \psi_k(z) \mu (z) $$
%
% $$ c_k = \frac{ k! \Gamma ( \alpha + 1) }{ \Gamma (k+ \alpha + 1) } = \frac{k}{k+\alpha} c_{k-1} $$
%
% $$ c_0 = 1 $$
% 
% $$ \Delta t = \frac{1}{252} \epsilon = 1 $$
% 
% $$ \alpha = 4 $$ 
% 
% $$ P_{l,l'} \propto \sum_{k=0}^{20} c_k e^{\frac{-k \Delta t}{\epsilon}} \psi(k)(z^l) \psi_k(z^{l'}) $$ 
% 
% $$ \sum_{l' =1}^{ 1000} P_{l,l'} = 1   $$
% 
% $$ \hat F (Z_i) = \frac{\#{j|Z_j \leq Z_i}}{10^4}  $$
%
% $$ Z_{i+1} = Z_i + (\alpha - Z_{i+1})\frac{\Delta t}{\epsilon} + \sqrt{2Z_{i+1}} \frac{\Delta W_i}{\sqrt{\epsilon}} $$
%
% $$ \epsilon = 10^{-3}    $$
%
%%
%% CIR Process and Laguerre Polynomials
close all;
clear all;

dt = 1/252;
l=1000;
k=20;
alpha = 4;
eps_seq = [1,0.001];


for e = 1:length(eps_seq)
epsilon = eps_seq(e);
z_l = gamrnd(alpha+1,1,l,1);



c_k = zeros(k+2,1);
multi_k = ones(1000,1000,k+1);
c_k(1) = 1;
for i = 1:k+1
    c_k(i+1) = (i/(i+alpha))*c_k(i);
end

for i = 1:k+1
    laguerre_poly_l = gLaguerre(i-1,z_l,alpha);
    laguerre_poly_l_dash= gLaguerre(i-1,z_l',alpha);
    multi_k(:,:,i) = laguerre_poly_l*laguerre_poly_l_dash*c_k(i)*exp(-(i-1)*dt/epsilon);

%     multi_k(i) = c_k(i)*exp(((-i-1)*dt)/epsilon)*laguerre_poly_l*laguerre_poly_l_dash;
%     c_k(i+1) = ((i)/(i+alpha))*c_k(i);
end

Pll_dash = sum(multi_k,3);

Pll_dash_threshold = max(Pll_dash,0);
Pll_dash_final = Pll_dash_threshold./sum(Pll_dash_threshold,2);

initial_state = randsample(1000,1);

X = simulate(dtmc(Pll_dash_final),10000);
X_final = X(2:end);


[cnt_unique,uniq] = hist(X_final, unique(X_final));

count_mapping = [cnt_unique',uniq];

for i = 1:1000
    count_mapping(i,2) = z_l(i); 
end

count_sorted = sortrows(count_mapping,2);

emp_cdf = zeros(l,1);
for i = 1:l
    emp_cdf(i) = sum(count_sorted(1:i,1));
end
emp_cdf = emp_cdf/10000;



%% Part 2

Z_i = zeros(10000,1);
Z_i(1) = gamrnd(alpha+1,1,1,1);


for i = 1:9999
    dW_i = randn(1)*sqrt(dt);
%     Z_i(i+1) = Z_i(i) + (alpha-Z_i(i))*dt/epsilon + sqrt(2*max(Z_i(i),0))*dW_i/sqrt(epsilon);
    sqrt_Z_next = (sqrt(2)*dW_i/sqrt(epsilon) + sqrt((2*power(dW_i,2)/epsilon)+4*(1+dt/epsilon)*(Z_i(i)+alpha*dt/epsilon)))/(2*(1+dt/epsilon));
    Z_i(i+1) = power(sqrt_Z_next,2);
end

Z_i_sorted = sort(Z_i);
[Z_i_bucketed,edges] = discretize(Z_i_sorted,1000);
edges = edges(2:end);

emp_cdf_implicit = zeros(1000,1);

for j = 1:l
for i = 1:10000
    if Z_i_sorted(i)<edges(j)
        emp_cdf_implicit(j)  = emp_cdf_implicit(j)+1;
    end
end

end
emp_cdf_implicit = emp_cdf_implicit/10000;
figure(e)
plot(count_sorted(:,2),emp_cdf);
hold on;
x_gamma = 0:0.1:16;
y_gamma = gamcdf(x_gamma,alpha+1,1);
plot(x_gamma,y_gamma);
plot(edges,emp_cdf_implicit);
legend('Markov Chain', 'Gamma', 'Monte Carlo Implicit Scheme')
t = "For epsilon =" + num2str(epsilon);
title(t);
hold off;


process_markov = zeros(10000,1);
for i = 1:10000
    process_markov(i) = z_l(X_final(i));
end

figure(2*e+1)
plot(1:10000,process_markov);
hold on;
plot(1:10000,Z_i);
legend('Markov process', 'Implicit Scheme');
t = "Process for epsilon =" + num2str(epsilon);
title(t);
hold off;
end

%% Answer for comparison
% When epsilon is set to 1, the convergence of the Monte Carlo method is 
% rapid, as epsilon represents the time scaling for the implicit scheme. 
% As a result, all three closely track the gamma cumulative distribution 
% function (CDF).
%% Answer for epsilon=10^-3
% When epsilon is set to 10^-3, the convergence of the Monte Carlo implicit
% scheme is slower compared to epsilon=1, making it take longer to converge
% to the gamma cumulative distribution function (CDF). As a result, the 
% Markov chain more closely follows the gamma CDF compared to the implicit 
% scheme.