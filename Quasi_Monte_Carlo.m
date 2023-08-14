%%
% $$\textbf{Problem\hspace{.05cm} 3:}$$
%
% $$\textbf{Estimate:}$$
%
% $$\theta = \mathrm{E}cos(||X||)$$ 
%
% $$where,\hspace{.1cm} X \in R^{100}\hspace{.1cm}normally\hspace{.1cm}distributed\hspace{.1cm}vector\hspace{.1cm}with,$$
%
% $$ X_i = \beta Z_0 + \sqrt{1-\beta^2} Z_i  $$
%
% $$where,\hspace{.1cm}Z_0\hspace{.1cm} is\hspace{.1cm} iid
% \hspace{.1cm}normal(0,1)$$
%
% $$ \beta \in (-1,1) $$
%
% $$\textbf{Estimation method:}$$
%
%
% $$ \hat{\theta} = \frac{1}{1000} \sum_{l=1}^{1000} cos(||X^l||)  $$
% 
% $$where,\hspace{.1cm} X^l \in
% R^{100}\hspace{.1cm}is\hspace{.1cm}a\hspace{.1cm}sample$$
%
% $$\textbf{Equation\hspace{.1cm}to\hspace{.1cm}reduce\hspace{.1cm}estimate\hspace{.1cm}using\hspace{.1cm}expectation\hspace{.1cm}of\hspace{.1cm}non-central\hspace{.1cm}chi-square:}$$
%
% $$ \theta =  \mathrm{E}\mathrm{E} \left[ cos \left(\sqrt{(1-\beta^2) \chi_{100}^2 (Z_0)}\right) | Z_0 \right]    $$
%
%
%
%%

%% Part 1 Standard Monte Carlo Simulation
d = 100;
n = 1000;
m = 1000;
theta_MC = zeros(m,2);
i = 1;
for beta  = [0, 0.8]

    
    for j = 1:1000
        z0 = randn(m,1);
        Zi = randn(m,d);
        X = zeros(n,d);
     
        X = beta*z0 + sqrt((1-beta.^2)).*Zi;
        
        norm_X = sqrt(sum(X.^2,2));
        Y = cos(norm_X);
        theta_MC(j,i) = mean(Y);

        
    end
    i=i+1;
end



%% Part 2 Monte Carlo Simulation with Chi - squared
dof = 100;

Y = zeros(1000);
theta_MC_2 = zeros(m,2);
i=1;
for beta =  [0, 0.8]
    
    for j = 1:1000
        z0 = randn(1000,1);
        nc = 100*beta.^2*z0.^2/(1-beta.^2);
        Y = cos(sqrt((1-beta.^2)*ncx2rnd(dof,nc,1000,1)));
        theta_MC_2(j,i) = mean(Y);
    end
    i=i+1;
end

%% Part 3 Quasi Monte Carlo Simulation


theta_RQMC_beta1 = zeros(m,1);
theta_non_RQMC_beta1 = zeros(m,1);


beta = 0;
p = sobolset(100,'Skip',1e3,'Leap',1e2);
for j = 1:1000
    p_scrambled = scramble(p,'MatousekAffineOwen');
    Sobol_unif_non_random = net(p,n);
    Sobol_unif_random = net(p_scrambled,n);
    Z_non_random = norminv(Sobol_unif_non_random);
    Z_random = norminv(Sobol_unif_random);
    Zi = Z_non_random;
    
    X = sqrt((1-beta.^2))*Zi;
    
    norm_X = sqrt(sum(X.^2,2));
    Y = cos(norm_X);
    theta_non_RQMC_beta1(j) = mean(Y);


    Zi = Z_random;
    
    
    
    X = sqrt((1-beta.^2))*Zi;
    

    norm_X = sqrt(sum(X.^2,2));
    Y = cos(norm_X);
    theta_RQMC_beta1(j) = mean(Y);
end



theta_RQMC_beta2 = zeros(m,1);
theta_non_RQMC_beta2 = zeros(m,1);


beta = 0.8;
p = sobolset(101,'Skip',1e3,'Leap',1e2);
for j = 1:1000
    p_scrambled = scramble(p,'MatousekAffineOwen');
    Sobol_unif_non_random = net(p,n);
    Sobol_unif_random = net(p_scrambled,n);
    Z_non_random = norminv(Sobol_unif_non_random);
    Z_random = norminv(Sobol_unif_random);
    Z0 = Z_non_random(:,1);
    Zi = Z_non_random(:,2:end);
    
    X = beta*Z0 + sqrt((1-beta.^2))*Zi;
    
    norm_X = sqrt(sum(X.^2,2));
    Y = cos(norm_X);
    theta_non_RQMC_beta2(j) = mean(Y);


    Z0 = Z_random(:,1);
    Zi = Z_random(:,2:end);
    
    
    
    X = beta*Z0 + sqrt((1-beta.^2))*Zi;
    

    norm_X = sqrt(sum(X.^2,2));
    Y = cos(norm_X);
    theta_RQMC_beta2(j) = mean(Y);
end

%% Plots

subplot(3,2,1);  
histfit(theta_MC(:,1),50)
title(sprintf('Std MC estimator beta=%g', 0))
hxl = xline(mean(theta_non_RQMC_beta1), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

subplot(3,2,2);
histfit(theta_MC(:,2),50)
title(sprintf('Std MC estimator beta=%g', 0.8))
hxl = xline(mean(theta_non_RQMC_beta2), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

subplot(3,2,3);  
histfit(theta_MC_2(:,1),50)
title(sprintf('1D MC estimatorbeta=%g', 0))
hxl = xline(mean(theta_non_RQMC_beta1), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

subplot(3,2,4);
histfit(theta_MC_2(:,2),50)
title(sprintf('1D MC estimatorbeta=%g', 0.8))
hxl = xline(mean(theta_non_RQMC_beta2), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

subplot(3,2,5);
histfit(theta_RQMC_beta1,50)
title(sprintf('RQMC estimator with beta=%g', 0))
hxl = xline(mean(theta_non_RQMC_beta1), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

subplot(3,2,6);
histfit(theta_RQMC_beta2,50)
title(sprintf('RQMC estimator with beta=%g', 0.8))
hxl = xline(mean(theta_non_RQMC_beta2), 'r', {'QMC'});
hxl.FontSize = 15;
hxl.LineWidth = 4;

  
%% Standard Error table

beta = [0,0.8];
error_table = table(beta',[std(theta_MC(:,1)),std(theta_MC(:,2))]',[std(theta_MC_2(:,1)),std(theta_MC_2(:,2))]',[std(theta_RQMC_beta1),std(theta_RQMC_beta2)]', 'VariableNames', {'Beta', 'Monte Carlo', '1-D Monte Carlo', 'RQMC'});
disp(error_table);


%% Unbiased Estimator:
% The standard Monte carlo and 1-Dimensional monte carlo appears to be
% unbiased as both are simulating the root equation with random sampling.
% Whereas the Random Quasi Monte Carlo Simulation is biased as it uses
% underlying sobol set grid and Matousek Affine Owen scrambling that can
% introduce bias in the estimators.

%% Least variance:
% The least variance is observed in RQMC. This could be due to the
% bias-variance trade off. And as the points are evenly disributed for this
% higher dimensionality problem.


