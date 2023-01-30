%% Question 4
%
% Use uniform samples to estimate the following integral by simulation:
%
% $$ {\int}{\int}  e^{\left(x+y\right)^{4}}dxdy $$ where 
% 
% $$ x {\in} (0,1) and y {\in} (0,1) $$
% 
%% Answer
clear
n_samples = [1000:1000:1000000];
integration_estimate = zeros(length(n_samples),1);
avg_integration_estimate = zeros(length(n_samples),1);

for j = [1:length(n_samples)]
    X = rand(n_samples(j),1);
    Y = rand(n_samples(j),1);
    integration_samples = exp((X+Y).^4);
    integration_estimate(j) = mean(integration_samples);
end

avg_integration_estimate_cumsum = cumsum(integration_estimate);

for j = [1: length(n_samples)]
    avg_integration_estimate(j) = avg_integration_estimate_cumsum(j)/j;
end




fig1 = figure(1);
% ax1 = axes('Parent', fig1);
title('Integration estimate vs Number of samples');
xlabel('Number of samples');
ylabel('Integration value');
hold on
plot(n_samples, integration_estimate);
plot(n_samples, avg_integration_estimate);
legend('Integration estimates','Averaged Estimates');
hold off

fprintf('The estimated value of integration is %d\n', avg_integration_estimate(length(n_samples)));


%% Actual integration value is



syms x y;

f = @(x,y) exp((x+y).^4);
q = integral2(f, 0,1,0,1);

fprintf('The Actual integrtaion value is %d \n', q)






