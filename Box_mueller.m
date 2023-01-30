%% Question 3
%
%
% Generate the pair of random variables $$ (X,Y) $$
%
% $$ (X,Y) = R(cos{\theta}, sin{\theta}) $$
%
% where we use a random number generator with poor lattice properties such as the generator
% $$ Xn+1 = (383Xn + 263) $$
% mod 10000 to generate our uniform random numbers. Use this
% generator together with the Box-Mueller algorithm to generate 5,000 pairs of independent
% random normal numbers. Plot the results. Do they appear independent?
%   
clear

% Setting parameters for the random sampling and number of samples as 5000
% as it generates normal distribution samples in pairs
a = 383;
b = 263;
m = 10000;
n = 5000;
x0 = 0.6;
y0 = 0.9;

% Preallocating for Normal samples as Z 
Z = zeros(n,2);

% Preallocating for u1 and u2
x = zeros(n,1);
y = zeros(n,1);

% Calculating the first sample 
x(1) = mod(a*x0+b,m);
y(1) = mod(a*y0+b,m);
u1 = x(1)/m;
u2 = y(1)/m;
z1 = sqrt(-2 * log(u1)) * cos(2 * pi * u2);
z2 = sqrt(-2 * log(u1)) * sin(2 * pi * u2);
Z(1,1) = z1;
Z(1,2) = z2;

% Loop to generate rest 4999 pairs of samples
for i = 2:n
    x(i) = mod(a*x(i-1)+b,m);
    y(i) = mod(a*y(i-1)+b,m);
    u1 = x(i)/m;
    u2 = y(i)/m;
    z1 = sqrt(-2 * log(u1)) * cos(2 * pi * u2);
    z2 = sqrt(-2 * log(u1)) * sin(2 * pi * u2);    
    Z(i,1) = z1;
    Z(i,2) = z2;
end

% Divide the Z vector of (5000,2) into X and Y
X = Z(:,1);
Y = Z(:,2);

% Plotting the samples
fig1 = figure(1);
plot(Z(:,1));
fig2 = figure(2);
plot(Z(:,2));

% Plotting the histogram of X and Y to show they are normal
fig3 = figure(3);
histogram(X);
fig4 = figure(4);
histogram(Y);


fprintf('The mean of X is %d \n', mean(X));
fprintf('The variance of X is %d \n', var(X));
fprintf('The mean of Y is %d \n', mean(Y));
fprintf('The variance of Y is %d \n', var(Y));

% For gauging independence calculating the correlation

corr = corrcoef(X,Y);
fprintf('The correlation between X and Y is %d', corr(1,2));

%% Independence
% As the approximate correlation is zero it implies that both X and Y are
% independent





