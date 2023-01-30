%% Question 2
% Suppose X1,.....,Xn are i.i.d with cdf 
%
% $$F(x) =\frac{1}{2}+\frac{1}{\pi}\tan^{-1}\left(x\right) $$
%
% Denote the order statistics by 
%
% $$ X_{(1)} \leq X_{(2)} \leq \dots \leq X_{(n)} $$
%
% Take n=100
%% Part 1
%
% Generate 1000 realizations of (X1,.....,X100) and use these samples to
% estimate     
% 
% $$ EX_{(i)} $$        
%
% for i = 30; 50; 70.
%
clear;

% Setting parameters for number of realizations, number of samples and
% order statistics required
N = 1000;
n = 100;
i = [30,50,70];

% Generating the samples, sorting and then calculating the mean of all the order
% statistics
U = rand(n,N);
X = tan(pi*(U-.5));
X = sort(X);
EXi = mean(X.');

% Printing the order statistics for 30th, 50th and 70th order statistic
for j = [1:3]
    fprintf('Method 1: EX[%d] = %1.6f\n',i(j),EXi(i(j)));
end


%% Part 2
%
% Show the density of X(i) is of the form 
% 
% $$ f(z) = cz^{a-1}(1-z)^{b-1} $$   for  $$ z {\in} (0, 1)  $$  and c a normalizing constant.
% Proof at the bottom
%% Part 3
%
% Use Matlab's beta distribution sampler to generate 1000 samples of  $$ X_{(1)} $$. 
% 
% Use these samples to estimate $$ EX_{(i)} $$ for i = 30, 50, 70.
%

% Generating 100 samples for required statistics and calculating the mean
for j =[1:3]
    U = betarnd(i(j),n+1-i(j),1,N);
    X = tan(pi*(U-.5));
    EXi = mean(X);
    fprintf('Method 2: EX[%d] = %1.6f\n',i(j),EXi);
end

%% Part 4
%
% What is the advantage gained by sampling for a beta rather than each time generating (X1,....,X100)?
%
%


% Setting parameters for the experiment 
N = 1000;
n = 100;
i = 50;

% Estimate the time required for order statistic 
tic


U = rand(n,N);
X = tan(pi*(U-.5));
X = sort(X);
EXi = mean(X);

toc

% Estimating the time required using beta distribution to sample required
% order statistic
tic

U = betarnd(i,n+1-i,1,N);
X = tan(pi*(U-.5));
EXi = mean(X);

toc

%% Advantage:1
% The advantage of using Beta random sampling is the number of samples
% needed are substantially lesser than the former method. Specifically, for
% 50th order statistic in a sample of 100 values we'll need to take only
% 1000 realizations. In contrast, for the former method to get the sample
% for a certain order statistic we'll need to take 1000x100 samples and
% sort them. Hence, it takes lesser time and compuattional power to sample using Beta random variable



   
