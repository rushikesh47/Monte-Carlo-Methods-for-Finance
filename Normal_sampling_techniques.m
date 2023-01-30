%% Question 5
%
%
% Assume that you can sample from uniform distribution on (0,1). Compare and report the
% computational time using following four methods by generating N = 100,000,000 indepen-
% dent samples from a standard normal distribution: a) Box-Muller; b) Marsaglia's polar
% method; c) rational approximation (for inverse of cdf); d) acceptance-rejection. Do not use
% vectorization, simply run a loop from 1 to N for each method. Provide the screenshot of
% your code.
%
%


%%
% clear the environment
clear
% initialize the number of random numbers to generate
N = 100000000
% preallocate vectors to store the results
X_boxmuller = zeros(1, N); 
X_bray = zeros(1, N);
X_rational = zeros(1, N); 
X_ar = zeros(1,N);

%% Box-Muller method

tic;
% generate N/2 pairs of normally distributed random numbers
for i = 1:N/2
    U1 = rand; %uniform random numbers
    U2 = rand;  % uniform random number
    sqrtln = sqrt(-2*log(U1)); %intermediate calculation
    X_boxmuller(2*(i-1)+1) = sqrtln*cos(2*pi*U2); %first normally distributed random number
    X_boxmuller(2*(i-1)+2) = sqrtln*sin(2*pi*U2); %second normally distributed random number
end
fprintf('Box Mullers ');
toc;

%% Bray_Marsaglia
tic;
% generate N/2 pairs of normally distributed random numbers
for i = 1:N/2
    while 1>0
        V1 = -1 + 2*rand;  % uniform random number
        V2 = -1 + 2*rand;   % uniform random number
        S = V1*V1 + V2*V2;  % intermediate calculation
        if S<=1  % check if the pair is inside the unit circle
            sqrtln = sqrt(-2*log(S)/S);  
            X_bray(2*(i-1)+1) = V1*sqrtln;  % first normally distributed random number
            X_bray(2*(i-1)+2) = V2*sqrtln;  % second normally distributed random number
            break  % exit the loop when a pair is found
        end
    end
end
fprintf('Bray_Marsaglias ');
toc;


%% Sampling Normals with Rational Approximation
tic;

%coeffs used for the method from  Beasley-Springer-Moro Algorithm
a0=2.50662823884;
a1=-18.61500062529;
a2=41.39119773534;
a3=-25.44106049637;

b0=-8.47351093090;
b1=23.08336743743;
b2=-21.06224101826;
b3=3.13082909833;

c0=0.3374754822726147;
c1=0.9761690190917186;
c2=0.1607979714918209;
c3=0.0276438810333863;
c4=0.0038405729373609;
c5=0.0003951896511919;
c6=0.0000321767881768;
c7=0.0000002888167364;
c8=0.0000003960315187;

% Loop to generate N samples of standard normal distribution using rational approximation
for i = 1:N
    U = rand();  % Generate a uniform random number
    Y = U - 0.5;
    if abs(Y) < 0.42
        r = Y^2;
        % Apply the rational approximation formula to generate a sample of standard normal
        X_rational(i) = Y * (a0 + r * (a1 + r * (a2 + r * a3))) / (1 + r * (b0 + r * (b1 + r * (b2 + r * b3))));
    else
        if Y > 0
            r = 1 - U;
        else
            r = U;
        end
        r = log(-log(r));
        % Apply the rational approximation formula to generate a sample of standard normal
        X_rational(i) = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
        if Y < 0
            X_rational(i) = -X_rational(i);
        end
    end
end
fprintf('Rational Approximations ');
toc;


%% Accept-Reject for Normal Samples
tic;
% Constant used in the acceptance-rejection method 
c = exp(.5);
% Loop to generate N samples of standard normal distribution
for i = 1:N
    while 1 > 0
        U1 = rand;
        if U1 < .5
            Y = log(2*U1);
        else
            Y = -log(2*(1-U1));
        end
        f = exp(-.5*Y^2)/sqrt(2*pi);
        g = .5*exp(-abs(Y));
        U2 = rand;
        if U2 < f/(c*g)
            X_ar(i) = Y;
            break
        end
    end
end
fprintf('Acceptance/Rejections ');
toc;

%% Histograms
figure;
subplot(2,2,1);
histogram(X_boxmuller);
title('Box Muller');

subplot(2,2,2);
histogram(X_bray);
title('Bray Marsaglia');

subplot(2,2,3);
histogram(X_rational);
title('Rational Approx');

subplot(2,2,4);
histogram(X_ar);
title('Acceptance/Rejection');

