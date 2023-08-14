%% Part 1
X0 = 0.05;
k = 1;
theta = 0.05;
sigma =  sqrt(2*theta*k);
T = 10;
n = T/dt;
dt = 1/(T*365);
M = 100;

X = zeros(M, n);
X(:,1) = X0;

counter = 0; % initialize the counter to zero
X2 = zeros(M, n);
X2(:,1) = X0;

counter2 = 0; % initialize the counter to zero
X3 = zeros(M, n);
X3(:,1) = X0;

counter3 = 0; % initialize the counter to zero


%% Euler Maruyama
for t=2:n
    X(:,t) = (1-k*dt) * X(:,t-1) + k*theta*dt + sigma*sqrt(max(X(:,t-1),0)) .* randn(M,1) .* sqrt(dt);
    if any(X(:,t-1) < 0) % check if X(:,t-1) has any negative values
        counter = counter + sum(X(:, t-1) < 0); % increment the counter by 1
    end
end

average = counter/(M) % calculate the average number of times X(:,t-1) turns negative


%% Milstien Scheme

for t=2:n
    dW=randn(M,1) .* sqrt(dt);
    X2(:,t) = X2(:,t-1) + k*theta*dt- k*X2(:,t-1).*dt +...
        sigma*sqrt(max(X2(:,t-1),0)) .*dW  +...
        0.25*(sigma^2)*(dW.^2-dt);

    if any(X2(:,t-1) < 0) % check if X(:,t-1) has any negative values
        counter2 = counter2 + sum(X2(:, t-1) < 0); % increment the counter by 1
    end
end
c2 = sum(X2(:, t-1) < 0)
average2 = counter2/(M) % calculate the average number of times X(:,t-1) turns negative

%% Ito's reduction through implicit scheme

for t=2:n
    dW=randn(M,1)*sqrt(dt);
    X3(:,t) = X3(:,t-1) + ((0.5*k*theta-0.125*sigma^2)./sqrt(max(X3(:,t-1),0)))*dt' - 0.5*k*sqrt(max(X3(:,t-1),0))*dt + (sigma/2).*dW;
    X3(:,t) = X3(:,t).^2;
    X3(:,t) = X3(:,t-1) - k*X3(:,t-1).*dt + k*theta*dt  + sigma*sqrt(X3(:,t)).*dW  - 0.5*(sigma^2)*(dW.^2);

    if any(X3(:,t-1) < 0) % check if X(:,t-1) has any negative values
        counter3 = counter3 + sum(X3(:, t-1) < 0); % increment the counter by 1
    end
end
c3 = sum(X3(:, t-1) < 0);
average3 = counter3/(M); % calculate the average number of times X(:,t-1) turns negative


fprintf('The average number of times the schemes go negative are\n Euler Maruyama: %f,\n Milstien: %f,\n and Itos implicit reduction: %f\n', average, average2, average3);

