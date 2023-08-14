%% Given Data

M = 50;
N = 200;

%% Data Specific to Part 1

m = 6;
EX = 3.1;

%% Part 1(a)

a = (m:-1:1) - EX;
rho = roots(a);
realroots = imag(rho)==0;
rho = rho(realroots);
p = rho.^(1:m);
p = p/sum(p);

%% Part 1(b)

p_acc = zeros(M,6);
x_acc = zeros(M,N);
count = 1;
while count<=50
    draws = randi([1,6],1,N);
    test = abs(mean(draws) - 3.1);
    if test<0.01
        x_acc(count,:) = draws;
        [unqcnt, uniq] = hist(draws,unique(draws));
        p_acc(count,:) = unqcnt/200;
        count = count+1;
    end
end
p_comp = mean(p_acc,1);

%% Data specific to Part 2

q = [0.1,0.1,0.2,0.1,0.2,0.3];
EXq = 3.8;
EXq_2 = 16;
var = EXq_2 - EXq^2;

%% Part 2(a)

negE = @(x) sum(x.*log(x./q));
p_opt =fmincon(negE,q,[(1:m);(1:m).^2],[EXq;EXq_2],[ones(1,m)],1,zeros(1,m));

% p_opt = fmincon(negE, q, ((1:m)-EXq).^2,EXq_2-EXq^2,1:m,1,zeros(1,m) );

%% Part 2(b)

p_acc2 = zeros(M,6);
x_acc2 = zeros(M,N);
count2 = 1;
cumulative_probs = cumsum(q);
while count2<=50
    random_numbers = rand(1, N);
    dice_rolls = zeros(1, N);
    lower_bound = 0;
    for face = 1:6
        dice_rolls(((lower_bound < random_numbers).*(random_numbers<=cumulative_probs(face)))==1) = face;
        lower_bound = cumulative_probs(face);
    end
    test1 = mean(dice_rolls) - 3.8;
    test2 = mean(dice_rolls.^2) - 16;
    if test1<0 && test2<0
        x_acc2(count2,:) = dice_rolls;
        [unqcnt2, uniq2] = hist(dice_rolls,unique(dice_rolls));
        p_acc2(count2,:) = unqcnt2/200;
        count2 = count2+1;
    end
end
p_comp2 = mean(p_acc2,1);

%% Part 2(c)

cnt = zeros(M,1);
for j = 1:M
cnt(j) = 0;
for i = 1:N-1
    if x_acc2(j,i) == 3
        if x_acc2(j,i+1) == 6
            cnt(j) = cnt(j)+1;
        end
    end
end
end

pr_cal = sum(cnt)/(199*50);

%% 
Prob_calc = p_opt(3)*p_opt(6);