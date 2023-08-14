%% Question 3: Part 2
close all;
clear all;
load HW5.mat;

%% Parameter estimation using middle noisy path



d=100;
sigma = zeros(d+2,d+2);
sigma(2:end-1,2:end-1) = Noisy_Binary_Wedding(221:320,131:230);
Y =sigma;
% T_g_initial_guess = [0.5 1;2.5 1;0.5 0.2;2.5 0.2];
% T_g_initial_guess = [0.5	0.5;1	0.5;1.5	0.5;2	0.5;2.5	0.5;0.5	1;1	1;1.5	1;2	1;2.5	1;0.5	2.5;1	2.5;1.5	2.5;2	2.5;2.5	2.5;0.5	3;1	3;1.5	3;2	3;2.5	3;0.5	5;1	5;1.5	5;2	5;2.5	5];
T_g_initial_guess = [0.5 0.5;1 0.5;1.5 0.5;2 0.5;2.5 0.5;2.5 1;2.5 2.5;2.5 3;2.5 5];
[row,col]=ind2sub([d,d],1:d^2);
T = max(1000000,d^2);
T_est = zeros(length(T_g_initial_guess),1);
g_est = zeros(length(T_g_initial_guess),1);


sigma_original = zeros(d+2,d+2);
sigma_original(2:end-1,2:end-1) = Binary_Wedding(221:320,131:230);

figure
imagesc(sigma_original)
colormap(gray)
title('Original middle patch')

figure
imagesc(Y)
colormap(gray)
title('Noisy middle patch')

for ini_index = 1: length(T_g_initial_guess)
% ini_index = 1;
 % Buffer present thats why 
% ini_index = 3;
sigma_post = Y;


g=T_g_initial_guess(ini_index,2);
b=1/T_g_initial_guess(ini_index,1);
b_prev = b;
g_prev = g;

fprintf("For initial guess T=" + num2str(T_g_initial_guess(ini_index,1))+ " g="+ num2str(T_g_initial_guess(ini_index,2)))
fprintf('\niter, T_est, g_est\n');
flag=true;
iter=0;
while flag
   
    C=0;
    res=0;
    k=0;
    for t=1:T
        i = 1+mod(t-1, d^2);
        ir = row(i):(row(i)+2);
        ic =col(i):(col(i)+2);
        E = sum(sum(sigma_post(ir,ic).*[0 1 0;1 0 1;0 1 0]));
        p = 1/(1+exp(-2*b*E));
        
        Yi = Y(1+row(i), 1+col(i));
        post = normpdf([Yi+1, Yi-1],0,g).*[1-p,p];
        post = post/sum(post);
        U = post(2)-1+rand;
        sigma_post(1+row(i), 1+col(i)) = sign(U);
        if t > 50*d^2 && mod(t/(d^2),1)==0
            res = res + mean((Y(:)-sigma_post(:)).^2);
            C = C + sum(sum(sigma_post(2:end-2,:).*sigma_post(3:end-1,:)))...
                + sum(sum(sigma_post(:, 2:end-2).*sigma_post(:,3:end-1)));
            k=k+1;
        end
    end
    g = sqrt(res/k);
    E1_iter = C/k;
    [~,idx] = min(abs(E1(1,:)-E1_iter));
    b = 1/(0.5+0.01*(idx-1));
    iter=iter+1;
    fprintf('%d %f %f\n',iter, 1/b,g );
    if (abs(g-g_prev)<0.02) && (b_prev==b) 
        flag = false;
        T_est(ini_index,1) = 1/b;
        g_est(ini_index,1) = g;
    else
        g_prev = g;
        b_prev = b;
    end
end
sig_restored=Y;
for t=1:T
    i=1+mod(t-1,d^2);
    ir=row(i):(row(i)+2);
    ic=col(i):(col(i)+2);
    E=sum(sum(sig_restored(ir,ic).*[0,1,0;1,0,1;0,1,0]));
    p=1/(1+exp(-2*b*E));
    Yi=Y(1+row(i),1+col(i));
    post=normpdf([Yi+1,Yi-1],0,g).*[1-p,p];
    post=post/sum(post);
    U=post(2)-1+rand;
    sig_restored(1+row(i),1+col(i))=sign(U);
end
% errind=sigma-sig_restored~=0;
% err=sum(errind(:))/(d^2);
% fprintf('prcntg increct: %1.2f\n',100*err)

%%
figure
% subplot(1,3,1)
% imagesc(sigma)
% colormap(gray)
% saveas(gcf, 'isingmodel_original', 'epsc')

% subplot(1,2,1)
% imagesc(Y)
% colormap(gray)
% saveas(gcf, 'isingmodel_noisy', 'epsc')

% subplot(1,2,2)
imagesc(sig_restored)
colormap(gray)
image_title = "Isingmodel restored with initial guess T=" + num2str(T_g_initial_guess(ini_index,1))+" g="+num2str(T_g_initial_guess(ini_index,2));
title(image_title);
% saveas(gcf, 'isingmodel', 'epsc')
end



%% Restoring full image using estimated parameters

[len, width] = size(Noisy_Binary_Wedding);
[row,col]=ind2sub([len, width],1:len*width);
sigma_full = zeros(len+2,width+2);
sigma_full(2:end-1,2:end-1) = Noisy_Binary_Wedding(:,:);
T = max(100*len*width,len*width);
Y =sigma_full;
sig_restored_full_image=Y;
b=1/mean(T_est);
g=mean(g_est);
for t=1:T
    i=1+mod(t-1,len*width);
    ir=row(i):(row(i)+2);
    ic=col(i):(col(i)+2);
    E=sum(sum(sig_restored_full_image(ir,ic).*[0,1,0;1,0,1;0,1,0]));
    p=1/(1+exp(-2*b*E));
    Yi=Y(1+row(i),1+col(i));
    post=normpdf([Yi+1,Yi-1],0,g).*[1-p,p];
    post=post/sum(post);
    U=post(2)-1+rand;
    sig_restored_full_image(1+row(i),1+col(i))=sign(U);
end

fprintf('The estimated T value is %f and estimated standard error is %f', mean(T_est), mean(g_est))


figure
imshow(Binary_Wedding)
colormap(gray)
title('Original Picture')
figure
imshow(Noisy_Binary_Wedding)
colormap(gray)
title('Noisy Picture')
figure
imshow(sig_restored_full_image)
colormap(gray)
title('Restored Picture')


%%

% The initial guesses of the T (entropy variable) and the g (standard error
% of the noise) have an impact on the convergence of the parameters in the
% model. The reason being the existence of a global and local maximas. It
% is not apparent and the difference between the two maximas is not
% distinguishable in this model. 
% 
% The estimates converge to T = 0.51 when the range of value of T and g
% satisfy both T<1 and g>1, which was observed by many simulations using a
% lot of different combinations, almost 30 combinations were used to check this. 
%
% The other minima occurs where the above two conditions not both occur,
% and the values converge at approximately T = 1.88 and g = 1.45. And we
% have used this combination as most of the initial guesses converge at
% this combination, for the restoration of the noisy wedding image.

%%