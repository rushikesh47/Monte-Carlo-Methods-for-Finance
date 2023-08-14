%% Part 1(a) 
close all
d = 100;

sigma = zeros(d+2,d+2);
T = [1 1.5 2];
sigma(2:end-1,2:end-1) = sign(randn(d,d));
b = 1/T(2) ;
[row,col]=ind2sub([d,d],1:d^2);
T = max(500000,d^2);
for t=1:T
    i = 1+mod(t-1,d^2);
    ir = row(i):(row(i)+2);
    ic =col(i):(col(i)+2);
    E = sum(sum(sigma(ir,ic).*[0 1 0;1 0 1;0 1 0]));
    p = 1/(1+exp(-2*b*E));
    U=p-1+rand;
    sigma(1+row(i),1+col(i))=sign(U);
end
%%
figure,imagesc(sigma)
colormap(gray)
saveas(gcf, 'isingmodel.eps', 'epsc')

%% Part 1(b)
%%GibbsSamplingfromHMMPosterior
g=2;
Y=sigma;
Y(2:end-1,2:end-1)=Y(2:end-1,2:end-1)+g*randn(d,d);
sigpost=Y;
for t=1:T
i=1+mod(t-1,d^2);
ir=row(i):(row(i)+2);
ic=col(i):(col(i)+2);
E=sum(sum(sigpost(ir,ic).*[0,1,0;1,0,1;0,1,0]));
p=1/(1+exp(-2*b*E));
Yi=Y(1+row(i),1+col(i));
post=normpdf([Yi+1,Yi-1],0,g).*[1-p,p];
post=post/sum(post);
U=post(2)-1+rand;
sigpost(1+row(i),1+col(i))=sign(U);
end
errind=sigma-sigpost~=0;
err=sum(errind(:))/(d^2);
fprintf('prcntg increct: %1.2f\n',100*err)



figure,imagesc(Y)
colormap(gray)
saveas(gcf, 'isingmodel.eps', 'epsc')

%% Part 1(c)
figure,imagesc(sigpost)
colormap(gray)
saveas(gcf, 'isingmodel.eps', 'epsc')

%% Part 1(d)

sigma_icm = Y;
% T = [1 1.5 2];
% b = 1/T(1) ;
[row,col]=ind2sub([d,d],1:d^2);
% T = max(500000,d^2);
T = 500000;
err= zeros(T,1);
for t=1:T
    i=1+mod(t-1,d^2);
    ir=row(i):(row(i)+2);
    ic=col(i):(col(i)+2);
    E=sum(sum(sigma_icm(ir,ic).*[0,1,0;1,0,1;0,1,0]));
    p=1/(1+exp(-2*b*E));
    Yi=Y(1+row(i),1+col(i));
    post=normpdf([Yi+1,Yi-1],0,g).*[1-p,p];
    [~, idx_post] = max(post);
    sigma_icm(1+row(i),1+col(i))=sign(idx_post-1.5);


% 
%     i = 1+mod(t-1,d^2);
%     ir = row(i):(row(i)+2);
%     ic =col(i):(col(i)+2);
%     E = sum(sum(sigma_icm(ir,ic).*[0 1 0;1 0 1;0 1 0]));
%     p = 1/(1+exp(-2*b*E));
%     U=p-1+rand;
%     sigma_icm(1+row(i),1+col(i))=sign(U);
    errind=sigma-sigma_icm~=0;
    err(t)=sum(errind(:))/(d^2);
    
end
%%
figure,imagesc(sigma_icm)
colormap(gray)
saveas(gcf, 'isingmodel.eps', 'epsc')

figure(5);
plot(err);
