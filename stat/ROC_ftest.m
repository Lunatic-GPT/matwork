function tpr=ROC_ftest(dof1,dof2,lambda,p)
%[p,tpr]=ROC_ftest(dof1,dof2,lambda,pmax)


tpr = zeros(size(p));

for i=1:length(p)
    f = p2f(p(i),dof1,dof2);
    disp(i);    
    tpr(i) = ftest_nc(dof1,dof2,f,lambda);
end

figure;plot(p,tpr);
xlabel('p value');
ylabel('true positive rate');