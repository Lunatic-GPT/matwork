function fitDCE(d,ainp,TR,b0,mask,prefix)
%fitDCE(d,ainp,TR,b0,mask,prefix)


a=BrikLoad(d);

sz=size(a);
if isempty(mask)
    mask=ones(sz(1:3));
end

fit=@(b,x) exp_decay_conv(b,x,TR);

options=optimset('MaxFunEvals',40,'Display','off');
b=zeros([sz(1:3),4]);
for i=1:size(a,1)
    disp(i);
    for j=1:size(a,2)
        for k=1:size(a,3)
            y=squeeze(a(i,j,k,:));
            beta=lsqcurvefit(fit,b0,ainp,y,[],[],options);
            b(i,j,k,1)=beta(1);
            b(i,j,k,2)=1/beta(2);
            b(i,j,k,3)=beta(1)*beta(2);
            mx=max(abs(fit(beta,ainp)));
            b(i,j,k,4)=mx*sign(beta(1));
             
            b(i,j,k,5)=1-mean((fit(beta,ainp)-y).^2)/mean((y-mean(y)).^2);
             
             
         
        end
    end
end


write_afni(b,prefix);