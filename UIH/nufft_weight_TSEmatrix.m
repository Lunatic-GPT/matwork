function dc=nufft_weight_TSEmatrix(m,w1,w2,w3)

% w: [x_l,y_u;x_r,y_d]
w=cat(3,w1,w2,w3);

m2=repmat(m*0,[1,1,3]);

for i=1:3
    area(i)=(diff(w(:,1,i))+1)*(diff(w(:,2,i))+1);

    sel_x=w(1,1,i):w(2,1,i);
    sel_y=w(1,2,i):w(2,2,i);

    n(i)=sum(vec(m(sel_y,sel_x)));
    m2(sel_y,sel_x,i)=1; 
end

den(1)=n(1)/area(1);
den(2:3)=diff(n)./diff(area);

dc=m*0;
dc(m2(:,:,1)==1)=1/den(1);

dc(m2(:,:,2)==1&m2(:,:,1)==0)=1/den(2);
dc(m>0&m2(:,:,2)==0)=1/den(3);

dc(m==0)=0;
