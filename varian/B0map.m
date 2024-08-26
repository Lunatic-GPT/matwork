function B0map(fname,tau,ref)
%B0map(fname,tau,ref)

prefix = strtok(fname,'+');
[d,info]=BrikLoad(fname);
%tau = [0,50,100,200,300,400,600,800,1000,1500];

ph_ref = d(:,:,:,ref)-d(:,:,:,1);
nwrap = round(ph_ref/2/pi);
ph_ref = ph_ref-nwrap*2*pi;

ph=d(:,:,:,2:end)-repmat(d(:,:,:,1),[1,1,1,size(d,4)-1]);

tau_rel = tau(2:end)/tau(ref);
tau_rel = repmat(tau_rel,[size(d,1),1,size(d,2),size(d,3)]);
tau_rel=permute(tau_rel,[1,3,4,2]);
ph_pred=repmat(ph_ref,[1,1,1,size(d,4)-1]).*tau_rel;

nwrap = round((ph_pred-ph)/2/pi);

ph = ph+nwrap*2*pi;


WriteBrikEZ(ph,info,'',[prefix,'_relph']);

b = zeros(size(d,1),size(d,2),size(d,3));
for i=1:size(d,1)
    for j=1:size(d,2)
        for k=1:size(d,3)
            y = squeeze(ph(i,j,k,:));
            x = tau(2:end);
          b(i,j,k)=regress(y,x(:));
        end
    end
end

b=b/2/pi*1000000;
disp('done');
WriteBrikEZ(b,info,'B0map',[prefix,'_B0map']);


