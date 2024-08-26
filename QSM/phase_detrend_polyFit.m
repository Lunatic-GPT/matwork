function a=phase_detrend_polyFit(dname,roi_name)

%a=phase_detrend_polyFit(dname[,roi_name])
% roi_name: can be a mask file/matrix or a dicom folder.  In the later, mask will
% % If not given, all voxels will be considered.



if isa(dname,'char')
    a=ri(dname);
else
    a=dname;
end


if exist('roi_name','var')
    if isa(roi_name,'char')
        m=ri(roi_name);
    else
        m=roi_name;
    end
else
    
    m=ones(size(a));
   % figure;imshow4(double(m),[],[1,size(m,3)]);
end

x=[];
y=[];

for k=1:size(a,4)
    
    for sl=1:size(a,3)
        pos=ind2subb(size(m(:,:,sl)),find(m(:,:,sl)>0));
        
        
        
        a_tmp=a(:,:,sl,k);
        y=a_tmp(m(:,:,sl)>0);
        pos=pos./repmat(size(m(:,:,1)),[size(pos,1),1]);
        x=[pos,ones(size(pos,1),1),pos.^2];
     %   norm=mean(x,1);
      %  x=x./repmat(norm,[size(pos,1),1]);
        b=x\y;
        %b=x\y;
        
        %%{
        %figure;plot(y,x*b,'o');
        %hold on;plot([min(y),max(y)],[min(y),max(y)],'r-');
        %ylim([min(x*b),max(x*b)]);
        %}
        sz=size(m(:,:,1));
               
        for i=1:size(a,1)
            for j=1:size(a,2)
               % disp([i,j]);
                 pos=[i,j]./sz;
                 x=[pos,1,pos.^2];
        
                a(i,j,sl,k)=a(i,j,sl,k)-x*b;
            end
        end
        a(:,:,sl,k)=a(:,:,sl,k);
    end
    
end

%%

if isa(dname,'char')
    
prefix=strtok(dname,'.');
    if strcmp(prefix(end-3:end),'.mat')
        prefix=prefix(1:end-4);
    else
        prefix=strtok(prefix,'.');
    end
    
    save([prefix,'_detrend.mat'],'a');
end
