function d=fl_fq_pc(dname,mask_factor,roi_name)
% fl_fq_pc(dname,mask_factor,roi_name)
% mask_factor: intensity threshold scale; threshold =
% mask_factor*max_intensity; default: 0.05.
% roi_name: can be a mask file/matrix or a dicom folder.  In the later, mask will
% be generated from 0.1*max(d(:)).  
% If not given, then assume the phase file is under dname/*_P_* and roi_name is under dname

if ~exist('mask_factor','var')
    mask_factor=0.1;
end

prefix=dname;
if ~exist('roi_name','var')

    roi_name=dname;
    dir_str=dir(fullfile(dname,'*_P_*'));
    dname=fullfile(dname,dir_str(1).name);
    
end
if isa(roi_name,'char')
    
    m=ri(roi_name);
    m=m(:,:,:,1)>mask_factor*max(m(:));
    figure;imshow4(double(m),[],[1,size(m,3)]);
    
else
    m=roi_name;
end


a=ri(dname);

a=double(a);
a_max=max(a(:));
a_min=min(a(:));

a=(a-(a_min+a_max)/2)/(a_max-a_min)*360;

x=[];
y=[];

for k=1:size(a,4)
    
    for sl=1:size(a,3)
        pos=ind2subb(size(m(:,:,sl)),find(m(:,:,sl)>0));
        a_tmp=a(:,:,sl,k);
        y=a_tmp(m(:,:,sl)>0);
        
        x=[pos,ones(size(pos,1),1),pos.^2];
        
        b=inv(x'*x)*x'*y;
        %b=x\y;
        
        %%{
        %figure;plot(y,x*b,'o');
        %hold on;plot([min(y),max(y)],[min(y),max(y)],'r-');
        %ylim([min(x*b),max(x*b)]);
        %}
        for i=1:size(a,1)
            for j=1:size(a,2)
                a(i,j,sl,k)=a(i,j,sl,k)-[i,j,1,i^2,j^2]*b;
            end
        end
    end
    
end
%%

if length(prefix)>4 && strcmp(prefix(end-3:end),'.mat')
    prefix=prefix(1:end-4);
else
    prefix=strtok(prefix,'.');
end

d=a;
try
[voxsize,center]=dcmDimCenter(dname);

save([prefix,'_detrend.mat'],'d','voxsize','center');

catch

save([prefix,'_detrend.mat'],'d');
    
    
end
