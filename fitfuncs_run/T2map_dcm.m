function T2map_dcm(dname,mask)
% T2map(fid_prefix[,ind_ex,stretch_exp])
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.

dname=num2str(dname);
d=dir(sprintf('%s/*.IMA',dname));
 ti=[];
ri(dname,1);
 for i=1:length(d)

 in=dicominfo([dname,'/',d(i).name]);
 ti(i)=in.EchoTime;

im(:,:,1,i)=dicomread([dname,'/',d(i).name]);

end

img2=ri(d,1);

[tmp,ind]=sort(ti);

im=im(:,:,1,ind);

sz=size(im);
if exist('mask','var')
    m=load(mask);
    m=m.roi;
else
   %m=ones(sz(1:3)); 
   m=im(:,:,1,1)>0.02*max(vec(im(:,:,1,1)));
end
T2 = zeros([sz(1:3),1]);

options=statset('FunValCheck','off','display','off');
for i=1:size(im,1)
    
    disp([i,size(im,1)]);
    for j=1:size(im,2)
        for k=1:size(im,3)
            
            y = double(squeeze(im(i,j,k,:)));
          
            if m(i,j,k)==0
                continue;
            end
            
%             if i==134 && j==161
%                 disp('');
%             else 
%                 continue;
%             end
            
         if length(ti)>2
            [beta,r]=nlinfit(ti(:),double(y(:)),@exp_decay,[max(y),mean(ti)],options); 
            if ~isnan(beta(1)) && ~isnan(beta(2))
             %ss = sum((y-mean(y)).^2);   
             T2(i,j,k,1) = beta(2);
 %            t2(i,j,k,2) = 1-sum(r.^2)/ss;
            end
         else
            T2(i,j,k)=1/log(y(1)/y(2))*(ti(2)-ti(1)); 
         end
        end
    end
end
T2(T2>10000)=0;

   name = ['T2map_',dname];
   nii=make_nii(T2);
    try
        orient=get_orient([dname,'.mat']);
        nii=nii_from_orient(nii,orient);
    catch      
        nii.untouch=1;
        nii.hdr.hist.magic='n+1';
    end
    save_untouch_niigz(nii,name);
    
%save(name,'T2');
%write_afni(t2(:,:,1,1),name);
   



