function T2map_dcm_multiScans(slist,mask)
% T2map_dcm_multiScans(slist,mask)
% ind_ex: subbriks to exclude.  1 based. default: [].
% stretch_exp: default false.

if isa(slist,'double')
    for i=1:length(slist)
      slist2{i}=num2str(slist(i));
    end
    slist=slist2;
end

 ti=[];

 
 for i=1:length(slist)
     
im(:,:,:,i)=ri(slist{i});
     
dir_str=dir(slist{i});

in=dicominfo(fullfile(slist{i},dir_str(end).name));
 ti(i)=in.EchoTime;


end


[ti,ind]=sort(ti);

im=im(:,:,:,ind);

sz=size(im);
if exist('mask','var')
    m=load(mask);
    m=m.roi;
else
    
   im0=im(:,:,:,1);
   
   m=im0>0.03*max(im0(:)); 
end
t2 = zeros([sz(1:3),2]);
options=statset('FunValCheck','off');
for i=1:size(im,1)
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
            
if length(ti)==2
    t2(i,j,k,1)=-(ti(2)-ti(1))/log(y(2)/y(1));
else
               [beta,r]=nlinfit(ti(:),double(y(:)),@exp_decay,[max(y),mean(ti)],options);
         
            if ~isnan(beta(1)) && ~isnan(beta(2))
             ss = sum((y-mean(y)).^2);   
             t2(i,j,k,1) = beta(2);
             t2(i,j,k,2) = 1-sum(r.^2)/ss;
            end
            
end
            
        end
    end
    disp(i);
end

   name = ['T2map_',slist{1}];

   save(name,'t2');
   



