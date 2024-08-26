function T2=T2roi_RARE_Siemens_dcm(dname,roi,showPlot,roi_label,te_index)

if isa(roi,'char') || isa(roi,'cell')
    
roi=str2cell(roi);

else
    roi={roi};
end

d=dir(sprintf('%s/*.IMA',dname));
 ti=[];
im=[];
 for i=1:length(d)

 in=dicominfo([dname,'/',d(i).name]);
 ti(i)=in.EchoTime;

im(:,:,1,i)=dicomread([dname,'/',d(i).name]);

end


[ti,ind]=sort(ti);

im=im(:,:,1,ind);
if exist('te_index','var')
    ti=ti(te_index);
    im=im(:,:,1,te_index);
end


%%
%roi={'roi_csf_T1','roi_gm_T1','roi_wm_T1'};

options=optimset('MaxIter',20,'Display','off');

T2=[];
for i=1:length(roi)
if isa(roi{i},'char')
    m=ri(roi{i});

else
    m=roi{i};
end
 y=mean_roi(im,m);
 [beta,r]=lsqcurvefit(@exp_decay,[max(y),200],ti(:),y(:),[],[],options);
  if showPlot
  figure;
  plot(ti,y,'o');
  hold on;plot(ti,exp_decay(beta,ti),'r-');     
  set(gca,'YScale','log');
 
  if isa(roi{i},'char')
     title(sprintf('%s: T2 = %5.2f\n',strrep(roi{i},'_','-'),beta(2)));
 
  else
       title(sprintf('%s: T2 = %5.2f\n',roi_label,beta(2)));
  end
 end
 if isa(roi{i},'char')
     fprintf('%s: T2 = %f\n',strrep(roi{i},'_','-'),beta(2));
  
  else
       fprintf('%s: T2 = %f\n',roi_label,beta(2));
 end

 
 T2(i)=beta(2);
 ylim([min(y)*0.9,max(y)*1.05]);
end


