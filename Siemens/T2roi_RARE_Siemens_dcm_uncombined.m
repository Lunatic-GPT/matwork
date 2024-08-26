function T2=T2roi_RARE_Siemens_dcm_uncombined(dname,roi,channels,roi_label,showPlot)

if isa(roi,'char') || isa(roi,'cell')
    
roi=str2cell(roi);

else
    roi={roi};
end

d=dir(sprintf('%s/*.IMA',dname));
 ti=[];
im=[];
ind=1;
 for i=1:length(d)

     if ~any((mod(i-1,32)+1)==channels)
         continue;
     end
 in=dicominfo([dname,'/',d(i).name]);
 ti(ind)=in.EchoTime;

im(:,:,:,ind)=dicomread([dname,'/',d(i).name]);
ind=ind+1;
 end


im=reshape(im,[size(im(:,:,1,1)),length(channels),size(im,4)/length(channels)]);

im=rms(im,3);
[ti,ind]=sort(ti);

im=im(:,:,1,ind);



%%
%roi={'roi_csf_T1','roi_gm_T1','roi_wm_T1'};

options=optimset('MaxIter',20,'Display','off');

T2=[];
for i=1:length(roi)
if isa(roi{i},'char')
    m=load(roi{i});
 m=m.roi;
else
    m=roi{i};
end
 y=mean_roi(im,m);
 [beta,r]=lsqcurvefit(@exp_decay,[max(y),mean(ti)],ti(:),y(:),[],[],options);
 if showPlot
  figure;
  plot(ti,y,'o');
  hold on;plot(ti,exp_decay(beta,ti),'r-');     
  set(gca,'YScale','log');
 
  if isa(roi{i},'char')
     title(sprintf('%s: T2 = %5.2f\n',roi{i},beta(2)));
 
  else
       title(sprintf('%s: T2 = %5.2f\n',roi_label,beta(2)));
  end
 end
 if isa(roi{i},'char')
     fprintf('%s: T2 = %f\n',roi{i},beta(2));
  
  else
       fprintf('%s: T2 = %f\n',roi_label,beta(2));
 end
  
 
 T2(i)=beta(2);
 ylim([min(y)*0.9,max(y)*1.05]);
end


