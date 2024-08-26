function path_to_mask(fmag, fmask, sig_scale)
% data_fname: the nifti file name
% sig_scale: if <0; then assume vessels are hypointense.  default 2.


if ~exist('sig_scale','var')
    sig_scale=-2;
end

nii=load_untouch_niigz(fmag);
d=double(nii.img);

m=ri(fmask);

sz=size(d);
roi=zeros(sz,'single');

m_background=zeros(sz,'single');

m=clusterize2(m,2);
tic;
for j=1:max(m(:))
    
        [all,ring]=m_ext_calc(m==j,1);
       
        d_sort=sort(d(ring>0));
        
       
        
        if sig_scale>0
            dtmp=d_sort(1:round(end*0.9));
            mn=mean(dtmp);
            sd=std(dtmp);
        
            mtmp=d>(mn+sig_scale*sd)&all>0;
            
        else
            dtmp=d_sort(round(end*0.1)+1:end);
            mn=mean(dtmp);
            sd=std(dtmp);
            mtmp=d<(mn+sig_scale*sd)&all>0;
        end
        mtmp=clusterize2(mtmp,2);
        roi(mtmp==1)=j;       
        time_left(j,max(m(:)),toc);

end
%roi=single(m==1);

roi=single(roi>0);
%roi_path=single(roi_path);
nii.img=roi;
prefix=strtok(filename(fmask),'.');
prefix=[prefix,'_sig',num2str(sig_scale)];
%prefix='DMV_biggest';
save_untouch_niigz(nii,prefix);

%nii.img=all;
%prefix='DMV_biggest_dilate1';
%save_untouch_niigz(nii,prefix);

%nii.img=roi;
%prefix=strtok(ftrace,'.');
%save_untouch_niigz(nii,prefix);

function [m,ring]=m_ext_calc(c,thk)
% m: the mask after expanding roi by thk voxels
% ring: the mask after excluding the original roi

ind=find(c>0);
l=2;
ind2=zeros(length(ind)*(2*l+1)^3,1);
count=0;
sz=size(c);
for i=1:length(ind)
    
    ijk=ind2subb(sz,ind(i));
    
   
    
    for i1=-thk:thk
        for i2=-thk:thk
            for i3=-thk:thk   
                 if ijk(1)+i1<1 || ijk(1)+i1>sz(1) || ijk(2)-i2<1 || ijk(2)+i2>sz(2) ||ijk(3)-i3<1 || ijk(3)+i3>sz(3)
                    continue;
                 end
    
                count=count+1;
                ind2(count)=ind(i)+i1+i2*sz(1)+i3*(sz(1)*sz(2));
                
                subtmp=ind2subb(sz,ind2(count));
                if ~any(subtmp~=[158,162,89])
                   disp(''); 
                end
            end
        end
    end
    
end

ind2(ind2==0)=[];
m=c*0;
m(ind2)=1;

ind2=setdiff(ind2,ind);


ring=c*0;
ring(ind2)=1;

function m=neighbors_mask(pos,dist,sz,dim)

m=zeros(sz);
for i=1:size(pos,1)
    
    for j1=-ceil(dist/dim(1)):ceil(dist/dim(1))
        for j2=-ceil(dist/dim(2)):ceil(dist/dim(2))
            for j3=-ceil(dist/dim(3)):ceil(dist/dim(3))
                
                
                ind=round(pos(i,:)./dim)+[j1,j2,j3];
                
                if any(ind<1) || any(ind>sz)
                    continue;
                end
                
                if m(ind(1),ind(2),ind(3))==1
                    continue;
                end
                
                if sos([j1,j2,j3].*dim)<=dist
                    m(ind(1),ind(2),ind(3))=1;
                end
            end
        end
    end
end
