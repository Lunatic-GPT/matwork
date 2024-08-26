function neurite_trace_to_mask(fimage,ftrace,sig_scale)
% data_fname: the nifti file name
% sig_scale: if <0; then assume vessels are hypointense.  default 2.
% if sig_scale is empty, then no extra voxels will be added to the mask,
% only the voxels that intersect with the path will be included in the ROI.
% assuming data_fname is the same file used for drawing the traces and that
% (x,y,z) in the trace file corresponds to the first, second, and third dimensions of the
% image.

if ~exist('sig_scale','var')
    sig_scale=2;
end


nii=load_untouch_niigz(fimage);
dm=nii.hdr.dime.pixdim(2:4);
d=double(nii.img);

sz=size(d);
roi=zeros(sz,'single');
roi_path=zeros(sz,'single');

m_background=zeros(sz,'single');

pos_all=read_swc(ftrace,dm);

for j=1:length(pos_all)
    
    pos=pos_all{j};
    pos=round(connect_pos(pos));
  
    pos=pos.*repmat(dm,[size(pos,1),1]);
    
    if ~isempty(sig_scale)
        me=neighbors_mask(pos,1.2,sz,dm);
        men=neighbors_mask(pos,0.8,sz,dm);
             
        me=me&~men;
        
        mn=mean(d(me));
        sd=std(d(me));
        
        if sig_scale>0
            roi(d>(mn+sig_scale*sd) & men)=j;
        else
            roi(d<(mn+sig_scale*sd) & men)=j;
        end
    end

    for i=1:size(pos,1)
        id=round(pos(i,:)./dm);

        roi(id(1),id(2),id(3))=j;
        roi_path(id(1),id(2),id(3))=j;       
     %   m_background(id(1),id(2),id(3))=2;     
    end
 
end

roi=single(roi>0);
%roi_path=single(roi_path);
nii.img=roi;
prefix=strtok2(filename(ftrace),'.');
save_untouch_niigz(nii,[prefix,'_sig',num2str(sig_scale)]);

%nii.img=roi;
%prefix=strtok(ftrace,'.');
%save_untouch_niigz(nii,prefix);

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

function pos2=connect_pos(pos)

pos2=[];
for i=1:size(pos,1)-1
    
    n=pos(i+1,:)-pos(i,:);
    
    npix=ceil(sqrt(sum(n.^2)));
    
    for j=0:npix-1
        tmp=round(pos(i,:)+j*n/npix);
        if isempty(pos2) || any(tmp~=pos2(end,:))
            pos2(end+1,:) = tmp;
        end
    end
    
    
end

if  any(pos(end,:)~=pos2(end,:))
    pos2(end+1,:) = pos(end,:);
end



function check_dir(pos)

for i=1:size(pos,1)-2
    
    n=pos(i+1,:)-pos(i,:);
    n2=pos(i+2,:)-pos(i+1,:);
    
    c=sum(n.*n2)/sos(n)/sos(n2);
    %disp([i,c]);
    if c<0
        disp([i,c]);
        disp('path error');
    end
    
    
end



