function m2=detectPVSCircleInMask(pc,m_file,r,thr_scale,neg_phase,mag_file,thr_scale_mag)
% m2=detectPVSCircleInMask(pc,m_file,r,thr_scale,neg_phase[,mag_file,thr_scale_mag])
% 2/21/2017: added the option of magnitude thresholding
%
if ~exist('thr_scale','var') || isempty(thr_scale)
    thr_scale=3;
end

if ~exist('r','var') || isempty(r)
    r=10;
end

if ~exist('thr_scale','var') || isempty(thr_scale)
    thr_scale=3;
end


if ~exist('thr_scale_mag','var') || isempty(thr_scale_mag)
    thr_scale_mag=1.3;
end


a=ri(pc);

if exist('mag_file','var')
mag=ri(mag_file);
end

prefix=strtok(pc,'.');

try
    m=ri(m_file,'','','d');
catch
m=ri(m_file);
end
[xx,yy]=meshgrid(1:size(m,2),1:size(m,1));


if exist('neg_phase','var') && ~isempty(neg_phase) && neg_phase
    a=-a;
else
    neg_phase=false;
end

m2=m*0;

for k=1:size(m,3)
    
        atmp=a(:,:,k);
 
    for i=1:size(m,1)
        fprintf('%d ',i);
        if mod(i,20)==0
            fprintf('\n');
        end
        for j=1:size(m,2)
            
            if m(i,j,k)==0 || i<r+1 || j<r+1 || j>size(m,2)-r || i>size(m,1)-r
                continue;
            end
            
            
            mloc=(xx-j).^2+(yy-i).^2<=r*r & m(:,:,k)>0;
       %     mloc=0*m(:,:,k);
       %     mloc(i-r:i+r,j-r:j+r)=1;
       %     mloc=mloc>0&m(:,:,k)>0;
            
            x=atmp(mloc>0);
            mn=mean(x);
            sd=std(x);
            
            if exist('mag','var')
                mag_tmp=mag(:,:,k,1,1);
                x=mag_tmp(mloc>0);
                mn_mag=mean(x);
                sd_mag=std(x);
                if atmp(i,j)>mn+thr_scale*sd && mag_tmp(i,j)>mn_mag+thr_scale_mag*sd_mag
                    m2(i,j,k)=1;
                end
            else
                if atmp(i,j)>mn+thr_scale*sd 
                    m2(i,j,k)=1;
                end      
            end
            
        end
    end
    
    tmp=clusterize2(m2(:,:,k));
    fprintf('\n%d clusters detected in slice %d\n',max(tmp(:)),k);
    
end

if neg_phase
neg_str='_negPhase';
else
    
neg_str='_posPhase';
end

d=m2;
[dir_name,fname]=fileparts(pc);
if isempty(dir_name)
    dir_name='.';
end
prefix=strtok(fname,'.');
cur_dir=cd(dir_name);

try
    voxsize=ri(m_file,[],[],'vxosize');
    
    center=ri(m_file,[],[],'center');
    
    save([prefix,'_Vssl.mat'],'d','voxsize','center');
    
catch
    save([prefix,'_Vssl.mat'],'d');
      
end
cd(cur_dir);