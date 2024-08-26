function [v2,sl_u2o]=reslice_extend(odata,o_voxsize,u_voxsize,udsz,center_u2o,use_nn)
 % udsz: matrix sizes (1st, 2nd, and 3rd dim) of underlay    
% center_u2o: the difference of center positions between underlay and
% overlay, i.e. center_u - center_o, assuming z increase with slice index.


if length(udsz)==2
    udsz(3)=1;
end


sz_v2=udsz;


%sz_v2(3)=min([udsz(3),ceil(o_voxsize(3)*size(odata,3)/u_voxsize(3))]);  %

if udsz(3)>ceil(o_voxsize(3)*size(odata,3)/u_voxsize(3))
    % assume overlay within underlay
    if size(odata,3)==1
        sz_v2(3)=1;
    end
    nvox2shft = (center_u2o./u_voxsize);
    center_u2o(3)=(round(nvox2shft(3))-nvox2shft(3)).*u_voxsize(3);
    sl_u2o=round(nvox2shft(3));  % difference of the center slice
    sl_u2o=sl_u2o+ceil((sz_v2(3)+1)/2)-ceil((udsz(3)+1)/2);
    
    if abs(u_voxsize(3)-o_voxsize(3))<0.1
        if  abs(floor(nvox2shft(3))-nvox2shft(3))<0.1 || abs(ceil(nvox2shft(3))-nvox2shft(3))<0.1
            twoD = true;
        else
            twoD = false;
        end
    else
        twoD=false;
    end
    
    if ~exist('use_nn','var')
        if ~any(u_voxsize~=o_voxsize)
            
            if ~any( abs(round(nvox2shft)-nvox2shft)>0.1)
                use_nn=true;
                
            else
                use_nn=false;
            end
            
        else
            use_nn=false;
            
        end
    end
    
else
    % assume underlay within overlay
    sl_u2o=0;
    twoD = false;
    if ~exist('use_nn','var')
      use_nn = true;
    end
end


v2=reslice(odata,eye(3),o_voxsize,center_u2o,sz_v2,u_voxsize,twoD,use_nn);

%v2(isnan(v2))=0;

%need to calculate twoD, sz_v2, and sl_u2o;

