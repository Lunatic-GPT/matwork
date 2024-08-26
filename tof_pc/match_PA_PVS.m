function match_PA_PVS(fPA,fDir,fPC,fT2,sep_max,oname)

%Inputs:
%fPA: nifti file for PA mask
%fDir: nifti file for resampled PVS dirs
%fPC: PC nifti file, needed for converting dir from image index space to physical space.
%fT2: 3D T2 mask, needed for converting dir from image index space to physical space.
%sep_max: max separaton in mm
%oname: save results as a .mat file which contains
% PAmatch: PA mask with clusters that can match PVS labelled as 1 .. nPA_match,
% while unmatched ones labelled as 1001,1002,....
% PVSmatch: PVS mask:  same as values as PA mask, Same index is used for
% matched PA-PVS.
% matched with 
% angle_match [nPA_match*2]: Angles of the matched PVS relative to the
% slice normal direction in degrees, the index matches with the PAmatch.
% first and second columns are calculated from dir_whole and dir_seg, respectively.
% angle_nomatch: [nPA_nomatch*2].the index matches with the PAmatch after
% subtracting 1000.
        
    

    nii_T2=load_untouch_niigz(fT2);
    nii_pa=load_untouch_niigz(fPC);
    m_pa=ri(fPA);
    pa=clusterize2(m_pa>0);
    dr=ri(fDir);
    pvs=clusterize2(dr(:,:,1,1)>0); %otherwise, different PVSs might be connected.
    
    pixdim=nii_pa.hdr.dime.pixdim;
    npix=round(sep_max/pixdim(2));  %maximum separation;
    
    [ipvs_ipa,~,ipvs_ipa2]=clusters_match(pvs,pa,npix,false);
    
    ipvs_ipa=cat(1,ipvs_ipa,ipvs_ipa2);
    PAmatch=get_mask_match(pa,ipvs_ipa(:,2));
    PVSmatch=get_mask_match(pvs,ipvs_ipa(:,1));
    
    [angle_match(:,1),angle_nomatch(:,1)]=get_angles(PVSmatch,dr(:,:,1,2:4),nii_T2,nii_pa);
    [angle_match(:,2),angle_nomatch(:,2)]=get_angles(PVSmatch,dr(:,:,1,5:7),nii_T2,nii_pa);
    npa_total=max(pa(:));
    
    npvs_total= max(pvs(:));

    
    npa_match=size(ipvs_ipa,1);
    
    save(oname,'ipvs_ipa','angle_match','angle_nomatch','PAmatch','PVSmatch','npa_total','npvs_total','npa_match');
    
    
    
    function  [angle_match,angle_nomatch]=get_angles(m,dr,nii_T2,nii_pa)
        
        npa_match=max(m(m<1000));
        npa_nomatch=max(m(:))-1000;
        angle_match=NaN(1,npa_match);
        angle_nomatch=NaN(1,npa_nomatch);
        
        orient_pa=get_orient_from_nii(nii_pa);
        orient_T2=get_orient_from_nii(nii_T2);
        
        for i=1:npa_match
            dr_mean=mean_roi(dr,m==i);
            x=orient_pa.rotmat(:,3);
            y=orient_T2.rotmat*dr_mean(:);
            
            if any(abs(y)>0)
              angle_match(i)=acos(sum(x.*y)/sos(x,1)/sos(y,1))*180/pi;
            end
            
        end
       for i=1:npa_nomatch
            dr_mean=mean_roi(dr,m==i+1000);
            x=orient_pa.rotmat(:,3);
            y=orient_T2.rotmat*dr_mean(:);
            
            if any(abs(y)>0)
              angle_nomatch(i)=acos(sum(x.*y)/sos(x,1)/sos(y,1))*180/pi;
            end
            
        end
    
    
    function res=get_mask_match(pa,index_match)
        
        res=zeros(size(pa));
        c_nomatch=0;
        for i=1:max(pa(:))
            id=find(index_match==i);
            
            if isempty(id)
                if any(pa(:)==i)
                 c_nomatch=c_nomatch+1;
                 res(pa==i)=c_nomatch+1000;
                end
            else
                res(pa==i)=id;
            end
            
        end
    
    
    
    
    