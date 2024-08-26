function csi_voxels_afni(fname,trnc,phase_adjust,baseline,prefix)
% csi_voxels_afni(fname,prefix,trnc,phase_adjust)
% trnc: 1*2; truncate the time series at the beginning and end of the fid
% phase_adjust: global phase adjustment in degree


if ~exist('prefix','var') || isempty(prefix)
    prefix=fname;
end

prefix_orig=prefix;
prefix=sprintf('%s_%d_%d_%d',prefix,trnc(1),trnc(2),phase_adjust);

[d,in]=dicomreads(fname);
% in.RealDwellTime

d=reshape(d,[cell2num(in.DataPointColumns(1)),...
             cell2num(in.Rows(1)),cell2num(in.Columns(1)),cell2num(in.NumberOfFrames(1))]);
d=permute(d,[2,3,4,1]);



save(prefix,'d');

d=d(:,:,:,trnc(1)+1:end-trnc(2));

d=d-repmat(mean(d(:,:,:,end-baseline+1:end),4),[1,1,1,size(d,4)]);

sz=size(d);

d=d*exp(1i*phase_adjust/180*pi);


fd=fft(d,[],4);
fd=fftshift(fd,4);


 info.DATASET_DIMENSIONS = size(d(:,:,:,1));
 info.DATASET_RANK=  [3, sz(4)];
 info.BRICK_TYPES=ones(1,sz(4))*3;
 
 d2=reshape(d,[prod(sz(1:3)),sz(4)]);
 min_val=min(d2,[],1);
 max_val=max(d2,[],1);
 stat=[min_val;max_val];
 info.BRICK_STATS=stat(:)';
 
 info.BRICK_FLOAT_FACS=zeros(1,sz(4));
 info.BYTEORDER_STRING= 'LSB_FIRST';
 
     orient(1)=vector_2_afni_orient_code(cell2num(in.ImageOrientationPatient(1:3)));
     orient(2)=vector_2_afni_orient_code(cell2num(in.ImageOrientationPatient(4:6)));
     orient(3)=vector_2_afni_orient_code(cell2num(in.VoiOrientation(1:3)));
     info.ORIENT_SPECIFIC = orient;
 
 
    info.ORIGIN = cell2num(in.ImagePositionPatient(1:3));
    info.DELTA = [cell2num(in.PixelSpacing(1:2)),cell2num(in.SliceThickness(1))];
        
    for i=1:3
        if any(orient(i)==[1,2,5])
         info.DELTA(i)=-info.DELTA(i);
        end
    end
    
    
    
    info.BRICK_LABS='';
    info.BRICK_KEYWORDS='';
    
    info.SCENE_DATA =[0 2 0];
    info.TYPESTRING= '3DIM_HEAD_ANAT';
    
    info.TAXIS_NUMS = [sz(4) 0 77002];
    
    info.TAXIS_FLOATS= [0 5 0 -18.1426 6];
    
    info.TAXIS_OFFSETS=0; % not used;
    opt.OverWrite='y';
    opt.Prefix=sprintf('%s_Real',prefix);
    WriteBrik(real(fd),info,opt);
    
    opt.Prefix=sprintf('%s_Imag',prefix);
    WriteBrik(imag(fd),info,opt);
    
    opt.Prefix=sprintf('%s_Mag',prefix);
    WriteBrik(abs(fd),info,opt);
    
    
%{
               IDCODE_STRING: 'XYZ_0j3bEsuOhT5fH4tA2zpHeA'
                 IDCODE_DATE: 'Thu Jan 28 13:13:51 2016'
               BRICK_STATAUX: []
                    STAT_AUX: []
                HISTORY_NOTE: '[xiaopeng@ubuntu: Thu Jan 28 13:13:51 2016] ...'
                 NOTES_COUNT: []
             NOTE_NUMBER_001: []
             TAGALIGN_MATVEC: []
           VOLREG_CENTER_OLD: []
          VOLREG_CENTER_BASE: []
     VOLREG_ROTPARENT_IDCODE: []
       VOLREG_ROTPARENT_NAME: []
    VOLREG_GRIDPARENT_IDCODE: []
      VOLREG_GRIDPARENT_NAME: []
         VOLREG_INPUT_IDCODE: []
           VOLREG_INPUT_NAME: []
          VOLREG_BASE_IDCODE: []
            VOLREG_BASE_NAME: []
           VOLREG_ROTCOM_NUM: []
           
           IJK_TO_DICOM_REAL: [1x12 double]
           
          IDCODE_ANAT_PARENT: []
                   TO3D_ZPAD: []
          IDCODE_WARP_PARENT: []
                   WARP_TYPE: []
                   WARP_DATA: []
                   MARKS_XYZ: []
                   MARKS_LAB: []
                  MARKS_HELP: []
                 MARKS_FLAGS: []
                  TAGSET_NUM: []
               TAGSET_FLOATS: []
               TAGSET_LABELS: []
                     LABEL_1: 'pcasl2+orig'
                     LABEL_2: 'Viggo!'
                DATASET_NAME: './pcasl2+orig'
            DATASET_KEYWORDS: []
                  WORSLEY_DF: []
               WORSLEY_NCONJ: []
                WORSLEY_FWHM: []
        }
                    TypeName: 'short'
                   TypeBytes: 2
                   ByteOrder: 'ieee-le'
                 Orientation: [3x2 char]
                  FileFormat: 'BRIK'
         %       Extension_1D: ''
        %}        
                
    function orient=vector_2_afni_orient_code(vect)
                
                [tmp,max_ind2]=max(abs(vect));
                
                if max_ind2==1 
                    if vect(max_ind2)>0 
                        orient=0;
                    else
                        orient=1;
                    end
                             
                elseif max_ind2==2
                    
                     if vect(max_ind2)>0 
                        orient=3;
                    else
                        orient=2;
                     end
                else
                    
                    if vect(max_ind2)>0 
                        orient=4;
                    else
                        orient=5;
                    end
                end
                
                
                
function res=cell2num(d)
                        
                        res=zeros(1,length(d));
                        for i=1:length(d)
                            
                            res(i)=str2num(d{i});
                        end
                        
                            
 
                
                