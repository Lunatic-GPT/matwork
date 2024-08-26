function FatNav_Ice()

global m_iAccFactor m_iNCh m_dFFTFactor;
global m_DataAs m_DataRefAs;
global m_iRawCol m_iTotalPar m_iNRefLin m_iNRefPar;
global m_iTotalLin m_iFirstRefLin m_iFirstRefPar;

m_iAccFactor=4;
m_iFirstRefLin = 11;
m_iFirstRefPar = 21;

EndInit ();

m_iNCh=32;
m_iRawCol=64;
m_iTotalPar=64;
m_iNRefLin=22;
m_iNRefPar=22;
m_iTotalLin=44;
m_DataAs=zeros(m_iRawCol, m_iNCh,m_iTotalLin,m_iTotalPar);
m_DataRefAs=zeros(m_iRawCol, m_iNCh,m_iNRefLin,m_iNRefPar);
m_dFFTFactor=100;


load meas_MID17_FN_44lines_centerLinParset_FID44154_FatNav.mat
Data=reshape(Data,[64,32,79200/32]);
Line=Line(1:32:end);
Partition=Partition(1:32:end);
Repetition=Repetition(1:32:end);
Flags=Flags(1:32:end);

for i=1:size(Data,3)
    pscan.lin=Line(i);
    pscan.par=Partition(i);
    pscan.rep=Repetition(i);
    pscan.bImageScan=isImagingScan(Flags(i));
    pscan.bPatRefScan=determineBitFields(Flags(i),'MDH_PATREFSCAN');
    pscan.bPatRefAndImageScan=determineBitFields(Flags(i),'MDH_PATREFANDIMASCAN');
    pscan.bLastScan=determineBitFields(Flags(i),'MDH_PHASEFFT')&determineBitFields(Flags(i),'MDH_D3FFT');
    ComputeScan(Data(:,:,i),pscan);
end

function EndInit ()

global m_lx m_ly m_kernelX m_kernelY nkernel m_iAccFactor;

m_lx = 2;
m_ly = 2;
m_kernelX=[-1,3;-2,2;-3,1;-4,0]';
m_kernelY=m_kernelX;
nkernel=m_iAccFactor*m_iAccFactor-1;

function ComputeScan(data,pscan)

global m_DataRefAs m_iFirstRefLin m_iFirstRefPar m_iNRefPar m_iNCh m_iRawCol nkernel m_lx m_ly m_iNRefLin;
global m_iAccFactor
global m_QA m_b m_DataAs m_coef

lin=pscan.lin;
par=pscan.par;
rep=pscan.rep;
bImageScan=pscan.bImageScan;
bPatRefAndImageScan=pscan.bPatRefAndImageScan;

if pscan.bPatRefScan%(~bImageScan||bPatRefAndImageScan)
    
    
    m_DataRefAs(:,:,lin-m_iFirstRefLin+1,par-m_iFirstRefPar+1)=data;
    
    
    if (par-m_iFirstRefPar == m_iNRefPar-1 && lin-m_iFirstRefLin == m_iNRefLin-1)
        
        m_QA=zeros(m_iNRefLin*m_iNRefPar,m_lx*m_ly*m_iNCh,m_iRawCol,nkernel);
        m_b=zeros(m_iNRefLin*m_iNRefPar, m_iNCh,m_iRawCol);
        m_coef=zeros(m_lx*m_ly*m_iNCh,m_iNCh,m_iRawCol,nkernel);
        for ikx=0:m_iAccFactor-1
            for iky=0:m_iAccFactor-1
                
                if (ikx==m_iAccFactor-1 && iky==m_iAccFactor-1)
                    break;
                end
                fillCalibrationMatrix(ikx,iky);
                
            end
        end
        
        %m_b dimension: m_iNRefLin*m_iNRefPar, m_im_iNCh,m_iRawCol
        % m_QA dimension: [m_iNRefLin*m_iNRefPar,m_lx*m_ly*m_im_iNCh,m_iRawCol,nkernel]);
        %  m_coef dimensions: [m_lx*m_ly*m_iNCh,m_iNCh,m_iRawCol,nkernel]);
        
        for ikx=0:m_iAccFactor-1
            for iky=0:m_iAccFactor-1
                
                if (ikx==m_iAccFactor-1 && iky==m_iAccFactor-1)
                    break;
                end
                
                ikernel=ikx+iky*m_iAccFactor;
                
                for iro=0:m_iRawCol-1
                    
                    m_QAtmp=m_QA(:,:,iro+1,ikernel+1);
                    m_QAtQA=m_QAtmp'*m_QAtmp;
                    
                    m_btmp=m_b(:,:,iro+1);
                    
                    m_QAtb=m_QAtmp'*m_btmp;
                    
                    m_coef(:,:,iro+1,ikernel+1)=m_QAtQA\m_QAtb;
                    
                end
            end
        end
    end
end
if (bImageScan)
    if rep==1
        disp('');
    end
    m_DataAs(:,:,lin+1,par+1)=data;
    if (pscan.bLastScan)
        
%%
load('../coef/coef.mat', 'coef')
figure;plot(real(coef(:)),real(m_coef(:)),'.');
load('../k_data_rep1/kdata.mat', 'k')
figure;plot(real(k(:)),real(vec(m_DataAs(:,4,:,:))),'.');
        
%%
do_Grappa();
sendImages(rep);
m_DataAs=0*m_DataAs;
    end
end


function sendImages(rep)

global m_DataAs m_dFFTFactor

m_ImageAs_ch=abs(ifft1c(ifft1c(m_DataAs,3),4))*m_dFFTFactor;

m_ImageAs_ch=m_ImageAs_ch.^2;

m_ImageAs=squeeze(sqrt(sum(m_ImageAs_ch,2)));


save(['m_ImageAs_rep',num2str(rep)], 'm_ImageAs');


function do_Grappa()
global m_coef m_lx m_ly m_iNCh m_iRawCol m_iTotalLin m_iTotalPar m_iAccFactor m_kernelX m_kernelY m_DataAs 


 %  m_coef dimensions: [m_lx*m_ly*m_iNCh,m_iNCh,m_iRawCol,nkernel]);

lin0=m_iTotalLin/2;
par0=m_iTotalPar/2;
for lin=0:m_iTotalLin-1
    disp(lin);
    for par=0:m_iTotalPar-1
        
        ikx=mod(lin-lin0-1,m_iAccFactor);
        iky=mod(par-par0-1,m_iAccFactor);
        
        if (ikx<0) ikx=ikx+m_iAccFactor; end
        if (iky<0) iky=iky+m_iAccFactor; end
        
        if (ikx==m_iAccFactor-1 && iky==m_iAccFactor-1) continue; end
        
         if abs(m_DataAs(m_iRawCol/2,1,lin+1,par+1))>0  %new
             warning('This should not happen');
         end
         
        ikernel=ikx+iky*m_iAccFactor;
        
        for iro=0:m_iRawCol-1
            m_grappaArray=zeros(1,m_iNCh);
              
            for i=0:m_lx-1
                for j=0:m_ly-1
                    
                    iy = par + m_kernelY(iky*m_ly+j+1);
                    ix = lin+ m_kernelX(ikx*m_lx+i+1);
                    if (ix<0 || ix>=m_iTotalLin) continue; end
                    if (iy<0 || iy>=m_iTotalPar) continue; end
                    
                  
                    
                    for ich=0:m_iNCh-1                   
                    %    m_value=m_DataAs(iro+1,ich+1,ix+1,iy+1);
%                         for ich_c=0:m_iNCh-1
%                            % ind=ich_c+m_iNCh*(i+j*m_lx+ich*m_lx*m_ly+iro*m_lx*m_ly*m_iNCh);
%                            % m_grappaArray(ich_c+1)=m_grappaArray(ich_c+1)+m_value*m_coef(ind+1,ikernel+1); 
%                         end
                        ind=i+j*m_lx+ich*m_lx*m_ly;
                        m_grappaArray=m_grappaArray+m_DataAs(iro+1,ich+1,ix+1,iy+1)*m_coef(ind+1,:,iro+1,ikernel+1);
                    end
                    
                   
                    
                    
                end
            end
             m_DataAs(iro+1,:,lin+1,par+1)=m_grappaArray;
        end
    end
end

function fillCalibrationMatrix(ikx,iky)
global m_iNRefLin
global m_iNRefPar
for lin=0:m_iNRefLin-1
    for par=0:m_iNRefPar-1
        
        
        FillQA(lin,par,ikx,iky);
        Fill_b(lin,par,ikx,iky);
    end
end


function Fill_b(iLin, iPar,ikx,iky)

global m_iNCh m_iRawCol m_iNRefLin m_iNRefPar m_b m_DataRefAs
OutOfBound = IsOutOfBound(iLin, iPar,ikx,iky);

%m_b dimension: m_iNRefLin*m_iNRefPar, m_iNCh,m_iRawCol
% m_QA dimension: [m_iNRefLin*m_iNRefPar,m_lx*m_ly*m_im_iNCh,m_iRawCol,nkernel]);
%  m_coef dimensions: [m_lx*m_ly*m_iNCh,m_iNCh,m_iRawCol,nkernel]);

for ich=0:m_iNCh-1
    for iro=0:m_iRawCol-1
        
        row=iLin+iPar*m_iNRefLin;
        offset=iro*m_iNCh*m_iNRefLin*m_iNRefPar;
        
        if ~OutOfBound
            m_value=m_DataRefAs(iro+1,ich+1,iLin+1,iPar+1);
            
            
            ind=row+ich*m_iNRefLin*m_iNRefPar+offset;  % row as the fast index to speep up calculation
            
            m_b(ind+1)=m_value;
            
        end
    end
end


function OutOfBound=IsOutOfBound(iLin, iPar,ikx,iky)

global m_lx m_ly m_kernelX m_kernelY m_iNRefLin m_iNRefPar
OutOfBound = false;

for i=0:m_lx-1
    for j=0:m_ly-1
        
        iy = (iPar + m_kernelY(iky*m_ly+j+1));
        ix = (iLin+ m_kernelX(ikx*m_lx+i+1));
        if ( iy<0 || iy>=m_iNRefPar || ix<0 || ix>=m_iNRefLin)
            
            OutOfBound = true;
            break;
        end
    end
end


function FillQA(iLin, iPar,ikx,iky)

global m_QA m_iNRefLin m_iNRefPar m_iRawCol m_lx m_ly m_iNCh nkernel
global m_kernelX m_kernelY m_DataRefAs m_iAccFactor
m_QA=reshape(m_QA,[m_iNRefLin*m_iNRefPar*m_iRawCol*m_lx*m_ly*m_iNCh,nkernel]);


OutOfBound = IsOutOfBound(iLin, iPar,ikx,iky);

for  i=0:m_lx-1
    for j=0:m_ly-1
        
        if (~OutOfBound)
            
            iy = iPar + m_kernelY(iky*m_ly+j+1);
            ix = iLin+ m_kernelX(ikx*m_lx+i+1);
            
            
            for ich=0:m_iNCh-1
                
                for iro=0:m_iRawCol-1
                    
                    
                    m_value=m_DataRefAs(iro+1,ich+1,ix+1,iy+1);
                    
                   % rowt=i+j*m_ly+ich*m_lx*m_ly+iro*(m_lx*m_ly*m_iNCh);
                   % colt=iLin+iPar*m_iNRefLin;
                   % indt=colt+rowt*m_iNRefLin*m_iNRefPar;
                    
                    row=iLin+iPar*m_iNRefLin;
                    col=i+j*m_lx+ich*m_lx*m_ly+iro*(m_lx*m_ly*m_iNCh);
                    ind = row+col*m_iNRefLin*m_iNRefPar;
                    ikernel=ikx+iky*m_iAccFactor;
                    
                    m_QA(ind+1,ikernel+1)=m_value;
                end
            end
        end
    end
end

m_QA=reshape(m_QA,[m_iNRefLin*m_iNRefPar,m_lx*m_ly*m_iNCh,m_iRawCol,nkernel]);




