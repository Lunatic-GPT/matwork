function Eve_4ROIs(m_eve,m_wm)

%1. BG
%2. MB
%3. TH
%4. WM

save_log('Eve_4ROIs',m_eve,m_wm);
dir_str=dir(m_eve);
thalamus=[149];  %
BG = [147,148,150];  %caudate, putaman, globus pallidus
midbrain=[144,151];

corona_radiata=[125,126,127];
ic=[122,123];
ec=[135];
cc=[140,141,142];

occip=[163,164,166,167,168];
parietal=[154,160,161,155,162,174];
temporal=[165,169,170,171];
frontal=[156,157,158,159,172,173,175];

WM=[corona_radiata,ic,ec,cc,occip,parietal,temporal,frontal];

    
    nii=load_untouch_niigz(fullfile(dir_str.folder,dir_str.name));
    wm=ri(m_wm);
    m_wm=wm>0.9;
    m=0*nii.img;
    m_eve=nii.img;
    m_th=m_eve==149|m_eve==149-88;
    m_bg=m_eve==147|m_eve==147-88|m_eve==148|m_eve==148-88|m_eve==150|m_eve==150-88;
    m_mb=m_eve==144|m_eve==144-88|m_eve==151|m_eve==151-88;
    m(m_th)=3;
    m(m_bg)=1;
    m(m_mb)=2;
   % m_wm=roi_from_number(WM,m_eve);
    m(m_wm>0)=4;
    m(m_wm>0&(m_th>0|m_bg|m_mb))=0;
    nii.img=m;
    prefix=strtok(dir_str.name,'.');
    out=fullfile(dir_str.folder,[prefix,'_4ROI']);
    save_untouch_niigz(nii,out);


function res=roi_from_number(n,roi)

res=0*roi;
for i=1:length(n)
    res=res|roi==n(i)|roi==n(i)-88;
end


