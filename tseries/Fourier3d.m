function Fourier3d(dname,lpf,hpf,out_name)
%Fourier3d(dname,lpf,hpf,out_name)
% use [] for lpf and hpf if you don't want to do lowpass or high pass filtering. 
[data,info]=BrikLoad(dname);
TR=info.TAXIS_FLOATS(2);

npt = size(data,4);
fd= fft(data,[],4);

if ~isempty(lpf)
  lw = ceil(lpf*(npt*TR))+1;
  up = npt-(lw-2);
  fd(:,:,:,lw:up)=0;
end

if ~isempty(hpf)
  hpf_ctof = round(hpf*npt*TR)+1;
  fd(:,:,:,1:hpf_ctof) = 0;
  if hpf_ctof > 1
    fd(:,:,:,npt-hpf_ctof+2:npt) = 0;
  end
  data = ifft(fd,[],4);
end

info.DATASET_RANK(2) = npt;  % the number of subbricks to save.
info.BRICK_TYPES= 3*ones(1,npt);  %1: short, 3: float
      
info.IDCODE_STRING = 'low passed filtered time course';
       
       info.TYPESTRING = '3DIM_HEAD_ANAT';
       info.SCENE_DATA = [0 2 0];  % first 0: orig view; second 2: ANAT_EPI_TYPE; 0: matches TYPESTRING 
       
       info.BRICK_STATS=[];  %minimum and maximum values of the subbrick;
       
       info.BRICK_FLOAT_FACS = [];
       
       info.IDCODE_DATE = date;
       info.TAXIS_NUMS=[];
       
info.BRICK_LABS = [];
opt.Prefix = out_name;
opt.OverWrite='y';
history=sprintf('Fourier3d(dname=%s,lpf=%3.2f,hpf=%3.2f,out_name=%s)',dname,lpf,hpf,out_name);
info.HISTORY_NOTE = [info.HISTORY_NOTE,'\n',history];
WriteBrik(data,info,opt);




