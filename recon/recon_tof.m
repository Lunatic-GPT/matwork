function recon_tof(filename,interp_factor)
% Example:
%
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70]);
% kSpaceData = readMeasDat('meas_yyT2.dat',[204 256 70],'Echo',2,'Rep',50);
%filename='meas_MID114_TOFZ2D_0_3cm_s_FID11616.dat';
%

%interp_factor=10;
%filename='meas_MID106_TOFZ2D_noflow_FID11608.dat';
[d,lin,par,sl]=readMeasDat(filename,inf,0,true);

prefix=strtok(filename,'.');

prefix_dir=fullfile(prefix,prefix);
nsl=readsPar([prefix_dir,'.pro'],'lConc');
nave=readsPar([prefix_dir,'.pro'],'lAverages');


lin2=lin(1:32:end);
lin2=lin2(3:end/nsl);
lin2=reshape(lin2,[length(lin2)/nave,nave]);
lin2=lin2(:,1);


nro=size(d,1);
d=reshape(d,[nro,32,length(d(:))/nro/32]);
dref=d(:,:,1:2);

ps=readsPar([prefix_dir,'.pro'],'ucPhaseStabilize');

if strcmp(ps,'0x01')
    d2=reshape(d(:,:,3:end),[nro,32,2,(size(d,3)-2)/nsl/nave/2,nave,nsl]);
else
    d2=reshape(d(:,:,3:end),[nro,32,(size(d,3)-2)/nsl/nave,nave,nsl]);
end


%%
if nsl>1
ucMode = readsPar([prefix_dir,'.pro'],'sSliceArray.ucMode');

if strcmp(ucMode,'0x4')  % interleaved
    d2=interleave2linear(d2,5);  % the resulting first slice on the feet side.
elseif strcmp(ucMode,'0x2')
    d2=flip(d2,5);
end
end

d3=mean(d2,4);

d3=permute(d3(:,:,1:2:end,1,:),[1,3,2,5,4]);

%%

%PixelSizeFactor = 1.275;

fd3=fft1c(d3,1);
sz0=size(fd3);
sz=sz0;
sz(1:2)=round(sz0(1:2)*interp_factor);

if length(sz)==3  
    sz(4)=1;
end

fd3=zpad(fd3,sz(1),sz(2),sz(3),sz(4));
im=squeeze(sos(ifft2c(fd3),3));
%cut=round((sz(1:2)-sz0(1:2)*interp_factor)/2);

%im=im(cut+1:end-cut,cut+1:end-cut,:,:);
save(sprintf('recon_%s_interp%d.mat',prefix,interp_factor),'im');




