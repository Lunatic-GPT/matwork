function epiref(ref)
% epiShaper32(fid_dir,triple_ref,format,nophase)
% nophase: do not save phase image.default true.
% format: 'a' or 's', default 'a';



sz=readbPar(fullfile(ref,'method'),'PVM_EncMatrix');
ref=readb_fid(ref);
ref=reshape(ref,[sz(1),sz(2)]);



ref=convertTraces(ref);
figure;imshow(abs(ref(:,:)),[]);

function z = convertTraces(z_tmp)
 % reverse the direction of even number traces
 z = zeros(size(z_tmp));
 
 z(:,1:2:end-1,:,:,:) = z_tmp(:,1:2:end-1,:,:,:);
 z(1:end,2:2:end,:,:,:) = z_tmp(end:-1:1,2:2:end,:,:,:);

        