function [roi_str,roi_name]=mask_str(mask)
%[roi_str,roi_name]=mask_str(mask)
% generate mask strings for rois generated with sroi2v function. 
[m,info] = BrikLoad(mask);
label = info.BRICK_LABS;

if label(end) ~='~'
    label(end+1) = '~';
end
tilde = find(label=='~');

ntd = length(tilde);

if (mod(ntd-3,3)~=0) 
    error('Wrong roi label format');
end
nrois = (ntd-3)/3;

fprintf('%d rois found\n',nrois);
roi_str = cell(1,nrois);
roi_name = cell(1,nrois);

fprintf('                      ');
for i=1:nrois

    len = 0;
    while label(tilde(3+i)-len-1)~=' ' ...
       && label(tilde(3+i)-len-1)~='/' ...
       && label(tilde(3+i)-len-1)~='\'
      len=len+1;
    end
    
    roi_name{i}= label(tilde(3+i)-len:tilde(3+i)-1);
    roi_str{i} = sprintf('%s[%d]',mask,ntd-nrois+i-1);
    
    fprintf('%7s',roi_name{i});
    
    
end

    fprintf('\n');
    
    for j=1:nrois-1
        nmap = length(find(m(:,:,:,j+3)>0));
        nuniq = length(find(m(:,:,:,j+3)>0 & m(:,:,:,3) ==1));
        nmax = length(find(m(:,:,:,j+3+2*(nrois))>0 )); 
       fprintf('%6s %4d (%3.2f,%3.2f) ',roi_name{j},nmap,nmax/nmap,nuniq/nmap);
        for i=1:nrois
            nol = length(find(m(:,:,:,j+3)>0 & m(:,:,:,i+3)>0));
            fprintf('%7.2f',nol/nmap);
        end
        fprintf('\n');
        
    end
    
    