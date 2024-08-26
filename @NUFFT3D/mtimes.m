function ress = mtimes(a,bb)
% performs the normal nufft

if a.adjoint
    nd=2;
    ress = zeros([a.imSize,size(bb,nd)],'single');
    
else
    nd=4;
    ress = zeros([a.dataSize(1),size(bb,nd)],'single');
    
end
% 
% if a.pools>1
%     parpool(a.pools);
% end


imSize=a.imSize;
w=a.w(:);
dataSize = a.dataSize;

ind=linspace(0,dataSize(1),a.nblock+1);
ind=round(ind);

Nc = size(bb,nd);
for i=1:a.nblock
    tic;
    st = nufft_init(a.om(ind(i)+1:ind(i+1),:), a.Nd, a.Jd, a.Kd, a.n_shift,'kaiser');
    tinit =toc;
    fprintf('nufft_init time = %d\n',round(tinit));
    if a.adjoint
        if length(w(:))>1
             bb=bb.*repmat(w,[1,size(bb,nd)]);
        else
             bb=bb.*w;
        end
            for m=1:Nc
                tic;
                b = bb(ind(i)+1:ind(i+1),m);
                % b = b(:).*w;
                res = nufft_adj(double(b), st)/sqrt(prod(imSize));
                ress(:,:,:,m)= ress(:,:,:,m)+reshape(res, imSize);
                fprintf('estimated total = %d s; remaining = %d s; nufft_adj per chan = %f s\n',...
                    round(toc*Nc*a.nblock+tinit*a.nblock),round(toc*(Nc-m+(a.nblock-i)*Nc)+tinit*(a.nblock-i)),toc);
            end
    else
            for m=1:Nc
                tic;
                b = bb(:,:,:,m);
                b = reshape(b,imSize);
                res = nufft(double(b), st)/sqrt(prod(imSize));
                ress(ind(i)+1:ind(i+1),m) = reshape(res,[ind(i+1)-ind(i),1]);
                
                fprintf('estimated total = %d s; remaining = %d s; nufft per chan = %f s\n',...
                    round(toc*Nc*a.nblock+tinit*a.nblock),round(toc*(Nc-m+(a.nblock-i)*Nc)+tinit*(a.nblock-i)),toc);
            end
%         end
        if length(w(:))>1
             ress=ress.*repmat(w,[1,size(ress,2)]);
        else
            
            ress=ress*w;
        end
    end
end
