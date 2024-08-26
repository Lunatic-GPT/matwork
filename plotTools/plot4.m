function plot4(d,shape,rg)

for i=1:size(d,2)
   
    subplot(shape(1),shape(2),i);
    
    plot(d(:,i));
    if exist('rg','var')
    ylim(rg);
    end
end