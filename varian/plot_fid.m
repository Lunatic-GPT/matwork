function plot_fid(a)
% plot_fid(a,ind)

uicontrol(gcf,'style','slider');

figure;

  subplot(2,1,1);
  plot(real(a),'r');
  hold on; plot(imag(a),'b');
  plot(abs(a),'k');





