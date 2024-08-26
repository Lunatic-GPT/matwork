 function plot_TimeCourse
 
 gh=gcf;
     h=datacursormode(gcf);

  %  set(h,'Enable','on');
        s=getCursorInfo(h);

        udata=getappdata(gcf,'udata');
        
        
       xdata=findobj(gcf,'Tag','xdata');
       xdata=get(xdata,'Value');
       
       logy=findobj(gcf,'Tag','logy');
       logy=get(logy,'Value');
       
        odata=getappdata(gcf,'odata');
        if ~isempty(odata)
        obrik=findobj(gcf,'Tag','obrik');
        obrik=get(obrik,'String'); 
        obrik=str2num(obrik);
        odata=odata(:,:,:,obrik);
        end
    slice=getappdata(gca,'slice');
    roi=getappdata(gcf,'roi');
    
       hlr=findobj(gcf,'Tag','croplr');
       croplr=get(hlr,'String');
       
       hup=findobj(gcf,'Tag','cropud');
       cropud=get(hup,'String');
       
       croplr=str2num(croplr);
       cropud=str2num(cropud);
       
             plot_ts=findobj(gcf,'Tag','plotts');
             
             hind=findobj(gcf,'Tag','ts_brik');
             tFT=findobj(gcf,'Tag','tFT');
             tFT=get(tFT,'Value');
             hind=get(hind,'String');
             hind=str2num(hind);
             if length(hind)==1
                 hind=ones(1,2)*hind;
             end
             x=s.Position(1);
             
             y=s.Position(2);

             
             plot_ts=get(plot_ts,'Value');
             if plot_ts
            figure(getappdata(gcf,'handle_ts'));
            
            disp('....................................................');
            
             for i=1:3
                 if (hind(2)-hind(1)>0)
                 fprintf('tSNR= ');
                 end
                 
                 for j=1:3
                  
                  x1=x-2+j;
                  y1=y-2+i;
                  x1=mod(x1-1,size(udata,1))+1+croplr(1);
                  y1=mod(y1-1,size(udata,2))+1+cropud(1);
                  subplot(3,3,(i-1)*3+j);
                  
                  ts=squeeze(udata(x1,y1,slice,hind(1):hind(2)));
                  if (hind(2)-hind(1)>0)
                      if ~isempty(odata)
                       fprintf('%4.1f (u:%4.1f+-%4.1f; o:%4.1f);   ',mean(ts)/std(ts),mean(ts),std(ts),odata(x1,y1,slice));
                      else
                          fprintf('%4.1f (u:%4.1f=-%4.1f );   ',mean(ts)/std(ts),mean(ts),std(ts));
                      end
                  end
                 
                  if tFT
                      ts=abs(fft(ts-mean(ts),[],1));
                      ts=ts(1:ceil(end/2));
                  end
                  
                  if xdata
                      
                      xd=load('xdata.mat');
                      xd=xd.x;
                      
                      plot(xd,ts,'ro-');
                      
                  xlim([min(xd),max(xd)]);
                  else
                  
                     plot(ts);
                     
                  xlim([0,length(ts)]);
                  end
                 
                  if logy
                      set(gca,'yscale','log');
                  else
                      set(gca,'yscale','linear');
                  end
                 %    title(num2str(obrik));
                     
                  rg=min_max(squeeze(ts));
                  dlt=rg(2)-rg(1);
                  if dlt==0
                      dlt=0.01;
                  end
                  ylim([rg(1)-0.01*dlt,rg(2)+0.01*dlt]);
                  
                 end
                 fprintf('\n');
             end
             end
             figure(gh);
             