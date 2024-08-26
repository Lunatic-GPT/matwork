function dim= showImages(handles)
    

   layout=str2num(get(handles.Layout,'String'));
   
   maxNumCols=layout(2);
   maxNumRows=layout(1);

    %set(handles.cancell,'value',0);

    brik=get(handles.cur_brik,'String');
    brik=str2double(brik);
   
    umdata=getappdata(handles.figure1,'umdata');    
    
    umbrik= str2double(get(handles.um_brik,'String'));
    
    if isappdata(handles.figure1,'udata2')
      udata=getappdata(handles.figure1,'udata2');
    else
      udata=getappdata(handles.figure1,'udata');
    end
    
    odata=getappdata(handles.figure1,'odata');
     
    if brik>size(udata,4)
        brik=size(udata,4);
        set(handles.cur_brik,'String',num2str(brik));
    end
    if brik<0
        brik=1;
        set(handles.cur_brik,'String',num2str(brik));
    end
    
 
    
    obrik=get(handles.obrik,'String');
    obrik=str2num(obrik);
    if obrik>size(odata,4)
        obrik=size(odata,4);
        
        set(handles.obrik,'String',num2str(obrik-1));
        
    end

    
    dim=size(udata);
    
   
    hActions = getappdata(handles.figure1,'actionItems');
    pMargin = 0.0;
    
  %  N=size(udata,3);
                       
    
    setappdata(handles.figure1,'panelMargin',pMargin);
    rPitch = (1-pMargin)/maxNumRows;
    setappdata(handles.figure1,'rPitch',rPitch);
    cPitch = 0.98/maxNumCols;
    axWidth = 0.95/maxNumCols;
    axHight = 0.88/maxNumRows;
    
    axesProp = {'dataaspectratio','Parent','PlotBoxAspectRatio','xgrid','ygrid'};
    axesVal = {[1,1,1],handles.uipanel1,[1 1 1],'off','off'};
    
    hAxes = getappdata(handles.figure1,'hAxes');
    
    for i=1:length(hAxes)
        j=1;
        slice=getappdata(hAxes(i),'slice');
        
        while ~isempty(findobj(hAxes,'Tag',roi_pos_name(slice,j)))
            hr=findobj(hAxes,'Tag',roi_pos_name(slice,j));
            api = iptgetapi(hr);
                
                 roi_pos = api.getPosition();
                 ud=get(handles.cropud,'String');
                 ud=str2num(ud);
                 lr=get(handles.croplr,'String');
                 lr=str2num(lr);
                 
                 roi_pos(:,1)=roi_pos(:,1)+lr(1);
                 
                 roi_pos(:,2)=roi_pos(:,2)+ud(1);
                 
                 setappdata(gcf,roi_pos_name(slice,j),roi_pos); 
                 delete(hr);  % will redraw later
                 j=j+1;
        end
    end

    %{
    sval=get(handles.slider1,'Value');
    sstep=get(handles.slider1,'SliderStep');
    if sstep(1)==0
        sstep(1)=3;
    end
    sos=round((1-sval)/sstep(1)*maxNumRows*maxNumCols);
    set(handles.first_slice,'String',sprintf('%d',sos+1));
    
    %}
    
    sos=get(handles.first_slice,'String');  % slice offset for the underlay
    sos=str2num(sos)-1;
    
       croplr=get(handles.croplr,'String');
       cropud=get(handles.cropud,'String');
       croplr=str2num(croplr);
       cropud=str2num(cropud);
       
    slice_step=get(handles.sliceStep,'String');
    
    slice_step=str2double(slice_step);
    
    for ind =1:maxNumRows*maxNumCols
    
        [c r] = ind2sub([maxNumCols maxNumRows],ind);
        x = (c-1)*cPitch;
        y = 1-pMargin-r*rPitch;
        
        if isempty(hAxes) || ind>length(hAxes)
            
        hAxes(ind) = axes( 'position', [x y axWidth axHight],axesProp,axesVal);
        
        elseif ind<=length(hAxes)
            set(hAxes(ind),'position', [x y axWidth axHight]);
        end
        
     current_slice=(ind-1)*slice_step+sos+1;
    
     odata_rs=getappdata(handles.figure1,'odata_rs');
     ind_u2o=getappdata(handles.figure1,'ind_u2o');
     
     o_voxsize=str2num(get(handles.o_voxsize,'String'));
     o_center=str2num(get(handles.o_center,'String'));
     
     u_voxsize=str2num(get(handles.u_voxsize,'String'));
     u_center=str2num(get(handles.u_center,'String'));
     ignoreVoxSize=get(handles.ignoreVoxSize,'Value');
        
 %% resample odata if necessary       
 odata2=[];
 
 o_orient=getappdata(handls.figure1,'o_orient');
 u_orient=getappdata(handls.figure1,'u_orient');
 
 if ~isempty(odata)
     %  doresample=getappdata(handles.figure1,'doresample');
  
     
     %if isempty(odata_rs)
         
         %if ~ignoreVoxSize && ((ndims(odata(:,:,:,1))~=ndims(udata(:,:,:,1)) || any(size(odata(:,:,:,1))~=size(udata(:,:,:,1)))) || any(o_voxsize~=u_voxsize) || any(o_center~=u_center))
                 
             %[odata_rs,ind_u2o]=reslice_extend(odata,o_voxsize,u_voxsize,size(udata(:,:,:,1)),u_center-o_center);
      if o_ref        
             
              odata_rs=reslice_with_orient(sz_ref,o_ref,o_in,d_in,use_nn);
             
             
         else
             odata_rs=odata;
             ind_u2o=0;
         end
         
         setappdata(handles.figure1,'odata_rs',odata_rs);
         setappdata(handles.figure1,'ind_u2o',ind_u2o);
     end
     
     if isempty(ind_u2o)
         ind_u2o=0;
     end
     
     if current_slice+ind_u2o>=1  && current_slice+ind_u2o<=size(odata_rs,3)
         odata2=odata_rs(:,:,current_slice+ind_u2o,obrik);
     end
      
 end  % end resampling odata
 
     matchOverlay = get(handles.ROIMatchOverlay,'Value');
     
     if matchOverlay
         roi=getappdata(handles.figure1,'oroi');  % roi match the overlay
         roi_rs=getappdata(handles.figure1,'oroi_rs');
         if isempty(roi_rs) && ~ignoreVoxSize && ((ndims(roi(:,:,:,1))~=ndims(udata(:,:,:,1)) || any(size(roi(:,:,:,1))~=size(udata(:,:,:,1)))) || any(o_voxsize~=u_voxsize) || any(o_center~=u_center))
             roi_rs=reslice_extend(roi,o_voxsize,u_voxsize,size(udata(:,:,:,1)),u_center-o_center,true);
             setappdata(handles.figure1,'oroi_rs',roi_rs);
         elseif isempty(roi_rs)
             roi_rs=roi;
         end
         
         roi=roi_rs;
       
     else
         roi=getappdata(handles.figure1,'uroi');   %roi match the underlay
     end
     
     
   if get(handles.checkbox2,'Value')  && exist('odata_rs','var')
       
       mdata=odata_rs;                 
    else
       mdata=getappdata(handles.figure1,'mdata');
    end    
    mbrik= str2double(get(handles.edit11,'String'));
   
    if mbrik>size(mdata,4)
        mbrik=size(mdata,4);
    end
       
        
        setappdata(hAxes(ind),'slice',current_slice);
        
        mdata2=[];
        if ~isempty(mdata)     %no resampling for mdata yet; can implement later if needed.
             if current_slice+ind_u2o>=1  && current_slice+ind_u2o<=size(odata_rs,3)
               mdata2=mdata(:,:,current_slice+ind_u2o,mbrik);
             end
        end
        
        
        roi2=[];
        if ~isempty(roi)
            
            if matchOverlay
                ind_u2o_roi=ind_u2o;
            else
                ind_u2o_roi=0;
            end
            
            if current_slice+ind_u2o_roi<=size(roi,3) && current_slice+ind_u2o_roi>=1
                roi2=roi(:,:,current_slice+ind_u2o_roi,:);
            end
            
        end
        
        if current_slice>size(udata,3) || current_slice <1
            udata2=0*udata(:,:,1,1);
        else
           udata2=udata(:,:,current_slice,brik);
        end
        
        if isempty(umdata)   
           plotImInAxis(hAxes(ind),udata2,odata2,mdata2,roi2,handles);
        else
            umrange=get(handles.um_range,'String');
            umrange=str2num(umrange);
            if current_slice<=size(udata,3) && current_slice >=1
                mask_tmp=double(umdata(:,:,current_slice,umbrik)>umrange(1) & umdata(:,:,current_slice,umbrik)<umrange(2));
            else
                mask_tmp=0.*udata2;
            end
            
             plotImInAxis(hAxes(ind),udata2.*mask_tmp,odata2,mdata2,roi2,handles);
            
        end
        j=1;
       while isappdata(handles.figure1,roi_pos_name(current_slice,j));
        roi_pos=getappdata(handles.figure1,roi_pos_name(current_slice,j));
          roi_pos(:,1)=roi_pos(:,1)-croplr(1);
          roi_pos(:,2)=roi_pos(:,2)-cropud(1);
          
          h=impoly(hAxes(ind),roi_pos);
          api=iptgetapi(h);
          api.setColor('green');
          rmappdata(  handles.figure1,roi_pos_name(current_slice,j));    
          set(h,'Tag',roi_pos_name(current_slice,j));
        j=j+1;
       end
    end
    
    for ind=maxNumRows*maxNumCols+1:length(hAxes)
        delete(hAxes(ind));
        
    end
    hAxes=hAxes(1:maxNumRows*maxNumCols);
    
    
   % set(hActions,'enable','on');
    setappdata(handles.figure1,'hAxes',hAxes);
    
    ufile=getappdata(handles.figure1,'ufile');
    if ~isempty(ufile)
        fig_name=sprintf('U:%s',ufile);
    end
  
  umfile=getappdata(handles.figure1,'umfile');
  if ~isempty(umfile)
      fig_name=sprintf('%s/UM:%s',fig_name,umfile);
  end
  
  ofile=getappdata(handles.figure1,'ofile');
  if ~isempty(ofile)
      fig_name=sprintf('%s/O:%s',fig_name,ofile);
  end
  
  mfile=getappdata(handles.figure1,'mfile');
  if ~isempty(mfile)
      fig_name=sprintf('%s/OM:%s',fig_name,mfile);
  end
  
  if ~matchOverlay
      uroifile=getappdata(handles.figure1,'uroifile');
      if ~isempty(uroifile)
          fig_name=sprintf('%s/uroi:%s',fig_name,uroifile);
      end
  else
      oroifile=getappdata(handles.figure1,'oroifile');
      if ~isempty(oroifile)
          fig_name=sprintf('%s/oroi:%s',fig_name,oroifile);
      end    
  end
  
  set(handles.figure1,'Name',fig_name);





function odata2=calc_odata(sz_udata,odata,handles)
        
        os=get(handles.centerShift,'String');
        os=str2num(os);
        odim=get(handles.odim,'String');
        odim=str2num(odim);
        if length(sz_udata)==2
            sz_udata(3)=1;
        end
        odata2=reslice(odata,eye(3),[1,1,1],-os,sz_udata,odim,1);
        
        
   function res=get_next_handle
        
        res=101;
        while ishandle(res)
            res=res+1;
        end
        
   function plotImInAxis(hAx,im,odata,mdata,roi,handles)
       
       if any(abs(imag(im(:)))>0)
           im=abs(im);
       end
       
        if any(abs(imag(odata(:)))>0)
           odata=abs(odata);
        end
       
         if any(abs(imag(mdata(:)))>0)
           mdata=abs(mdata);
         end
       
       if isempty(im)
           return;
       end
       axes(hAx);

      sf= get(handles.showFull,'Value');
       
       xlcur=xlim;
       ylcur=ylim;
   imageProp = { 'ButtonDownFcn'};
    
   scaling=get(handles.listbox2,'value');
     
  %  imageVal = { sprintf('openSpecificImage(''%s'',%d,%f)',fullfile(dir_name,fname),slice,scale(2)) };
   imageVal = { 'openSpecificImage2' };
       switch scaling 
        case 2
            tmp_range = get(handles.edit3,'String');
            
            scale = str2num(tmp_range);
        case 1
            
            scale = [min(im(:)),max(im(:))];
     
        otherwise
            error('unknown image intensity scaling type');
       end 
    
      nvox_thr=get(handles.nvox_thr,'String');
      nvox_thr=str2num(nvox_thr);
      
       nclr=50;
       under=(im-scale(1))*nclr/diff(scale);
    
       
              cm=get(handles.listbox4,'Value');
             
            switch cm
                case 1
                    cm_u=gray(nclr);
                case 2
                    cm_u=gray(nclr);
                    cm_u=flipud(cm_u);
                case 3
                  cm_u=jet(nclr);
                case 4
                    cm_u=jet(nclr);
                    cm_u=flipud(cm_u);
                case 5
                    cm_u=autumn(nclr);
                case 6
                    cm_u=autumn(nclr);
                    cm_u=flipud(cm_u);
                case 7
                    cm_u=winter(nclr);
                    cm_u=flipud(cm_u);
                case 8
                    cm_u=winter(nclr);
               
                otherwise
            end
            
        if get(handles.showOverlay,'Value') && ~isempty(odata)
             
             oscale=get(handles.edit12,'String');
             
             oscale=str2num(oscale);
             if length(oscale)~=2
                 warning('range should be two values');
                 oscale=[0,100];
             end
             
                 nclr=50;
             ov =nclr*(odata-oscale(1))/diff(oscale);
              if ~isempty(mdata)
               mscale=get(handles.mrange,'String');
               mscale=str2num(mscale);
               if length(mscale)==2
                 mask=(mdata>=mscale(1) & mdata<=mscale(2));
               elseif length(mscale)==4
                   mask=(mdata>=mscale(1) & mdata<=mscale(2)) | (mdata>=mscale(3) & mdata<=mscale(4));
               else
                  error('wrong threshold size'); 
               end
               
             else
                mask=true(size(under));
              end
             
            
             cm=get(handles.listbox3,'Value');
         
            switch cm
                case 1
                  cm_o=jet(nclr);
                case 2
                    cm_o=jet(nclr);
                    cm_o=flipud(cm_o);
                case 3
                    cm_o=autumn(nclr);
                case 4
                    cm_o=autumn(nclr);
                    cm_o=flipud(cm_o);
                case 5
                    cm_o=winter(nclr);
                    cm_o=flipud(cm_o);
                case 6
                    cm_o=winter(nclr);
                case 7
                    cm_o=gray(nclr);
                case 8
                    cm_o=gray(nclr);
                    cm_o=flipud(cm_o);
                otherwise
            end
            if nvox_thr>1
             mask=clusterize2(mask>0,nvox_thr);
            end
           [out,cm_out]=combine_over_under(under,ov,cm_u,cm_o,mask);
                     
        else
            out=under;
            cm_out=cm_u;
        end
       
               out2=out;
       nclr=get(handles.nroicolor,'String');
       nclr=str2num(nclr);
              cm=get(handles.listbox5,'Value');
             
           
            switch cm
                case 1
                  clr=jet(nclr);
                case 2
                    clr=jet(nclr);
                    clr=flipud(clr);
                case 3
                    clr=autumn(nclr);
                case 4
                    clr=autumn(nclr);
                    clr=flipud(clr);
                case 5
                    clr=winter(nclr);
                    clr=flipud(clr);
                case 6
                    clr=winter(nclr);
                case 7
                    clr=gray(nclr);
                case 8
                    clr=gray(nclr);
                    clr=flipud(clr);
                case 9
                    clr=[1 0 0;0 0 0];
                otherwise
                    clr=[0 0 0;1 0 0];
            end
            

   
       %[1,0,0;1,1,0;0,1,0;0,0,1;0,1,1;1,0,1];
       
       showROI=get(handles.showROI,'Value');
       roi_brik=get(handles.roi_brik,'String');
       roi_brik=str2num(roi_brik);
       
       ind=unique(roi(:));
       ind(ind==0)=[];
       ind(isnan(ind))=[];
       
       croplr=get(handles.croplr,'String');
       cropud=get(handles.cropud,'String');
       croplr=str2num(croplr);
       cropud=str2num(cropud);
       
       
       nd=get(handles.magn4roi,'String');
       nd=str2num(nd);
       

           
           
           
           
       out2 = out2(croplr(1)+1:end-croplr(2),cropud(1)+1:end-cropud(2));
       if nd>1
       out2=repmat2(out2,nd);
       end 
       if ~isempty(roi)&&showROI
           m2=out2*0;
       for i=1:length(ind)
           if nd>1  
               tmp=repmat2(roi(croplr(1)+1:end-croplr(2),cropud(1)+1:end-cropud(2),:,roi_brik)==ind(i),nd);
           else
               tmp=roi(croplr(1)+1:end-croplr(2),cropud(1)+1:end-cropud(2),:,roi_brik)==ind(i);
           end
               
         m=bwmorph(tmp,'remove');
        % cind=mod(ceil(ind(i)/nclr)-1,size(clr,1))+1;
          cind=mod(ind(i)-1,size(clr,1))+1;
         
         m2(m>0)=cind;
       end

        [out2,cm_out]=combine_over_under(out2,m2,cm_out,clr,m2);
       
       end
       
       
       
       %imshow(uint16(out2'),cm_out); % this somehow does not work anymore;
       %has to change to the following code to make it work.
       out2=uint16((out2));     
       out2(out2<1)=1;
       out2(out2>size(cm_out,1))=size(cm_out,1);    
       out3=cat(3,cm_out(out2',1),cm_out(out2',2),cm_out(out2',3));  
       out3=reshape(out3,[size(out2'),3]);
       imshow(out3);
        
        
        setappdata(hAx,'colormap',cm_out);
 
 %    set(h,imageProp{1},imageVal{1});
    axis(hAx,'image');
    axis(hAx,'off');   
    
    setappdata(hAx,'image',out2');
    
    tmp=out(croplr(1)+1:end-croplr(2),cropud(1)+1:end-cropud(2));
    setappdata(hAx,'rect_pos',[croplr(1)+1,cropud(1)+1,size(tmp,1),size(tmp,2)]);
    slice=getappdata(hAx,'slice');
    title(sprintf('%d',round(slice)));
    xlim([0,size(out2,1)]);
    ylim([0,size(out2,2)]);
    if ~sf
     xlim(xlcur);
     ylim(ylcur);
    end

  % drawnow;
    