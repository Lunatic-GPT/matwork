function set_data_img(par_fname,par_vname,h)
% set_data_img(par_fname,par_vname,h)
% par_fname contains:
% (ufile,ofile,uroifile,oroifile)
% par_vname contains: % variable name to use inside each file
% (ufile,ofile,uroifile,oroifile)

% h is the handle;


if isfield(par_fname,'ufile') && exist(par_fname.ufile,'file')  
    ufile=par_fname.ufile;
    if isfield(par_vname,'ufile')
     udata=ri_d1(ufile,'','',par_vname.ufile);
    else
     udata=ri_d1(ufile);
    end
    
    setappdata(h,'udata',double(udata));
    setappdata(h,'ufile',ufile);
else
     setappdata(h,'udata',[]);
    setappdata(h,'ufile',[]);
end


if isfield(par_fname,'ofile') && exist(par_fname.ofile,'file')  
    ofile=par_fname.ofile;
    if isfield(par_vname,'ofile')
     odata=ri_d1(ofile,'','',par_vname.ofile);
    else
     odata=ri_d1(ofile);
    end
    
    setappdata(h,'odata',double(odata));
    setappdata(h,'ofile',ofile);
else

    if isappdata(h,'odata')  
     rmappdata(h,'odata');
   end

   if isappdata(h,'ofile')  
     rmappdata(h,'ofile');
   end
end



if isfield(par_fname,'uroifile') && exist(par_fname.uroifile,'file')  
    uroifile=par_fname.uroifile;
    if isfield(par_vname,'uroi')
     uroi=ri_d1(uroifile,'','',par_vname.uroi);
    else
     uroi=ri_d1(uroifile);
    end
    
    setappdata(h,'uroi',double(uroi));
    setappdata(h,'uroifile',uroifile);
else
    if isappdata(h,'uroi')
        rmappdata(h,'uroi');
    end
    if isappdata(h,'uroifile')
        rmappdata(h,'uroifile');
    end
end



if isfield(par_fname,'oroifile') && exist(par_fname.oroifile,'file')  
    oroifile=par_fname.oroi;
    if isfield(par_vname,'oroi')
     oroi=ri_d1(oroifile,'','',par_vname.oroi);
    else
     oroi=ri_d1(oroifile);
    end
    
    setappdata(h,'oroi',double(oroi));
    setappdata(h,'oroifile',oroifile);
else
    if isappdata(h,'oroi')
        rmappdata(h,'oroi');
    end
    if isappdata(h,'oroifile')
        rmappdata(h,'oroifile');
    end
end


showImages(getappdata(h,'handles'));


