function sdt2afni(fname)

   prefix = strtok(fname,'.');
   brikdata=rdSdt(prefix);
   fid=fopen([prefix,'.spr'],'r');
   a=textscan(fid,'%s%s%s%s%s');
   id = find(strcmp(a{1},'interval:'));
   info.ORIENT_SPECIFIC = [0  2 4];
   info.ORIGIN = [0,0,0];
   info.DELTA = str2double([a{2}(id),a{3}(id),a{4}(id)]);
   info.TYPESTRING = '3DIM_HEAD_ANAT';
   info.SCENE_DATA = [0,2,0];
   info.HISTORY_NOTE = 'sdt2afni';
   WriteBrikEZ(brikdata,info,'',prefix);
   
   