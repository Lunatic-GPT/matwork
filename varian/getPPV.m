function vals = getPPV( ppName, ppPath )
% Get the values for parameter name in procpar file path
% Usage: vals = getPPV( ppName, ppPath )

% fn = 'I_t.fid/procpar'

fp = fopen( ppPath, 'r');
done = 0;
vals = [];

while( done == 0 )
   line = fgetl(fp);
   if (line == -1)
      done = 1;
   else
      if ~isletter(line)
         error( 'bad format')
      else
         [name, attr] = strtok(line);
         attr = str2num(attr);

         % Read in the values
         line = fgetl(fp);
         [cnt, parm] = strtok(line);
         cnt = str2num(cnt);

         % REAL_VALS
         if (attr(2) == 1) 
            vals = str2num( parm );
         % STRING_VALS
         else 
            vals = dbl_quote_extract( parm );
            while( size(vals,1) ~= cnt )
               line = fgetl(fp);
               vals = char(vals, dbl_quote_extract( line ) );
            end
         end

         if (strcmp(name, ppName))
            break;
			else
				vals=[];
			end

         % Read in the enums
         enum_line = fgetl(fp);
      end
   end
end

fclose( fp );


function outStr = dbl_quote_extract( inStr )
%dbl_quote_extract - Extract String from between pair of Double Quotes

dqIdx = findstr(inStr, '"');
dqCnt = size(dqIdx,2);

if ((dqCnt == 0) | (1 == mod(dqCnt,2)))
   error( 'Bad string double quote balance')
end 

for idx = 1:2:dqCnt
   off = [dqIdx(idx) + 1, dqIdx(idx+1) - 1];
   if (idx == 1)
      outStr = inStr(off(1):off(2));
   else
      outStr = char(outStr, inStr(off(1):off(2)));
   end
end

