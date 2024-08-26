function par = parsejcamp(filename)
% PARSEJCAMP: Bruker JCAMP-DX Parameter File Parser
%      Usage: par = parsejcamp(filename);
%             par : structure with parsed values
%             filename : filename (with path from pwd) of parameter file
%
%      Parsejcamp is meant to extract most parameter values from a Bruker
%      JCAMP-DX version 4.24 file.  All the values are returned in a
%      structure with elements named identically to those in the parameter
%      file.  Currently, parsejcamp will skip multi-dimensional arrays and
%      structures in the JCAMP file.  This functionality may be added in
%      the future.


% Author: J. Luci, jeffrey.luci@utexas.edu, 26 March, 2013.
%
% License and terms of use:
% License to use and freely distribute this code is granted on the
% following conditions.  All comments are to be left unedited.
% All changes to the code should be thoughtfully commented, 
% described in the VERSIONS section of comments, and that the
% original author be made aware of all changes and improvements.

%VERSIONS:
% v1.0 (26 March, 2013): Released

% KNOWN BUGS
%
% v1.0:
% Any attempt to parse a file with XY data points saved in a parameter will
% result in an error.  This mostly only affects the RF files from MRI
% scans, but spectrum and peaklist information from NMR data will also
% yield the same error.


fid = fopen(filename, 'rt');
if fid == -1
    error('parsejcamp:fileNotFound', ['File: ' filename ' could not be found.']);
end
fileText = fread(fid, inf, 'uint8=>char')';
fclose(fid);

curNumElem=0;

poundPoundDollar = strfind(fileText, '##$');
poundPound = strfind(fileText, '##');
poundPound = setdiff(poundPound, poundPoundDollar);  %remove ##$ lines
poundPound = poundPound(1:numel(poundPound)-1);  %remove ##END
%dollarDollar = strfind(fileText, '$$'); %included only for future use of
                                         %comments, but not used now
returns = strfind(fileText, char(10));
equals = strfind(fileText, '=');

%parse required JCAMP-DX values
for ii = 1:numel(poundPound)
    %prepare to parse a single line of text
    %Step 1: Determine variable name
    startOfVarName = poundPound(ii)+2;
    eqLocation = find(equals > startOfVarName);
    endOfVarName = equals(eqLocation(1))-1;
    varName = fileText(startOfVarName:endOfVarName);
    
    %Step 2: Get value
    returnLocation = find(returns > startOfVarName);
    value = fileText(equals(eqLocation(1))+1:returns(returnLocation(1))-1);
    
    %Step 3: Determine value datatype and assign to element in struct
    if istxtnum(value)
        eval(['par.' varName '=' value ';']);
    else
        eval(['par.' varName '=''' value ''';']);
    end
end

if par.JCAMPDX ~= 4.24
    error('parsejcamp:IncorrectVersion', ['JCAMP file version is ' num2str(par.JCAMPDX) '.' ...
                                           char(10) 'Only v.4.24 is supported.']);
end

%parse ParaVision values
for ii = 1:numel(poundPoundDollar)
    %Step 1: Determine variable name
    startOfVarName = poundPoundDollar(ii)+3;
    eqLocation = find(equals > startOfVarName);
    endOfVarName = equals(eqLocation(1))-1;
    varName = fileText(startOfVarName:endOfVarName);
    returnLocation = find(returns > startOfVarName);
    value = fileText(equals(eqLocation(1))+1:returns(returnLocation(1))-1);
    
    %Step 2: Determine if single value, if so, assign
    if ~strcmp(fileText(equals(eqLocation(1))+1), '(')
        if istxtnum(value)
            eval(['par.' varName '=' value ';']);
        else
            eval(['par.' varName '=''' value ''';']);
        end
        
        %Step 3: Determine if 1D Array (including strings), if so, assign
    elseif ~sum(strfind(value, ','))
        numElemInArray = str2double(value(2:numel(value)-1));  %strip off parenth
        value = fileText(returns(returnLocation(1))+1:returns(returnLocation(2))-1);
        if istxtnum(value(1)) || strcmp(value(1), '-') %then it's a numeric array
            eval(['par.' varName '=[' value '];']);
            eval(['curNumElem = numel(par.' varName ');']);
            step = 2;
            while curNumElem < numElemInArray  %loop to account for multi-line arrays
                appendedValue = fileText(returns(returnLocation(step))+1:returns(returnLocation(step+1))-1);
                eval(['par.' varName '=[par.' varName ' ' appendedValue '];']);
                step=step+1;
                eval(['curNumElem = numel(par.' varName ');']);
            end
        else
            if strcmp(value(1), '<')     %strip off <> if present
                value = value(2:numel(value)-1);   
            end
            eval(['par.' varName '=''' value ''';']);
        end   
    end
end

    function trueFalse = istxtnum(str)
        trueFalse = logical(~isnan(str2double(str)));
    end


end