%function [perf variable] = behavPerf();

% This function reads the stimulus output files generated by gk_DCOHmap()
% and evaluate the behavioral performance

% Xinmiao, 03/05/09

clear all;
expName = 'Sub10_CTR_RDK';
%expdir = ['/Users/Xinmiao/My_File_Folder/Smirnakis_lab/fMRI/CorticalDamage/fMRI_DATA/RDK_out/'];
%expdir = ['/Users/Xinmiao/Desktop/Sub10_CTR_RDK_102510/Sub10_CTR_RDK_102510_Mat/'];
expdir = pwd;
%expdir = 'temp/';
flist=dir('*RDKMotDetect*.mat');
load(flist(1).name);
expMode = params.mode;

nrun = size(flist,1);  %number of scans
TR=2;
ntrTrl = 12;% number of TRs per trial.
ncoh = params.nTestedLevel;
response = [];  % first column, right(1) or wrong(-1); second, variable

count = 0;
for ifl = 1:nrun
    load(flist(ifl).name);
    frametimes = params.frametimes;
    frametimes = frametimes(frametimes(:,2)>0, :);
    conditions = params.conditions;
    %condlist = params.condlist;
    keyTimes = params.keyTimes;
    keys = params.keys;
    keyPressCurDirCoh = params.keyPressCurDirCoh(params.keyPressCurDirCoh(:,1)~=-1, :);
    
    active = false;
    for it = 1:size(frametimes,1)
        if frametimes(it,2)>2 && active == false
            active = true;           
            count = count + 1;
            if expMode == 1
                variable = conditions(frametimes(it,2)).fractionSignalDots;
            else
                variable = conditions(frametimes(it,2)).directionRangeSignalDots;
            end
            stiDir = conditions(frametimes(it,2)).dotDirection;
            response(count, 2) = variable;
            timeBegin = frametimes(it, 1);  
            ind = find((keyTimes - timeBegin)<TR & (keyTimes - timeBegin)>0);
            if isempty(ind) 
                response(count, 1) = -1;
            else
                keyCodes = {};
                nKey = 0;
                for iKey = 1:length(ind)
                    if (keys(ind(iKey)) ~='5')
                        nKey = nKey + 1;
                        keyCodes{nKey} = keys(ind(iKey));
                    end
                end
                if isempty(keyCodes) || (length(keyCodes)>=2 && ~strcmp(keyCodes{1}, keyCodes{length(keyCodes)}) )
                    response(count, 1) = -1;
                elseif (stiDir ==0 && (keyCodes{1}=='3' || keyCodes{1}=='4') || ...
                        (abs(stiDir-3.1416)<0.02 && (keyCodes{1}=='1' || keyCodes{1}=='2')))
                    response(count, 1) = 1;
                else
                    response(count, 1) = -1;
                end
                
            end
        else
            if (frametimes(it,2))<=2
                active = false;
            end
        end
    end
end

vals = response(:, 2)
if isempty(strfind(expName, 'Coh'))
    vals = sortrows(vals);
else
    vals = sortrows(vals);
end
variable = [];
count = 0;
while ~isempty(vals)
    variable(count+1) = vals(1);
    vals = vals(vals~= variable(count+1));
    count = count+1;
end


for ind = 1:length(variable)
    perf(ind) = sum(response(:,2)==variable(ind) & response(:,1)==1)/sum(response(:,2)==variable(ind))*100;
end
%perf = (sortrows(perf'))';

h=figure;
plot(variable, perf, 'o-', 'LineWidth', 2);
fname = flist(1).name;
fname = fname(length(fname)-13:length(fname)-4);
%title(fname, 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
if expMode == 1
    xlabel('coherence', 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
else
    xlabel(['dir range ',expName], 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
end
ylabel('percent correct', 'FontName', 'Arial', 'FontSize', 20, 'FontWeight', 'bold');
set(gca,'FontSize',18,'LineWidth',2,'TickLength',[0.02 0.02])
%axis([-0.05 1.05 48 100])
title(num2str((sortrows(variable',1))'))
box on

if isempty(strfind(expName, 'Coh'))
    variable = sortrows(variable, -1);
end

resultName = ['perfResult', expName];
eval(['save ', resultName, ' variable perf']);

    

