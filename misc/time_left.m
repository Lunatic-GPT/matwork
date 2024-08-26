function time_left(ind,total,dur,ind2,total2)
%  time_left(ind,total,dur,ind2,total2)
% dur is the duration up to the current iteration.
% ind: 1 based
% ind2: the slow loop count

if ~exist('ind2','var')
    dur=dur/ind;
    fprintf('Time: %d/%d - passed/remaining/total time: %d/%d/%d s\n',...
            ind,total,round(dur*ind),round((total-ind)*dur),round(total*dur));
        
else
    
    ind=ind+(ind2-1)*total;
    dur=dur/ind;    
    total=total*total2;
    fprintf('Time: %d/%d - passed/remaining/total time: %d/%d/%d s\n',...
            ind,total,round(dur*ind),round((total-ind)*dur),round(total*dur));
    
end
  