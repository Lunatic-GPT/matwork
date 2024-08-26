
function iro=get_iro(data,Line,Par,Rep)

% 24:125

data=reshape(data,[size(data,1),32,length(data(:))/32/size(data,1)]);
Line=Line(1:32:end);
Par=Par(1:32:end);
Rep=Rep(1:32:end);

    nro=size(data,1);
d2=zeros(nro,max(Line)+1,max(Par)+1,32,'single');

for i=1:length(Line)

    if Rep(i)>min(Rep)
        break;
    end
    d2(:,Line(i)+1,Par(i)+1,:)= data(:,:,i);
               
end

tmp=sos(sos(sos(d2,4),3),2);


iro=find(tmp>0.05*max(tmp));
