function [vec,b_unique]=read_diffvectors(fname)

fid=fopen(fname,'r');


n=fscanf(fid,'[directions=%d]\n');

d=fgetl(fid);
d2=fgetl(fid);


    vec=fscanf(fid,'Vector[%f] = ( %f, %f, %f )\n');
    fclose(fid);
    
vec=reshape(vec,[4,n]);

vec=vec(2:end,:)';

figure;




b=sos(vec,2);
b=b.^2;

b_unique=uniquetol(b,0.001);

for i=1:length(b_unique)

    ind=abs(b-b_unique(i))<0.001;
    fprintf('b=%f: %d\n',b_unique(i),sum(ind));

end

b_unique(b_unique<0.01)=[];

for iview=1:2
sym={'r.','b.','k.','c.','g.','y.'};
subplot(1,2,iview);

[x,y,z] = sphere(128);
surf(x,y,z,'FaceColor', 'none','EdgeColor',0.7*[1,1,1]);

hold on;
    
for i=1:length(b_unique)

    ind=abs(b-b_unique(i))<0.05;
    
    if b_unique==0
        continue;
    end
    
    
    x=vec(ind,:)/sqrt(b_unique(i));
    
     [theta,phi]=unitVec2thetaPhi(x');
     
     if iview==1
       sel= theta<=pi/2;
     else
         sel= theta>pi/2;
     end
   
    plot3(x(sel,1),x(sel,2),x(sel,3),sym{mod(i-1,6)+1},'MarkerSize',8);
    
end

view([0,0,(-1)^iview]);
end

