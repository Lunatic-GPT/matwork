function plot_retro_time_course_1voxroi(phase_file,mask_file,venc)
%m=ri('roi_vessel_mask_scan26.mat');
%m=ri('roi_vessel_mask_scan26_anteriorOnly.mat');
m=ri(mask_file);

b=ri(phase_file);
b=b*venc/180;


%% assume different baseline


vox_count=0;
v_all=[];    % assume different baseline
v_all_same=[];  % same baseline
for i=1:size(m,1)
    for j=1:size(m,2)
        for k=1:size(m,3)
            if m(i,j,k)==0
                continue;
            end
            
            m_bg=0*m;
            m_bg(i-4:i+4,j-4:j+4,k)=1;
            m_bg(i-1:i+1,j-1:j+1,k)=0;
            
            
            v=squeeze(b(i,j,1,:))-vec(mean_roi(b,m_bg));
            v_all=cat(1,v_all,v);
            
            v_all_same=cat(1,v_all_same,b(i,j,1,:));
            vox_count=vox_count+1;
        end
    end
end

v_all=reshape(v_all,[size(b,ndims(b)),vox_count]);

v_all_same=reshape(v_all_same,[size(b,ndims(b)),vox_count]);

v_all=circshift(v_all,[3,0]);

figure;errorbar(1:size(b,ndims(b)),-mean(v_all,2),std(v_all,[],2)/sqrt(vox_count));%

 xlim([1,10]);
set(gca,'FontSize',12);
    ylabel('Velocity (cm/s)');
    xlabel('Cardiac Phase');


figure;errorbar(1:size(b,ndims(b)),mean(v_all_same,2),std(v_all_same,[],2)/sqrt(vox_count));

 xlim([1,10]);
set(gca,'FontSize',12);
    ylabel('Velocity (cm/s)');
    xlabel('Cardiac Phase');


    
