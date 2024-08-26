function line_label(clr,spacing,hor_ver)

figure;
for i=1:size(clr,1)

    if hor_ver==1 %horizontal
       x= spacing*i;
       y=1;
    else
        x=1;
        y=spacing*i;
    end
    plot([x,x+0.3],y*[1,1],'LineWidth',3,'Color',clr(i,:));

hold on;
axis off;
end

set(gca,'YDir','reverse');

savetiffc('Line_label');
