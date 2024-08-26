function plot2(x,layout)

figure;
for i=1:prod(layout)
    if i>size(x,2)
        continue;
    end
subplot(layout(1),layout(2),i);


plot(x(:,i));
end

