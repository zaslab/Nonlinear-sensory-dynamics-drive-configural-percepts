

color = colormap;
color(256,:) = [0.1 0.8 0.2];
color([1:30, 40:60, 125:145],:) = color([125:145,1:30, 40:60],:);
ind = 170:200;
color(ind,:) = repmat([1,0.5,0], length(ind),1);

ind = 80:110;
color(ind,:) = repmat([1,0.1,0.1], length(ind),1);

ind = 200:240;
color(ind,:) = repmat([0.9,0.9,0.0], length(ind),1);


for i = 1:2
    if i == 1 
        var = explained;
    else
        var(1) = 100*(1-R(1));
        var(2) = 100*(R(1)-R(2));
        var(3) = 100*(R(2)-R(3));
    end
subplot(1,2,i)
xlabel(['PC1 ', num2str(round(var(1),2)),'%'])
ylabel(['PC2 ', num2str(round(var(2),2)),'%'])
zlabel(['PC3 ', num2str(round(var(3),2)),'%'])
set(gca, 'Colormap',color)
ax = get(gca);
ax.Children.SizeData = 50;
end

