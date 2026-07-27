% Change scRNA data
T = readtable(['102_DamID_AlisonEala.txt']);
list = [T.inLADs];
lia1 = find(list == 1);
X = T{lia1,["RNA"]};

lia2 = find(list == 0);
Y = T{lia2,["RNA"]};

%heatmap plot
figure(1)
% subplot(1,2,1)
imagesc(sort(X,'descend'));
ColorMap = load('viridis.txt');
colormap(ColorMap);
% ColorMap = load('RedBlue.txt');
% ColorMap = flipud(ColorMap);
% colormap(ColorMap/255);
clim([25 100])
xticks([])
title('active genes ON LADs')
set(gca, 'LineWidth', 2, 'FontSize', 12, 'FontName', 'Arial');
colorbar

figure(2)
% subplot(1,2,2)
imagesc(sort(Y,'descend'));
ColorMap = load('viridis.txt');
colormap(ColorMap);
% ColorMap = load('RedBlue.txt');
% ColorMap = flipud(ColorMap);
% colormap(ColorMap/255);
clim([25 100])
xticks([])
title('active genes OFF LADs')
set(gca, 'LineWidth', 2, 'FontSize', 12, 'FontName', 'Arial');
colorbar

%pie chart
on = sort(X,'descend');
off = sort(Y,'descend');
on100 = length(on(on<=100));
on200 = length(on(on<=200 & on>100));
on_large = length(on(on>200));

off100 = length(off(off<=100));
off200 = length(off(off<=200 & off>100));
off_large = length(off(off>200));

figure(3)
pie([on100 on200 on_large]);

figure(4)
pie([off100 off200 off_large]);
