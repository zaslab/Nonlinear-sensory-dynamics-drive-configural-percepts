
sample_num = 20;
range = 7*sample_num;
[all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained,combo_proj] = ...
                   GetSynthWormScores(data_labels, data, feat_score,range,'sample_num', sample_num);
valence = [zeros(1,40), ones(1,20), zeros(1,20), ones(1,40), ones(1,20)*2];
feat_worms = feature_worms;
feat_worms(:,sum(feat_worms == 0) > 10) = [];
dat_mat = [raw_acts];%, feat_worms];
norm_dat = NormalizeForPCA(dat_mat);
[score,combo_proj, combo_ids,coefforth , mu, explained] = PCAandProjection(dat_mat,ids,range);

options.verbose = 0;
options.dims = 1:26;
options.display = 0;
inds = [141:160
        161:180
        181:200
        201:220
        221:240
        241:260];

for i = 1
X_data = NormalizeForPCA([raw_acts(1:range,:)]');%[1:120,inds(i,:)],:)])';%, feat_worms]';%([1:140,inds(i,:)],:)';
D = L2_distance(X_data, X_data, 1); 
figure

[Y, R, E] = Isomap(D, 'k', 15, options); 

rows = [1:120, inds(i,:)];
rows = 1:range;
% scatter3(Y.coords{3,1}(1,rows), Y.coords{3,1}(2,rows), Y.coords{3,1}(3,rows), 50, ids(1:length(rows)), 'filled');
% scatter(Y.coords{3,1}(1,rows), Y.coords{3,1}(2,rows), 50, ids(1:length(rows)), 'filled');
subplot(1,2,1)
scatter3(all_score(rows,1), all_score(rows,2), all_score(rows,3), 50, ids(rows), 'filled');
axis equal
subplot(1,2,2)
scatter3(Y.coords{3,1}(1,rows), Y.coords{3,1}(2,rows), Y.coords{3,1}(3,rows), 50, ids(rows), 'filled');
axis equal
end

SetScatterColors

%%


plotflag = 1;
range = 140;
% [all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained,combo_proj] = ...
%         GetSynthWormScores(data_labels, data, feat_score,range,'sample_num', sample_num);
y = all_score;
% y(idx3,:) = shuf_proj;
combonames = data_cell(8:13,3);
figure
color = [repmat([1,0 0],20,1); repmat([0 0 1],20,1); repmat([0 0 0],20,1)];
for i = 1:size(combo_list,1)
    subplot(2,3,i)   
   
           % scatter3(y(:,1), y(:,2),y(:,3),50, ids, 'filled') 
           % alpha(0.2)
hold on
    c = combo_list(i,:);
            idx1 = ((c(1)-1)*20 + 1):c(1)*20;
            s1 = y(idx1,:);
            idx2 = ((c(2)-1)*20 + 1):c(2)*20;
            s2 = y(idx2,:);
            idx3 = ((c(3)-1)*20 + 1):c(3)*20;
            s3 = y(idx3,:);
plotdata = [s1;s2;s3];


% colormap(allcolor)
caxis([1 7])
hold on

[dist_rat(i), deg(i), radi(i)] = DrawLines(s1,s2,s3,color,plotflag);


title(combonames(i))
 grid on

 % PlotRotatedComboScatter
% xlabel(['PC1 ', num2str(explained(1)),'%'])
% ylabel(['PC2 ', num2str(explained(2)),'%'])
% zlabel(['PC3 ', num2str(explained(3)),'%'])

% xlim([min(score(:,1)) - 0.4, max(score(:,1)) + 0.4])
% ylim([min(score(:,2)) - 0.4, max(score(:,2)) + 0.4])
% zlim([min(score(:,3)) - 0.4, max(score(:,3)) + 0.4])
% 
% xlim([min(data(:,1)) - 0.4, max(data(:,1)) + 0.4])
% ylim([min(data(:,2)) - 0.4, max(data(:,2)) + 0.4])
% zlim([min(data(:,3)) - 0.4, max(data(:,3)) + 0.4])

end
