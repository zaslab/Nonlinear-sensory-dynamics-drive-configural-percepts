%% Figure 1   

neur_nums = unique(data_labels(:,1));
all_neur_acts = nan(15,26);
for n = 1:length(neur_nums)
    for cond = 1:15
        for step = 1:2
                inds= find(data_labels(:,1) == neur_nums(n) & ...
                data_labels(:,2) == cond & ...
                data_labels(:,3) == step);
                mean_act = nanmean(data(inds));
                all_neur_acts(cond, (step-1)*13 + n) = mean_act;
        end
    end
end

colormap('default')
colors = colormap('parula');
jetcolors = colormap('jet');
ind_c = 220:256;
ind_jet = round(linspace(155, 246, length(ind_c)));
colors(ind_c,:) = jetcolors(190:226,:);
neurons = relevant_neurons(unique(data_labels(:,1)));
neurons = [neurons, neurons];
%%

mse = [];
range = 140;
feat_score = [];
sample_num=20;
acts_by_combo = GenerateComboActs(all_neur_acts, 1);
[all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained,combo_proj] = ...
    GetSynthWormScores(data_labels, data, feat_score,range,'sample_num', sample_num);


pts = acts_by_combo(1:end,1:(end-2));
idn = acts_by_combo(1:end,(end-1:end));
% idn(length(idn)/2+1:end,1) = idn(length(idn)/2+1:end,1) - 13;
pts(idn(:,1) == 1,:) = [];
idn(idn(:,1) == 1,:) = [];
pts(idn(:,1) == 14,:) = [];
idn(idn(:,1) == 14,:) = [];
idx = max(abs(pts(:,1:2)),[],2) < 0.1;
pts(idx,:) = [];
idn(idx,:) = [];
shuffle = 0;
[p1,p2,e1,e2,mu,alph,beta,mse1,mse2,ind_s,pts] = FitParabola(pts,shuffle);
mse(1,:) = [mse1,mse2];
shuffle = 1;
for i = 2:1000
[p1,p2,e1,e2,mu,alph,beta,mse1,mse2,ind_s,pts] = FitParabola(pts,shuffle);

mse(i,:) = [mse1,mse2];
end


%%
t = linspace(min(alph), max(alph), 200)';
y1 =  p1(1).*alph + p1(2);

linepts  =  p1(1).*t + p1(2);

% mse = rmse(beta, y1);
mean(mse)


y2 = p2(1).*alph.^2 + p2(2).*alph+ p2(3);

curvepts = p2(1).*t.^2 + p2(2).*t+ p2(3);
% mse = rmse(beta, y2);
mean(mse)
% pts_flat is N×3, and lies exactly in the PCA plane
figure;
plot(t, linepts,'m--', 'LineWidth',3)
hold on
plot(t, curvepts,'k--', 'LineWidth',3)
scatter(alph, beta,50, 'filled')
colormap(colors)
grid on
xlabel('PC1')
ylabel('PC2')
ylim([-0.45 0.45])

%%



t_hat = linspace(min(alph), max(alph), 200)';
curve3D = mu' + t_hat*e1' + (p2(1)*t_hat.^2 + p2(2)*t_hat + p2(3))*e2';
Line3D  = mu' + t_hat*e1' + (p1(1)*t_hat + p1(2))*e2';
figure; hold on; axis equal; grid on;
view([-100, 16])
scatter3(pts(:,1), pts(:,2), pts(:,3), 70, 'filled');

plot3(Line3D(:,1), Line3D(:,2), Line3D(:,3), 'm--', 'LineWidth',3);
plot3(curve3D(:,1), curve3D(:,2), curve3D(:,3), 'k-', 'LineWidth',3);
xlabel('X'); ylabel('Y'); zlabel('Z'); 
title('Direct 3D Parabola Fit');
colormap(colors)

%%
figure; histogram(mse(:,1), 'Normalization', 'probability',BinEdges=0.05:0.001:0.3)
hold on
histogram(mse(:,2),'Normalization', 'probability', BinEdges=0.05:0.001:0.3)
xlim([0.1, 0.25])
scatter(mse(1,1),0.01, 50,'b','*')
scatter(mse(1,2),0.01, 50,'r','*')


ylabel('Frequency')

xlabel('RMSE')
