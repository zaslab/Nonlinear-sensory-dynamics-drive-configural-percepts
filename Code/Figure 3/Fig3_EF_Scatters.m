%% shift all points so that stim A is around [0,0] and stim B is around [1,0]

range = 140;
[all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained1,combo_proj] = ...
    GetSynthWormScores(data_labels, data, feat_score,range);

combo_list = [2 4 8
    2 1 9
    4 5 10
    2 6 11
    4 6 12
    1 4 13];
neuron_pairs = [6,15;
    16,nan;
    1,nan;
    2,20;
    7,12;
    8,11;
    3,18;
    5,13;
    10,17;
    14,22;
    4,19;
    9,nan;
    21,nan];
mean_pos = {};
all_dist = [];
num_to_sample = 20;
dist_rat = [];
deg = [];
radi = [];
ab_dist_rat = [];
ab_deg = [];
ab_radi = [];
range = 140;
which_features = []; 

PCs = 1:3;

p = nan(6,27);

dist = [];
org_dist = [];

for combo = 1:size(combo_list,1)
    figure
    subplot(2,3,combo)
    hold on

    for n = 1

        neurons = [n, n+13];
        for ablated = 1
            conds = combo_list(combo,:);
            a = neurons;
            b = [];
            for i = 1:3
                b = [b, ((conds(i) - 1)*20 + 1) : (conds(i)*20);] ;
            end

            mean_pos{1} = [mean(all_score(b(1:num_to_sample),PCs),1);
                mean(all_score(b((num_to_sample+1):(num_to_sample*2)),PCs),1);
                mean(all_score(b((num_to_sample*2+1):(num_to_sample*3)),PCs),1)];
            [c3, points] = DrawRotatedCombos(all_score, PCs,combo_list,combo);
            plot3([0.5,c3(1,3)], [0,c3(2,3)], [0,c3(3,3)], 'k', 'LineWidth',2)

        end
        p1 = GetGMPval(points(1:20,:), points(21:40,:), points(41:60,:), PCs);
        p(combo,1) = p1(3);

    end
    close
    %%
    figure;
    prob_range = prob_ranges(combo,:);

    color = [repmat([1,0 0],20,1); repmat([0 0 1],20,1); repmat([0 0 0],20,1)];
    maxval(1) = DrawGaussianSphere(points(1:20,:),'red',limfactor1(:,:,i),prob_range );
    maxval(2) = DrawGaussianSphere(points(21:40,:),'blue',limfactor1(:,:,i),prob_range );
    maxval(3) = DrawGaussianSphere(points(1:40,:),'magenta',limfactor1(:,:,i),prob_range );
    c3 = c3';
    PlotAlignedCombo(points, c3, color)

    xlim(xlimit(combo,:))
    ylim(ylimit(combo,:))
    zlim(zlimit(combo,:))
    axis equal
    title(combo)
   
    comb_center = mean(c3(:,1:2),2);
    centers = [mean(points(1:20,1:3),1); mean(points(21:40,1:3),1); mean(points(41:60,1:3),1)];
    plot3([comb_center(1);centers(3,1)], [comb_center(2);centers(3,2)], [comb_center(3);centers(3,3)],'k','LineWidth',2)
    view([0,0])
    lightangle(40,35)
    %%
end



prob_ranges(1,:) = [1.75 3.5 7];
prob_ranges(2,:) = [0.4 0.8 1.2];
prob_ranges(3,:) = [0.2 1.6 9.5];
prob_ranges(4,:) = [0.3 1.3 2.3];
prob_ranges(5,:) = [1.5 3 4.5];
prob_ranges(6,:) = [1 2 6];

xlimit(1,:) = [-0.5 1.8];
ylimit(1,:) = [-0.5 0.5];
zlimit(1,:) = [-1 1];

xlimit(2,:) = [-1.2 1.5];
ylimit(2,:) = [-1.2 1.2];
zlimit(2,:) = [-1.2 1.2];

xlimit(3,:) = [-0.5 1.8];
ylimit(3,:) = [-0.5 0.5];
zlimit(3,:) = [-1 1];

xlimit(4,:) = [-0.5 1.8];
ylimit(4,:) = [-0.5 0.5];
zlimit(4,:) = [-1 1];

xlimit(5,:) = [-0.5 1.8];
ylimit(5,:) = [-0.5 0.5];
zlimit(5,:) = [-1 1];

xlimit(6,:) = [-0.8 1.8];
ylimit(6,:) = [-0.5 0.5];
zlimit(6,:) = [-1 1];
%
limfactor1 = cat(3,[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2]);
limfactor2 = cat(3,[3,3;2,2;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2]);
limfactor3 = cat(3,[7,7;2,2;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2],[2,2;3,3;2,2]);


function PlotAlignedCombo(points, c3, color)
center_size = 120;
scatter3(points(:,1), points(:,2),points(:,3),50,color,'filled')
hold on
scatter3(c3(1,1), c3(2,1),c3(3,1),center_size,[1 0 0],'filled')
scatter3(c3(1,2), c3(2,2),c3(3,2),center_size,[0 0 1],'filled')
scatter3(c3(1,3), c3(2,3),c3(3,3),center_size,[0 0 0],'filled')
plot3([1,-1], [0,0], [0,0],'k')
plot3([0,0], [-1,1], [0,0],'k')
plot3([0,0], [0,0], [1,-1],'k')
end
