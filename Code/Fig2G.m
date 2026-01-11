
%%
load('dasds_data.mat')

%%
n1 = 2;
n2 = 20;
step_num = 2;
figure
   counter = 1;
   asj_acts = nan(50,9);
   area= nan(50,9);
for i = 1:size(dasds_data,1)

subplot(3,3,counter)
    [asj_acts(:,counter), area(:,counter)]= plotheatmap(...
                dasds_data{i,1},dasds_data{i,4}, ...
                dasds_data{i,3},n1,n2,step_num);
    counter = counter +1;
end


figure
plot_titles = {'WT', 'unc13', 'unc31'};
counter = 1;
all_ps = [];
for i = 1:3
 
subplot(3,1,i)
   
act_data = asj_acts(:, counter:counter+2);
% act_data = area(:, counter:counter+2);

UnivarScatterOriginal(act_data );
ylim([-0.1 1.5])
title(plot_titles(i))
    counter = counter +3;
    xticklabels({'DA', 'SDS','combo'})
end

%%
all_ps = [];
counter = 1;
for i = 1:3
    act_data = asj_acts(:, counter:counter+2);
    act_data = act_data(:);

[~,all_ps(i,1)] = ttest2(asj_acts(:,counter), asj_acts(:,counter+1));
[~,all_ps(i,2)] = ttest2(asj_acts(:,counter), asj_acts(:,counter+2));
[~,all_ps(i,3)] = ttest2(asj_acts(:,counter+1), asj_acts(:,counter+2));
counter = counter + 3;
end


%%
function [mean_step, area] = plotheatmap(data, steps, strain,n1,n2,step_num)


curr = permute(cat(3,data(n1,:,:),data(n2,:,:)),[3 2 1]);
% curr = permute(data(2,:,:),[3 2 1]);
step = permute(cat(3,steps(n1,step_num,:), steps(n2,step_num,:)), [3,2,1]);
% step = permute(steps(2,2,:), [3,2,1]);
traces = [];
for i = 1:size(step,1)
    traces(i,:) = smooth(curr(i,(step(i)-20):(step(i)+100)));
    traces(i,:) = traces(i,:) - mean(traces(i,12:20));
end
traces(isnan(curr(:,2)),:) = [];
mean_step = mean(traces(:, 22:40),2);
[~,ind] = sort(mean_step, 'descend');
imagesc(traces(ind(1:min([30,size(traces,1)])),:))

title(strain)
yticklabels([])
colormap jet
caxis([0,1])
[max_step, max_loc] = max(traces(:, 19:(end-40)),[],2);
area = [];
for t = 1:size(traces,1)
    mean_step(t) = mean(traces(t, (18 + max_loc(t)):(19 + max_loc(t)+2)))...
                    - mean(traces(t, 12:20));
    area(t) = trapz(1:40, traces(t,20:59)); %length(traces)-19
end
mean_step = max(traces(:, 20:end-40),[],2)-mean(traces(:, 12:20),2);
mean_step = unique(mean_step);
% mean_step = mean(traces(:, 20:60),2)-mean(traces(:, 12:20),2);
mean_step(end+1:50) = nan;
area(end+1:50) = nan;
end
