data = readtable('BehaviorData.csv');

%%

day_means = groupsummary(data,["Date","Stim","ConcNum","Conc",...
                              "Component1","Component1Conc",...
                              "Component2","Component2Conc"],"mean", "CI");

%%
day_iaa = day_means(strcmp(day_means.Stim, 'IAA'),["Conc", "mean_CI","Date"]);
day_da = day_means(strcmp(day_means.Stim, 'DA'),["Conc", "mean_CI","Date"]);
day_na = day_means(strcmp(day_means.Stim, 'NaCl'),["Conc", "mean_CI","Date"]);


%%


day_means.Properties.VariableNames{day_means.Properties.VariableNames=="mean_CI"} = 'CI';
day_means.PlateNum = ones(size(day_means,1),1);


results = analyze_chemo_integration(day_means, 'DoBootstrap', true, 'B', 2000, 'Seed', 13);


%%
results = compare_pairs_loewe_lme1(results, 'UseDayMeans', true);
head(results.loewe_pairwise_lme)

%%


mask = string(results.combo_detail.Pair)=="DA — IAA —";
overlay = results.combo_detail(mask,:);

idx = mask;
results.combo_detail.p_pred_loewe(idx) = arrayfun(@(a,b,x,y) ...
    loewe_effect_at_dose(results.calib_models, a, b, x, y), ...
    results.combo_detail.A(idx), results.combo_detail.B(idx), ...
    results.combo_detail.cA_half(idx), results.combo_detail.cB_half(idx));
% Then call the plot again (it will show white prediction markers).
di_dist_day = plot_loewe_zero_surface_auto(results, overlay, ...
                    'UseDayMeans', true, 'PadP', 0.03, 'PlotScatter', true);
view(18,18)
colorbar("Ticks", 0:0.2:1)
caxis([0 1])
h = gca;
h.YLim = [10^-5, 10^-3];
h.ZLim = [0.0, 1];
h.XLim = [10^-6, 10^-2];
zticks(0:0.2:1)
zticklabels([])
yticklabels([])
xticklabels([])

mask = string(results.combo_detail.Pair)=="DA — NaCl —";
overlay = results.combo_detail(mask,:);


idx = mask;
results.combo_detail.p_pred_loewe(idx) = arrayfun(@(a,b,x,y) ...
    loewe_effect_at_dose(results.calib_models, a, b, x, y), ...
    results.combo_detail.A(idx), results.combo_detail.B(idx), ...
    results.combo_detail.cA_half(idx), results.combo_detail.cB_half(idx));
% Then call the plot again (it will show white prediction markers).
dn_dist_day = plot_loewe_zero_surface_auto(results, overlay,...
    'UseDayMeans', true, 'PadP', 0.03, 'PlotScatter', true);
view(18,18)
colorbar("Ticks", 0:0.2:1)

caxis([0 1])
g = gca;
zticks(0:0.2:1)
zticklabels([])
yticklabels([])
xticklabels([])

g.ZLim = [0.0, 1];
g.XLim = [10^-1, 10^1];
g.YLim = [10^-5, 10^-3];
%%

figure
colors = [0, 0.4470, 0.7410;
          0.8500, 0.3250, 0.0980];
b = bar([1,2], [mean(di_dist_day(:,1)), mean(dn_dist_day(:,1))]);
b.FaceColor = 'flat';
b.CData(2,:) = colors(2,:);
alpha(0.3)
hold on 
length_di = size(di_dist_day,1);
length_dn = size(dn_dist_day,1);

jitter_di = randn(1,length_di)*0.1;
jitter_dn = randn(1,length_dn)*0.1;

scatter(ones(1,length_di)+jitter_di,di_dist_day(:,1),[], colors(1,:),'filled')
scatter(2*ones(1,length_dn)+jitter_dn,dn_dist_day(:,1), [],colors(2,:), 'filled')

di_std = std(di_dist_day(:,1));
dn_std = std(dn_dist_day(:,1));
errorbar([mean(di_dist_day(:,1)), mean(dn_dist_day(:,1))], ...
    [di_std, dn_std],LineStyle="none")

xticklabels({'DA + IAA', 'DA + NaCl'})
yticks(0:0.1:0.2)
ylabel('Distance from Loewe surface (CI)')



