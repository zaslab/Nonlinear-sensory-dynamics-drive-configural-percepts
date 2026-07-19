function [all_y, all_paired_ps,diffs,res]= RunComboPvals(data_labels, rawdata,...
    feat_score, feat_pcs,combos,range,maxPC, max_loop,varargin)
ablate_neurons = 0;
for i=1:1:length(varargin)
    if strcmp(varargin{i},'ablate_neurons')
        ablate_neurons = 1;
        ablate_struct = varargin{i+1};
        neurons = ablate_struct.neurons ;
        combo_list = ablate_struct.combo_list;
        combo = ablate_struct.combo;
    end

end
error_counter = 1;
reg = 0.001;
all_mean_vals = [];
combo_sigs  = [];
all_y = [];
res = [];
class_count= {};
neuron = 0;

all_means = [];

all_ps = [];
all_paired_ps = [];
PCs = [1:maxPC];

mean_vals = [];
all_most_like = [];
diffs = [];
for loop = 1:max_loop
    goagain = 1;

    while goagain
        try
            goagain = 0;
            [~, osm6_acts, ~, ids,feature_worms,~,~,~] = ...
                GetSynthWormScores(data_labels, rawdata, feat_score,range,feat_pcs);
            feature_worms(:, sum(feature_worms == 0,1) > 10) = [];

            data = [osm6_acts];

            mean_data = nanmean(data, 1);
            std_data = nanstd(data, [], 1);
            for i = 1:size(data,1)
                data(i, :) = (data(i, :) - mean_data) ./ std_data;
            end
            [coefforth,score,~,~,~,mu]= pca(data(1:range,:));
            combo_proj = (data((range+1):end,:)-mu) * coefforth;
            combo_ids = ids((range+1):end);

            all_score = [score; combo_proj];
            if ablate_neurons
                all_score= AblateNeurons(neurons,data,coefforth,mu, combo_list, combo,range);
            end
            all_combo_probs = [];
            all_stim_probs = [];
            pdfValues = [];
            pdfValues2 = [];
            pdfValues3 = [];
            other_pdfValues = [];
            most_like = [];
            curr_mean_vals = [];

            for i = 1:size(combos,1)


                stim1 = all_score (ids == combos(i,1),PCs);
                stim2 = all_score (ids == combos(i,2),PCs);
                stim3 = all_score (ids == combos(i,3),PCs);

                gm1 = fitgmdist([stim1(:,PCs)], 1, "SharedCovariance",true,"CovarianceType","diagonal", "RegularizationValue",reg );
                gm2 = fitgmdist([stim2(:,PCs)], 1, "SharedCovariance",true,"CovarianceType","diagonal", "RegularizationValue",reg );


                new_combined = [];
                for pt = 1:20
                    new_combined(pt,:) = mean([stim1(randi(20),:); stim2(randi(20),:)],1);
                end
                gm4 = fitgmdist(new_combined(:,PCs), 1,"SharedCovariance",true,"CovarianceType","diagonal", "RegularizationValue",reg);

                combined_stims = [stim1(:,PCs); stim2(:,PCs)];
                gm3 = fitgmdist(combined_stims(:,PCs), 1,"SharedCovariance",true,"CovarianceType","diagonal", "RegularizationValue",reg);
                pdfValues3 = [ pdf(gm1, stim3(:,PCs)),  pdf(gm2, stim3(:,PCs)), pdf(gm3, stim3(:,PCs))];
                pdfValues4 = [ pdf(gm1, stim3(:,PCs)),  pdf(gm2, stim3(:,PCs)), pdf(gm3, stim3(:,PCs)), pdf(gm4, stim3(:,PCs))];


                pdfValues(:,1) = pdf(gm1, stim1(:,PCs));
                pdfValues(:,2) = pdf(gm2, stim1(:,PCs));
                pdfValues(:,3) = pdf(gm1, stim2(:,PCs));
                pdfValues(:,4) = pdf(gm2, stim2(:,PCs));
                pdfValues(:,5) = pdf(gm1, stim3(:,PCs));
                pdfValues(:,6) = pdf(gm2, stim3(:,PCs));
                pdfValues(:,7) = pdf(gm3, stim1(:,PCs));
                pdfValues(:,8) = pdf(gm3, stim2(:,PCs));
                all_stim_probs = [all_stim_probs ; pdfValues ];
                all_combo_probs = [all_combo_probs; pdfValues4];%(:,5:6)];
                all_combo_probs(all_combo_probs == 0) = 10e-200;
                curr_mean_vals(i,:) = mean(log(all_combo_probs),1);

                [~,most_like(:,i)] = max(pdfValues4,[],2);

                res(i,:,loop) = [mean(residuals_to_line(stim1, mean(stim1,1), mean(stim2,1))),...
                    mean(residuals_to_line(stim2, mean(stim1,1), mean(stim2,1))),...
                    mean(residuals_to_line(stim3, mean(stim1,1), mean(stim2,1)))];
            end
            all_most_like = cat(1, all_most_like, most_like);

        catch ME
            disp('Error')

            goagain = 1;
            error_counter = error_counter + 1;
            rethrow(ME)
        end
    end
    pvals = [];
    paired_p = [];
    mean_vals = cat(3,mean_vals, curr_mean_vals);
    all_combo_probs = log10(all_combo_probs);
    all_stim_probs = log10(all_stim_probs);
    for i = 1:size(combos,1)
        idx = ((i-1)*20 + 1):i*20;
        ylims = [min([all_stim_probs(idx,:),all_combo_probs(idx,:)],[],'all'),...
            max([all_stim_probs(idx,:),all_combo_probs(idx,:)],[],'all')];

        [~,pvals(i,1)] = ttest2(all_stim_probs(idx,1),all_stim_probs(idx,2));

        data1 = max(all_combo_probs(idx,1:2),[],2);
        data2 = max( all_combo_probs(idx,3),[],2);
        diffs(i,loop) = mean(data1) - mean(data2);
        [~,pvals(i,3)] = ttest2(data1, data2);
        [~,pvals(i,2)] = ttest2(all_stim_probs(idx,3),all_stim_probs(idx,4));

        [paired_p(i,1)] = signrank(all_stim_probs(idx,1),all_stim_probs(idx,2), "Tail","left");
        [paired_p(i,2)] = signrank(all_stim_probs(idx,3),all_stim_probs(idx,4), "Tail","left");
        [paired_p(i,3)] = signrank(data1, data2, "Tail","left");

 
    end
    all_ps = cat(3,all_ps, pvals);
    all_paired_ps = cat(3,all_paired_ps , paired_p);
end

FDR_pvals = mafdr(all_paired_ps(:),'BHFDR', true);
FDR_pvals = reshape(FDR_pvals',size(all_paired_ps,1),size(all_paired_ps,2),size(all_paired_ps,3));

all_paired_ps = FDR_pvals;


all_means = cat(3, all_means, sum(all_paired_ps< 0.05,3));
combo_sigs = cat(2,combo_sigs,all_means(:,3));
all_mean_vals = cat(3, all_mean_vals , mean(mean_vals,3));
method = 2;

[combo_pts , y] = PlotClusters(all_paired_ps(:,3,:), all_most_like,method );
all_y(:, neuron(1)+1) = y;
class_count{neuron(1)+1} = CountClasses(combo_pts);


end


function [combo_pts , y] = PlotClusters(ps, most_like,method )
x = [];
y = [];
combo_names = {'DANaClcombo','DAIAAcombo','QuinNaClcombo','SDSDAcombo','NaClSDS','IAANaCl',[],[],[]};
combo_pts  = [];

for i = 1:size(ps,1)
    subplot(2,3,i)
    a = ps(i,:);
    if method == 1
        y(i) = (100 - sum(a<0.05));
    else
        y(i) = (mean(a));
    end
    b = most_like(:,i);
    x(i,:) = [sum(b == 1), sum(b == 2), sum(b == 3), sum(b == 4)];

    x(i,:) = x(i,:)/sum(x(i,:));
    pts = rand(1,100);
    pts(pts < x(i,1)) = 1;
    pts(pts < sum(x(i,1:2))) = 5;
    pts(pts < sum(x(i,1:3))) = 2;
    pts(pts < sum(x(i,1:4))) = 3;
    if x(i,1) < x(i,2)
        pts(pts == 2) = 4;
    end
    xnoise = randn(1,100)/4;
    ynoise = y(i) + randn(1,100)*0.05;
    combo_pts = cat(1, combo_pts, pts);



end
end

function class_count = CountClasses(pts)
class_count(:,1) = sum(pts == 1,2);
class_count(:,2) = sum(pts == 2,2) + sum(pts == 4,2);
class_count(:,3) = sum(pts == 3,2);
class_count(:,4) = sum(pts == 4,2);


end