function acts_by_combo = GenerateComboActs(y, sample_num)
combo_list = [2 4 8
    2 1 9
    4 5 10
    2 6 11
    4 6 12
    1 4 13
    ];

%     [all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained,combo_proj] = ...
%         GetSynthWormScores(data_labels, data, feat_score,range,'sample_num', sample_num);
% 
% y = raw_acts;


acts_by_combo = [];

for n = 1:26
    neuron_data = [];
    for i = 1:size(combo_list,1)
        c = combo_list(i,:);
        idx1 = ((c(1)-1)*sample_num + 1):c(1)*sample_num;

        idx2 = ((c(2)-1)*sample_num + 1):c(2)*sample_num;

        idx3 = ((c(3)-1)*sample_num + 1):c(3)*sample_num;
        
        s1 = y(idx1,n);
        s2 = y(idx2,n);
        s3 = y(idx3,n);
        if length(c) > 3
            idx4 = ((c(4)-1)*sample_num + 1):c(4)*sample_num;
            s4 = y(idx4,n);
        else
            s4 = [];
        end

        plotdata = [s1,s2,s3,s4];
       

        meandata = [mean(plotdata,1), n,i];
        [~,idx] = sort(abs(meandata([1,2])));
        % [~,idx] = sort((meandata([1,2])));
        meandata([1,2]) = meandata(idx);
        neuron_data = cat(1, neuron_data, meandata);
   
    end
 
    acts_by_combo = cat(1, acts_by_combo, neuron_data );

end
%%
% figure; 
% pts = acts_by_combo(3:end,:);
% scatter3(pts(:,1), pts(:,2), pts(:,3), 36, pts(:,4), 'filled');