%% Load and clean data
range = 140;
data_cell = load('data_cell.mat').data_cell;
load('relevant_neurons.mat')
%% Fixes

data_cell{6,1}(6,:,12) = nan;

%%
% ; 3,1,3,7
% 'cross_reads', [10,3,3,7; 3,2,3,7],
[data,data_labels,feat_score,raw_features,data_cell,time_series] = ...
      PrepData(data_cell,range,'cross_reads', [10,3,3,7; 3,2,3,7;3,1,3,7],'trace_features');
%%
% rm_neurons = [7,8];
% [data, data_labels] = remove_neurons(data, data_labels, rm_neurons);
combo_list = [2 4 8
    2 1 9
    4 5 10
    2 6 11
    4 6 12
    1 4 13
    ];

range = 140;
[all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained1,combo_proj] = ...
                   GetSynthWormScores(data_labels, data, feat_score, range, 'time_traces', time_series);
% figure; imagesc(raw_acts)

%%

% num_to_sample = 20;
% [new_data, new_ids] = CompleteActData(data_labels, data,num_to_sample);
% [new_features, ~] = CompleteActData(data_labels, data,num_to_sample);
% [all_score,raw_acts,acts,ids,feature_worms,coefforth,mu,explained,combo_proj] = ...
%                    GetSynthWormScores(new_ids, new_data, new_features,range,'feat_pcs',1);
% % figure; imagesc(raw_acts)
% 
% [data,data_labels,~,~,data_cell] = ...
%       PrepData(data_cell,range);