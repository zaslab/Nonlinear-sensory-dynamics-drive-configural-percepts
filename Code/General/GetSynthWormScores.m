function [all_score, raw_acts, acts, ids,feature_worms,coefforth,mu,explained,combo_proj,time_traces] = GetSynthWormScores(data_labels, data, feat_score,range,varargin)
feat_pcs = 1:3;
sample_num = 20;
for i=1:1:length(varargin)
    if strcmp(varargin{i},'feat_pcs')
        feat_pcs = varargin{i+1};
    end
    if strcmp(varargin{i},'sample_num')
        sample_num = varargin{i+1};
    end
    if strcmp(varargin{i},'time_traces')
        trace_data = varargin{i+1};
    end
end

[acts, ids,feature_worms] = ...
    MakeSynthWorms(data_labels, data, feat_score, sample_num, feat_pcs);

% [time_traces, ids,feature_worms] = ...
    % MakeSynthWormTimeTraces(data_labels, trace_data, feat_score, sample_num, feat_pcs);
% [acts] = ...
%     TestMakeSynthWorms(data_labels, data, feat_score, sample_num, feat_pcs,'completion_method', 'IALM');

mean_data = nanmean(acts, 1);
std_data = nanstd(acts, [], 1);
raw_acts = acts;
for i = 1:size(acts,1)
    acts(i, :) = (acts(i, :) - mean_data) ./ std_data;
end
[coefforth,score,~,~,explained,mu]= pca(acts(1:range,:));
combo_proj = (acts((range+1):end,:)-mu) * coefforth;
all_score = [score; combo_proj];
end
