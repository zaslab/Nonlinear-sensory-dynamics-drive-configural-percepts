function [wt_data,data_labels,feat_score,raw_features,data_cell,data] = PrepData(data_cell, range,varargin)
get_features = 0;
combine_all = 0;
cross_reads = 0;
feat_score = [];
raw_features = [];
for i=1:1:length(varargin)
    if strcmp(varargin{i},'remove_conds')
        data_cell([varargin{i+1}],:) = [];
    end
    if strcmp(varargin{i},'cross_reads')
        cross_reads = varargin{i+1};
    end
    if strcmp(varargin{i},'trace_features')
        get_features = 1;
    end
     if strcmp(varargin{i},'combine_all')
        combine_all = 1;
    end
    
end

load('relevant_neurons.mat')
dbstop if error
data_cell = SortAWCs(data_cell);
[all_neurons_data, variable_matrix] = CutStepTraces(data_cell);

for u = 1:size(all_neurons_data,1)
    all_neurons_data(u,:) = smooth(all_neurons_data(u,:));
end

data = Reallign(all_neurons_data, variable_matrix,'window_end', 120);
[data, data_labels] = GetSecondStepForLightResponders(data, variable_matrix);
% data_labels = variable_matrix;
data(data_labels(:,3) == 3,:) = [];
data_labels(data_labels(:,3) == 3,:) = [];
[unique_ids] = unique(data_labels, 'rows', 'first');
[data_labels] = CombineNeuronSides(data_labels,unique_ids,combine_all);
if cross_reads
    data = CleanCrossReads(data,data_labels,cross_reads);
end
if get_features
    [raw_features,feat_score] = GetTraceFeatures(data, range);
end

[~,wt_data] = GetActivationMagnitude(data_labels, data,0, ...
    'all_window', 36:60, 'ASH_window', 36:50);

end

function cleaned_data = SortAWCs(data)
% Sort AWCs by strength

for i = 1:size(data,1)
    if ismember(i,[4,14])
        step = 1;
    else
        step = 2;
    end
    data(i,:) = sortByStrength(data(i,:), 1,{'AWC'},step);
end
cleaned_data = data;
end

function [all_trace_features,all_feat_score] = GetTraceFeatures(data, range)
currpath = cd;
% cd('D:\Google Drive\code\catch22-main\wrap_Matlab')
cd('G:\My Drive\code\catch22-main\wrap_Matlab')
all_trace_features = zeros(size(data,1), 22);
for i = 1:size(data,1)
    all_trace_features(i,:) = catch22_all(smooth(data(i,:)), false);
end

[coefforth,feat_score,~,~,~,mu]= pca(all_trace_features(1:range,:));
feat_combo_proj = (all_trace_features((range+1):end,:)-mu) * coefforth;
all_feat_score = [feat_score; feat_combo_proj];
cd(currpath)
end