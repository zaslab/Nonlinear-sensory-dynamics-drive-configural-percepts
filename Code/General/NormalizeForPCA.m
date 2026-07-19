function data = NormalizeForPCA(data)

mean_data = nanmean(data, 1);
std_data = nanstd(data, [], 1);
for i = 1:size(data,1)
    data(i, :) = (data(i, :) - mean_data) ./ std_data;
end
end