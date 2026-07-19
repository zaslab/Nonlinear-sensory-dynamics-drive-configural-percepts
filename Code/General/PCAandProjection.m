function [score,combo_proj, combo_ids,coefforth , mu, explained] = PCAandProjection(data,ids,range)

data = NormalizeForPCA(data);
[coefforth,score,~,~,explained,mu]= pca(data(1:range,:));
combo_proj = (data((range+1):end,:)-mu) * coefforth;
combo_ids = ids((range+1):end);
end