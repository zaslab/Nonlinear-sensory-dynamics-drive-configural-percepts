function aligned_trace = CheckIfASH(data, max_change,  varargin)
window_end = 138;
for i=1:1:length(varargin)
   
    if strcmp(varargin{i},'window_end') 
        window_end = varargin{i+1};
    end
end

[peak_size, peak] = max(data(40:60));
start = mean(data(1:10));
post_step = mean(data((max_change + 40):(max_change + 50)));
if  peak_size > (start + 0.08) & post_step < (start - 0.1)
    max_change = max([peak - 8,1]);

end
aligned_trace = data((max_change) : (max_change + window_end - 1));
end
