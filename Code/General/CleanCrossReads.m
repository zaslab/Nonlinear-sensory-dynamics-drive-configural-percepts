function [data] = CleanCrossReads(data,data_labels,cross_reads)




for i = 1:size(cross_reads,1)
    ref_cond = cross_reads(i,4);
    ref_data = data(data_labels(:,2) == ref_cond,:);
    ref_labels = data_labels(data_labels(:,2) == ref_cond,:);
    neuron_num = cross_reads(i,1);
    cond_num = cross_reads(i,2);
    if cross_reads(i,3) == 3
        step = [1,2];
    else
        step = cross_reads(i,3) ;
    end

    for step_num = 1:length(step)
        curr_step = step(step_num);
        ref_traces = ref_data(ref_labels(:,1) == neuron_num& ...
                              ref_labels(:,3) == curr_step,:);
        replace_locs = find(data_labels(:,1) == neuron_num &...
                            data_labels(:,2) == cond_num &...
                            data_labels(:,3) == curr_step);
        
        data(replace_locs,:) = ref_traces(randi(size(ref_traces,1),1,length(replace_locs)),:);
    end
end

end