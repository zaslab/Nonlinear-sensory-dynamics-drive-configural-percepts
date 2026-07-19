function [active_variable_matrix,unique_ids] = CombineNeuronSides(active_variable_matrix,unique_ids,combine_all)

if combine_all
neuron_pairs = [6,15;
    1,16;
    2,20;
    7,12;
    8,11;
    3,18;
    5,13;
    10,17;
    14,22;
    4,19;
    9,21];
else
    neuron_pairs = [6,15;
    16,nan;
    1,nan;
    2,20;
    7,12;
    8,11;
    3,18;
    5,13;
    10,17;
    14,22;
    4,19;
    9,nan;
    21,nan];
end
for neuron = 1:size(neuron_pairs,1)

    if isnan(neuron_pairs(neuron,2))
        continue

    else
        active_variable_matrix(active_variable_matrix(:,1) == neuron_pairs(neuron,2)) = neuron_pairs(neuron,1);
        unique_ids(unique_ids(:,1) == neuron_pairs(neuron,2),:) = [];
    end


end
end
