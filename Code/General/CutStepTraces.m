function [data_for_pca, variable_matrix] = CutStepTraces(data)
end_idx = 140;
variable_matrix = [];
data_for_pca = [];
counter = 1;

current_data = data;
step_window = [40, 119];

for cond = 1:size(data,1)
    current_cond = current_data{cond,1};
    current_steps = current_data{cond, 4};
    for worm = 1:size(current_cond,3)
        current_worm = current_cond(:,:,worm);
        for step = 1:3
            earliest_response = current_steps(:,step,worm);
            for neuron = 1:22
                % clean miss-marked step starts to avrg of all other
                % neurons
                [earliest_response, current_worm(neuron,:)] = CheckStepValidity(earliest_response,neuron, current_worm(neuron,:),step_window);
                current_step(neuron,:) = current_worm(neuron, ...
                    (earliest_response(neuron) - step_window(1)):(earliest_response(neuron) + step_window(2)));
                current_step_sub(neuron,:) = current_step(neuron,1:(end_idx));
            end
            data_for_pca(counter:(counter + 21),:) = current_step_sub;

            variable_matrix(counter:(counter + 21),:) = cat(2,(1:22)', repmat([cond, step],22,1));

            counter = counter + 22;
        end
    end
end

end

function fixed_step_starts = clean_step_starts(all_step_starts, neuron_num)
fixed_step_starts = all_step_starts;
all_step_starts(neuron_num) = [];
fixed_step_starts(neuron_num) = round(nanmean(all_step_starts));
end

function [earliest_response, current_worm] = CheckStepValidity(earliest_response,neuron, current_worm,step_window)
if earliest_response(neuron) > (length(current_worm) - step_window(2))
    earliest_response = clean_step_starts(earliest_response, neuron);
end

if earliest_response(neuron) > (length(current_worm) - step_window(2)) | isnan(earliest_response(1))
    earliest_response(:) = 100;
    current_worm = nan;
end
if (earliest_response(neuron) - step_window(1)) < 1
    earliest_response(neuron) = 1 + earliest_response(neuron) + (abs(earliest_response(neuron) - step_window(1)));
end
end

