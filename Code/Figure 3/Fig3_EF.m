warning off
load('combo_list.mat')
max_loop = 100;
combonames = data_cell(combo_list(:,3),3);
num_to_sample = 20;
range = 13*num_to_sample ;

maxPC = 3;
feat_pcs = 1:3;
[all_y, all_paired_ps,diffs,res]= RunComboPvals(data_labels, data, feat_score,...
               feat_pcs,combo_list,range,maxPC, max_loop);


 figure; 
 counter = 1;
 for i = 1:6
    subplot(2,3,i);
    dat = permute(all_paired_ps(i,3,:),[2,3,1]);
    histogram(dat,'BinEdges',[0:0.05:1]);
    ylim([0 100]);
    title(combonames(i))
    xlabel('p-values')
    ylabel('Count')
 end