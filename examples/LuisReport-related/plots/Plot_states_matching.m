index = reshape(1:9, 3,3)';
figure
for i = 1:9
    subplot(3,3,index(i))
    heatmap(sum(results(i).HMM.true_to_inferred_mapping,2));
    colorbar('off');
    title('HMM')
end
sgtitle('States matching HMM (true x modelled)')

for i = 1:9
    subplot(3,3,index(i))
    heatmap(sum(results(i).HSMM.true_to_inferred_mapping,2));
    colorbar('off');
    title('HSMM')
end
sgtitle('States matching HSMM (true x modelled)')