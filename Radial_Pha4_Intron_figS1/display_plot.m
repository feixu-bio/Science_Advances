%% plot the 3d reconstruction
function display_plot(spot,T_alpha,emb)
    intronx = [spot.x];
    introny = [spot.y];
    intronz = [spot.z];   
    figure
    hold on
    % plot the alphashape of the nuclei
    for i = 1:height(T_alpha)
        h{1} = plot(T_alpha.alphaShape{i,1} ,'FaceColor','blue','FaceAlpha',0.05,'LineStyle', 'none');
    end
    % plot the intron dots
    h{2} = scatter3(intronx, introny, intronz , '.','red');
    % plot the center of the nuclei
    legend([h{1}, h{2}], {'nuclei', 'pha-4_intron'} );
    title(['Emb_' num2str(emb) 'Reconstruction']);
    hold off
