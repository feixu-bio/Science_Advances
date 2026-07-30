function[radial,attached] = radial_distribution_slurm(SampleNum)

%%
clc

SampleNum

segment_folder = 'segmentation/';
intron_spot_folder = 'spot/';
intron_Im_folder = 'intron_im/';
LADs_thresh = 0.2; % in um if the distance to the boundary is less 0.20 um.
% Then consider as association.
sz = 2304;
scale3D = [0.034 , 0.034 , 0.13]; % pixel size
mkdir savefiles1

load('DeltaZ.mat');
load('tform.mat');

filelist = dir([segment_folder '*.mat']);
tic
% main call
[Age,T_alpha,emb,TransImg] = Generate_NucTable(SampleNum,segment_folder,sz,tform,scale3D);
disp('segmentation_done');
[spot,background_intensity] = reshape_spot(emb,intron_spot_folder,intron_Im_folder,DeltaZ,...
    TransImg,sz);
disp('table_reshape_done')
[radial] = distance_boundary(T_alpha,spot,background_intensity);
[attached] = attach_boundary(radial,LADs_thresh);
disp('calculation_done')
display_plot(spot,T_alpha,emb)
disp('saving_results')
saveas(gcf,fullfile(['savefiles1/Emb_' num2str(emb) 'Reconstruction.fig']));
save(['savefiles1/radial_intensity_' num2str(emb) '.mat'],'radial','attached','Age');
save(['savefiles1/alphashape_' num2str(emb) '.mat'],'T_alpha');
toc
end
