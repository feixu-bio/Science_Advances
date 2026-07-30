clear all 
close all
tic

% mkdir sequential_cellpose
sz = 500;

%% primary probe file
filelist = dir('cellpose_terri/*.tif');
for ii = 1:length(filelist)
    data = bfopen(['cellpose_terri/' filelist(ii).name]);
    T = extractBetween(filelist(ii).name,'RNA-DNA-','Terr_');
    emb = T{1,1};
    ImageStack = zeros(sz,sz,151);
    for i = 1:151 %planes are organized by z position first then channel
        ImageStack(:,:,i) = data{1, 1}{i,1};
    end
    

    %channel 1 = 647
    ImageStack647 = zeros(sz,sz,151);
    ImageStack647 = ImageStack(:,:,1:1:151);
    ImageMax = max(ImageStack647,[],3);
    figure
    imagesc(ImageMax)
    axis equal
    colormap gray
    title(['DNA_647_emb_', emb])
    saveas(gcf,['sequential_cellpose/DNA-647-emb-' emb])
    save(['sequential_cellpose/DNA-647-emb-' emb '.mat'],'ImageStack647', '-v7.3');
    clear ImageMax ImageStack647

    close all  
    clear data
end
toc
