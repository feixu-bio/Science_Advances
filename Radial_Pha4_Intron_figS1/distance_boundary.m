function [spot,background_intensity] = reshape_spot(emb,intron_spot_folder,intron_Im_folder,DeltaZ,...
    TransImg,sz)
%load csv
    if emb < 10
        T = readmatrix([intron_spot_folder '20230401pha4smFISHintron-experiment2_0' num2str(emb) '_intron.csv']);
    elseif emb >= 10
        T = readmatrix([intron_spot_folder '20230401pha4smFISHintron-experiment2_' num2str(emb) '_intron.csv']);
    end
    disp('csv loaded')
    %remove NaN
    T(:,[1,3,12])=[];
    T(1,:)=[];
    % reassign matrix into structure
    % 3d position of each ID
    disp('Nan removed')
    %culculate the median intensity in intron channel within the nuclei
    if emb < 10
        load([intron_Im_folder 'segmentation-' num2str(emb) '.mat']);
    elseif emb >= 10
        load([intron_Im_folder 'segmentation-' num2str(emb) '.mat']);
    end
    Mask = TransImg > 0;
    background_intensity = median(ImageStackseg(Mask));
    disp('BGintensity done')
    k = 0 ;
    for i = 1:height(T)
        if T(i,17)>1.0 && T(i,16)>0.2
            %This is based on the examination of individual spots fitting.
            %And ideal spots fitting will be SNR over 1 and contrast over
            %0.2
            k=k+1;
            spot(k).ID = i;
            spot(k).QUALITY = T(i,2);
            spot(k).x = T(i,3); % um
            spot(k).y = sz*0.0338542 - T(i,4); % need to consider the flip in um
            %align z_560 to z_405. deltaz = z_560-z_405; so z_560 =
            %deltaz-z_405;
            spot(k).z = T(i,5)-DeltaZ*0.13; % um
            spot(k).mean_inty = T(i,10);
            spot(k).median_inty = T(i,11);
            spot(k).max_inty = T(i,13);
            spot(k).std_inty = T(i,15);
            spot(k).contrast = T(i,16);
            spot(k).SNR = T(i,17);
        end
    end
    disp('forloop done')
