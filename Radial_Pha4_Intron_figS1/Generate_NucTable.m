%% Generate alphashape of each nuclei in embyro
function[Age,T_alpha,emb,TransImg] = Generate_NucTable(SampleNum,segment_folder,sz,tform,scale3D)
    alpha_radius = sqrt(sum(scale3D.^2)); 
    % Use a maximum alpha of the longest possible diagonal
    filelist = dir([segment_folder '*.mat']);
    load([segment_folder filelist(SampleNum).name]);
    emb = str2double(extractBetween(filelist(SampleNum).name,'segmentation-','.mat'));
    % here imageStackseg is labeled images
    slices = size(ImageStackseg,3);
    %align DAPI to intron
    TransImg = zeros(sz,sz,slices);
    for m = 1:slices
        TransImg(:,:,m) = imtransform(ImageStackseg(:,:,m), tform, 'XData', [1 sz], 'Ydata', [1 sz]);
    end
    CC2 = regionprops3(TransImg,'Volume');
    V = [CC2.Volume];
    %     C = 1:height(CC2);
    %     nucID = C.';
    %     T = addvars(CC2,nucID,'Before',"Volume");
    %     nuclist = table2array(CC2);
    % filter out small volume
    % V = [T.Volume];
    % Here xy resolution is 0.033um/pix z resolution is 0.13um/pix
    % take diameter of 2um minimum as example. Since it is anisotrophic, so
    % I assume the nuclei as a cubic volume, xyz 60*60*15 voxel.
    % in total 54000.
    %     small_ind = find(V<40000);
    %     nuclist(small_ind,:) = [];
    %     L = 1:height(nuclist);
    %     nucID = L.';
    %     nucTable = array2table(nuclist,...
    %     'VariableNames',{'Volume'}); %get rid of small objects
    %     nucTable = addvars(nucTable,nucID,'Before',"Volume");
    indx = find(V>=40000);% consider size over 54000 as true nuclei. 
    % the label in raw seg image is stored in indx
    Age = length(indx); % embryonic stage
    % need to consider the preallocation.   
    construct = cell(Age,3);
    for label = 1 : Age
        single_Mask = TransImg == indx(label);
        [Y, X, Z] = ind2sub(size(single_Mask),find(single_Mask));
        Y_new = sz - Y; 
        %convert to um into isotrophic
        X_um = X*0.0338542; % um
        Y_um = Y_new*0.0338542; % um
        Z_um = Z*0.13; % um
        shp = alphaShape(X_um, Y_um, Z_um,alpha_radius); % alphashape in um
        v = volume(shp);
        %create the boundary in um. Using boundary instead of alphashape is
        %sometimes there is a shrinkage of the shape.
%         P = [X_um Y_um Z_um]; % set of points
%         [~, v] = boundary(P,0); % No shrinkage
        % calculate the volume in um
        radius = ((0.75*v)/3.1416)^(1/3); %radius of the volume in um
        construct(label,:) = {label,shp,radius};
    end
    T_alpha = cell2table(construct,"VariableNames",["label" "alphaShape" "radius"]);
