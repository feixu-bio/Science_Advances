%%
% This script is to calculate the distance between pha-4 region and the edge
% of central territories.

clc
clear all
close all
imtool close all

global SampleNum
global RNASampleNum

%% setttings
SampleNum = '6';
RNASampleNum = '6';
% -------------------------------------------------------------------------
% if you SAVE the aligned Tiff file for both RNA and DNA. Then this not need 
% T = readtable('RNA-DNA.txt');
% xStartPoint = T.xStartPoint(str2num(RNASampleNum));
% yStartPoint = T.yStartPoint(str2num(RNASampleNum));
% FOV = T.FOV(str2num(RNASampleNum));
%--------------------------------------------------------------------------
xStartPoint = 0;
yStartPoint = 0;
FOV = 1;

Set.sz = 500; % the length of the ROI
ROI.x1 = xStartPoint + 1;
ROI.x2 = xStartPoint + Set.sz;
ROI.y1 = yStartPoint + 1;
ROI.y3 = yStartPoint + Set.sz;

% set boundary for remove RNA outside embryo generated from alignment 
Set.rnax1 = 50;
Set.rnax2 = 475;
Set.rnay1 = 50;
Set.rnay2 = 475;

Set.remove_small_object = 1;
Set.radius = 5; %pixel

%for reconstruction
Set.dilateThresh = 5; % 5
Set.erodeThresh = 5; % 5

% do you want to filter out territories in the nuclei over the number of 4
% over segmentation
filter_over4terri = 0;

%for binarization and watershed strength%%%change terrmin if traces are too
%many/few per territory
% Set.TerrMin = 0.25; %0.25 lower is more seperated
Set.NucMin = 0.25; %0.4
% Set.TerrDia = 1;
Set.dot_Nuc_thresh = 0.5; %distance of dot to terri is 0.5 um
Set.mRNA_dis_thresh = 1; %distance thresh to nuclei is 1 um
Set.mRNA_Num_thresh = 10; %minimum number of mRNA dots
Set.nRNA_dis_thresh = 0.5; %distanct of nRNA to nuc is 0.5 um
%% main call
if ~(isfolder('3D_Reconstruction'))
    mkdir('3D_Reconstruction')
end
FigurePath = '3D_Reconstruction'; %relative path to save figure
load(strcat('Sequential/RNA_405_emb_', RNASampleNum, '.mat'));
%DAPI in DNA FISH is not good so use RNA FISH signal
load(['sequential_cellpose/DNA-647-emb-' SampleNum]); 
load('tform.mat');
load(['DeltaZ.mat']);
load(['RNA_fit_result_' num2str(RNASampleNum) '.mat']);
load(['spotdetection_' num2str(SampleNum) '.mat']);
[nucLall,age,Set] = seg_territory(ImageStack405,405,tform,ROI,Set);
% [L2,NumTerr,Set] = seg_territory(ImageStack647,647,tform,ROI,Set);
L2 = ImageStack647;
nucL =  excludeMitotic(nucLall, ImageStack405,FigurePath,SampleNum);
clear ImageStack405 ImageStack560 ImageStack647 
foci = reshapeDNAfoci(spot.Xfit,spot.Yfit,spot.Zfit, ROI,Set); %reshape region 28 dot 
[mRNA, nRNA] = qualityRNA(ListOfmRNA,ListOfnascentRNA,Set);
Nuc = maskshape(nucL);%create nuc mask 
Terr2Struct = maskshapeDeltaZ(L2,DeltaZ); %create terr2 mask plus chromatic abbration
NucStruct = DoQualityControl(Nuc,405); %do quality control 
% Terr2Struct = DoQualityControl(Terr2,647); %do quality control 
Dot = DotQualityControl(foci,NucStruct); %do quality control of dots only in the L1 or less than 1 um 
[pha4Nuc, nonpha4Nuc] = classify_mRNA_nuc(mRNA, NucStruct,Set); %classify nuc
List_Nuc_Intron = classify_intron_nuc(nRNA, NucStruct,Set.nRNA_dis_thresh);%associate intron with nuclei
List_Nuc_Terr2 = NucTerrAssociation(NucStruct, Terr2Struct,filter_over4terri);%assign Terr2 with pha4 nuc
List_Nuc_Dot = NucDot(NucStruct, Dot, Set);
List_Nuc_Dot_Terri = concatlist(pha4Nuc,List_Nuc_Dot,List_Nuc_Terr2,List_Nuc_Intron);
[pha_dis, nonpha_dis, BothON_dis, OneON_dis, BothOFF_dis] = calculateDis(List_Nuc_Dot_Terri, Dot, NucStruct);
display(mRNA,nRNA,pha4Nuc,nonpha4Nuc,Dot,List_Nuc_Dot_Terri,SampleNum,FOV,FigurePath,List_Nuc_Intron)
%% save data
save(['FOV-' num2str(SampleNum) '-Emb-' num2str(FOV) '-Distance to permissive domain.mat'],...
    'pha_dis', 'nonpha_dis','BothON_dis', 'OneON_dis', 'BothOFF_dis','age');
save(['FOV-' num2str(SampleNum) '-Emb-' num2str(FOV) '-Reconstruction.mat'],...
    'Set','Dot','mRNA','nRNA','NucStruct','Terr2Struct','pha4Nuc','nonpha4Nuc','List_Nuc_Dot_Terri');
%% main function
% segmention
function [Terri_seg,numvol,Set] = seg_territory(ImageStack,channel,tform,ROI,Set)
switch channel
    case 405
        %         I = ImageStack(ceil(ROI.y1):ceil(ROI.y3),ceil(ROI.x1):ceil(ROI.x2),:);
        A = ImageStack;
        I = imgaussfilt3(A,2);
        %correct for uneven background illumination
        se = strel('disk',50);
        background = imopen(I,se);
        I2 = I - background;
        Im = max(I2,[],3);
        imtool(Im)
        disp("please select tools and adjust constrst")
        disp("sliding the maximum to where region inside is saturated")
        disp("press anykey to continue")
        pause
        prompt = {'Enter minThresh:','Enter maxThresh:'};
        dlgtitle = 'Input';
        answer = inputdlg(prompt,dlgtitle);
        Set.minThresh405 = str2num(answer{1,1});
        Set.maxThresh405 = str2num(answer{2,1});
        Set.TerrbinaryThresh405 = Set.minThresh405*1.25
        %set min and max values to even out territory intensity
        I2(I2<Set.minThresh405) = 0;
        I2(I2>Set.maxThresh405) = Set.maxThresh405;
        %3D binary
        BW = imbinarize(I2, Set.TerrbinaryThresh405);
        %remove small object
        if Set.remove_small_object == 1
        p = round(1.34*3.14*10^3);
        BW = bwareaopen(BW,p,26);
        end
        % fillin the small holes
        BW2 = imfill(BW, 'holes');
        figure
        imshowpair(max(I2,[],3),max(BW2,[],3),'montage')
        title('before and after 3D binary')
        %3D dist transform
        D = -bwdist(~BW2); 
        %set background to its own catchment basin - i.e. where BW is black
        D(~BW2) = -Inf;
        %3D watershed
        D = imhmin(D,Set.NucMin); %the height threshold for suppressing shallow minima
        L = watershed(D);
        L = imdilate(L,true(5));
        Terri = L-1;
        num = unique(Terri(Terri>0));
        Lmax = max(Terri,[],3);
        Lrgb = label2rgb(Lmax,'jet','w','shuffle');
        figure
        imagesc(Im)
        colormap gray
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.3;
        title('volume labels superimposed transparently on raw image')
        % question if you want to restart segmentation
        answer = questdlg('Are you satisfied with segmentation?','Yes','No');
        if strcmp(answer,'No')
           error('change the threshold for segmentation.')    
           return
        end
        Terri_seg = Terri;
        numvol = num;
    case 560
        I = ImageStack(ceil(ROI.y1):ceil(ROI.y3),ceil(ROI.x1):ceil(ROI.x2),:);
        %correct for uneven background illumination
        se = strel('disk',10);
        background = imopen(I,se);
        I2 = I - background;
        Im = max(I2,[],3);
        imtool(Im)
        disp("please select tools and adjust constrst")
        disp("sliding the maximum to where region inside is saturated")
        disp("press anykey to continue")
        pause
        prompt = {'Enter minThresh:','Enter maxThresh:'};
        dlgtitle = 'Input';
        answer = inputdlg(prompt,dlgtitle);
        Set.minThresh560 = str2num(answer{1,1});
        Set.maxThresh560 = str2num(answer{2,1});
        Set.TerrbinaryThresh560 = Set.minThresh560*1.25
        %set min and max values to even out territory intensity
        I2(I2<Set.minThresh560) = 0;
        I2(I2>Set.maxThresh560) = Set.maxThresh560;
        %3D binary
        BW = imbinarize(I2, Set.TerrbinaryThresh560);
        %remove small object
        if Set.remove_small_object == 1
        p = round(1.34*3.14*5^3);
        BW = bwareaopen(BW,p,26);
        end
        % fillin the small holes
        BW2 = imfill(BW, 'holes');
        figure
        imshowpair(max(I2,[],3),max(BW2,[],3),'montage')
        title('before and after 3D binary')
        %3D dist transform
        D = -bwdist(~BW2); 
        %set background to its own catchment basin - i.e. where BW is black
        D(~BW2) = -Inf;
        %3D watershed
        D = imhmin(D,Set.TerrMin); %the height threshold for suppressing shallow minima
        L = watershed(D);
        L = imdilate(L,Set.TerrDia);
        Terri = L-1;
        num = unique(Terri(Terri>0));
        L(find(~Terri))=0;
        Lmax = max(Terri,[],3);
        Lrgb = label2rgb(Lmax,'jet','w','shuffle');
        figure
        imagesc(Im)
        colormap gray
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.3;
        title('volume labels superimposed transparently on raw image')
        % question if you want to restart segmentation
        answer = questdlg('Are you satisfied with segmentation?','Yes','No');
        if strcmp(answer,'No')
           error('change the threshold for segmentation.')    
           return
        end
        Terri_seg = Terri;
        numvol = num;
    case 647        
        for m = 1:size(ImageStack,3)
            ImageStack(:,:,m) = imtransform(ImageStack(:,:,m), tform, 'XData', [1 500], 'Ydata', [1 500]);
        end
        I = ImageStack(ceil(ROI.y1):ceil(ROI.y3),ceil(ROI.x1):ceil(ROI.x2),:);
        %correct for uneven background illumination
        se = strel('disk',10);
        background = imopen(I,se);
        I2 = I - background;
        Im = max(I2,[],3);
        imtool(Im)
        disp("please select tools and adjust constrst")
        disp("sliding the maximum to where region inside is saturated")
        disp("press anykey to continue")
        pause
        prompt = {'Enter minThresh:','Enter maxThresh:'};
        dlgtitle = 'Input';
        answer = inputdlg(prompt,dlgtitle);
        Set.minThresh647 = str2num(answer{1,1});
        Set.maxThresh647 = str2num(answer{2,1});
        %set min and max values to even out territory intensity
        I2(I2<Set.minThresh647) = 0;
        I2(I2>Set.maxThresh647) = Set.maxThresh647;
        Set.TerrbinaryThresh647 = Set.minThresh647*1.25
        %3D binary
        BW = imbinarize(I2, Set.TerrbinaryThresh647);
        %remove small object
        if Set.remove_small_object == 1
        p = round(1.34*3.14*5^3);
        BW = bwareaopen(BW,p,26);
        end
        % fillin the small holes
        BW2 = imfill(BW, 'holes');
        figure
        imshowpair(max(I2,[],3),max(BW2,[],3),'montage')
        title('before and after 3D binary')
        %3D dist transform
        D = -bwdist(~BW2); 
        %set background to its own catchment basin - i.e. where BW is black
        D(~BW2) = -Inf;
        %3D watershed
        D = imhmin(D,Set.TerrMin); %the height threshold for suppressing shallow minima
        L = watershed(D);
        L = imdilate(L,Set.TerrDia);
        Terri = L-1;
        num = unique(Terri(Terri>0));
        L(find(~Terri))=0;
        Lmax = max(Terri,[],3);
        Lrgb = label2rgb(Lmax,'jet','w','shuffle');
        figure
        imagesc(Im)
        colormap gray
        hold on
        himage = imshow(Lrgb);
        himage.AlphaData = 0.3;
        title('volume labels superimposed transparently on raw image')
        % question if you want to restart segmentation
        answer = questdlg('Are you satisfied with segmentation?','Yes','No');
        if strcmp(answer,'No')
           error('change the threshold for segmentation.')    
           return
        end
        Terri_seg = Terri;
        numvol = num;
end
end
% exclude mitotic nuclei
function nucL_small_volume =  excludeMitotic(nucLall, ImageStack405,FigurePath,SampleNum)
nucLID = unique(nucLall(nucLall>0)); %fetch all labeled nuclei
nucL = nucLall;
% select the nuclei with small volume and measure
ImDAPI = max(ImageStack405,[],3);
Nucmaxall = max(nucLall,[],3);
nucLrgb = label2rgb(Nucmaxall,'jet','w','shuffle');
figure
imagesc(ImDAPI)
colormap gray
axis equal
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
title('pleasae select mitotic nuclei on DAPI image')
saveas(gcf,fullfile(FigurePath,['FOV-' num2str(SampleNum) '-WithMitotic']));
CC2 = regionprops3(nucL,'all'); %get all volumes of the nuclei.
[B I] = sort(CC2.Volume); % arrage the volume from small to large.
NucVolume = [I B]; % the labels of the nuclei and their volumes
ROI_selection = questdlg('select mitotic nuclei?');
if(strcmp(ROI_selection, 'Yes'))
    roiList.rect = [];
    roiList.index = [];
    roi.rect = [];
    roi.index = [];
    h = [];
    x1 = [];
    x2 = [];
    x3 = [];
    x4 = [];
    y1 = [];
    y2 = [];
    y3 = [];
    y4 = [];
    xi = [];
    yi = [];
    moreROI = [];
    nroi = [];
    selectROI % using selectROI fuction to select mitotic nuclei
    ID_SmallNuclei = [];
        for i = 1:length(roiList)
            % get the positions from roiList
            x1 = round(roiList(i).rect(1));
            x2 = round(roiList(i).rect(1)+roiList(i).rect(3));
            y1 = round(roiList(i).rect(2));
            y2 = round(roiList(i).rect(2)+roiList(i).rect(4));
            % cropout the individual nucleus, make sure the crop region a little bit larger
            % than the nucleus itself to include all the region it has.
            Label_Image_Crop = nucL(ceil(y1):ceil(y2),ceil(x1):ceil(x2),:);
            Nuc_Crop_max = max(Label_Image_Crop,[],3);
            nucLrgb = label2rgb(Nuc_Crop_max,'jet','w','shuffle');
            ImDAPI_crop = ImDAPI(ceil(y1):ceil(y2),ceil(x1):ceil(x2),:);
            figure
            imagesc(ImDAPI_crop)
            colormap gray
            hold on
            txt = ['No.' num2str(i) ' crop nuclei DAPI are shown' ]; 
            title([txt]) %For accuracy, it is better to compare with the previous superimposed image in figure 2
            CC = regionprops3(Label_Image_Crop,'all'); %get the volume of selected nuclei.
            Volume_of_selected_nucleus = max(CC.Volume); %find the max number to exclude some small volume by uncorrect selection
            nucLID = unique(nucL(nucL>0));
            nucL_small_volume = nucL;
            for j = 1:length(nucLID)
                if NucVolume(j,2) == Volume_of_selected_nucleus %find the selected nuclei
                   ID_SmallNuclei = [ID_SmallNuclei NucVolume(j,1)];
                end
            end
            for j = 1:length(ID_SmallNuclei)
            nucL_small_volume(nucL_small_volume == ID_SmallNuclei(j)) = 0; %exclude the prophase labeled nuclei
            end
        end
else
    nucL_small_volume = nucL;
    ID_SmallNuclei = [];
end    
Nucmax = max(nucL_small_volume,[],3);
nucLrgb = label2rgb(Nucmax,'jet','w','shuffle');
figure
imagesc(ImDAPI)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.3;
title('nuclei exclude all mitotic labels superimposed transparently on DAPI')
saveas(gcf,fullfile(FigurePath,['FOV-' num2str(SampleNum) '-ExcludeMitotic']));
end
% Rename the RNAs
function [mRNA, nRNA] = qualityRNA(ListOfmRNA,ListOfnascentRNA,Set)
mRNA.x = ListOfmRNA.XfitMRNA;
mRNA.y = ListOfmRNA.YfitMRNA;
mRNA.z = ListOfmRNA.ZfitMRNA;
nRNA.x = ListOfnascentRNA.XfitNRNA;
nRNA.y = ListOfnascentRNA.YfitNRNA;
nRNA.z = ListOfnascentRNA.ZfitNRNA;

k = 0;
for i = 1:length(mRNA.x) 
    if mRNA.x(i)>Set.rnax1 && mRNA.x(i)<Set.rnax2 && mRNA.y(i)>Set.rnay1 && mRNA.y(i)<Set.rnay2
       k = k+1;
       mRNA_new(k).x = mRNA.x(i);
       mRNA_new(k).y = mRNA.y(i);
       mRNA_new(k).z = mRNA.z(i);
    end
end

mRNA.x = [mRNA_new.x];
mRNA.y = [mRNA_new.y];
mRNA.z = [mRNA_new.z];

k = 0;
for i = 1:length(nRNA.x) 
    if nRNA.x(i)>Set.rnax1 && nRNA.x(i)<Set.rnax2 && nRNA.y(i)>Set.rnay1 && nRNA.y(i)<Set.rnay2
       k = k+1;
       nRNA_new(k).x = nRNA.x(i);
       nRNA_new(k).y = nRNA.y(i);
       nRNA_new(k).z = nRNA.z(i);
    end
end

nRNA.x = [nRNA_new.x];
nRNA.y = [nRNA_new.y];
nRNA.z = [nRNA_new.z];
end
% reshape DNA focifitting from um to pixel
function Dot = reshapeDNAfoci(Xfit,Yfit,Zfit,ROI,Set)

%rescale the dot
x = [];
y = [];
z = [];
for i = 1:length(Xfit) %get the region28
    x = [x Xfit(1,i)/0.11];
    y = [y (Set.sz*110/1000-Yfit(1,i))/0.11];
    z = [z Zfit(1,i)/0.2];
end
%exclude the dot outside of the ROI
Xfitnew = [];
Yfitnew = [];
Zfitnew = [];
q = 0;
for i=1:length(x)
    if x(1,i) >= ROI.x1 && x(1,i) <= ROI.x2 && y(1,i) >= ROI.y1 && y(1,i) <= ROI.y3 %inside FOV coordinates
        q = q+1;
        Xfitnew(1,q) = x(1,i);
        Yfitnew(1,q) = y(1,i);
        Zfitnew(1,q) = z(1,i);
    end
end
% normalize the dot position to the crop
dotx = [];
doty = [];
dotz = [];
for i = 1:length(Xfitnew)
    dotx = [dotx Xfitnew(1,i)-ROI.x1];
    doty = [doty Yfitnew(1,i)-ROI.y1];
    dotz = [dotz Zfitnew(1,i)];
end
Dot.x = dotx;
Dot.y = doty;
Dot.z = dotz;
Dot.x = x;
Dot.y = y;
Dot.z = z;
end
% alphashaping Territory
function struct = maskshape(Territory)

ID = unique(Territory(Territory>0));

for label = 1:length(ID)
    struct(label).label = label;
    single_Mask = Territory == label;
    struct(label).nuc_props = regionprops3(single_Mask, 'all');
    [Y, X, Z] = ind2sub( size(single_Mask),...
                                    find(single_Mask) );
    struct(label).mask_Pix = [X Y Z];
    struct(label).alphaShape = alphaShape(X,Y,Z);
end

end
% alphashaping territory correcting the DeltaZ
function struct = maskshapeDeltaZ(Territory,DeltaZ)

ID = unique(Territory(Territory>0));

for label = 1:length(ID)
    struct(label).label = label;
    single_Mask = Territory == label;
    struct(label).nuc_props = regionprops3(single_Mask, 'all');
    [Y, X, Z] = ind2sub( size(single_Mask),...
                                    find(single_Mask) );
    Z = Z-DeltaZ;
    struct(label).mask_Pix = [X Y Z];
    struct(label).alphaShape = alphaShape(X,Y,Z);
end

end
% do quality control for empty shape
% volume cutoff is based on channel
function ReshapeStruct = DoQualityControl(maskshape,channel)

%add the condition segmented by cellpose

switch channel
    case 405
        volumecutoff = 500;
    case 647
        volumecutoff = 100;
end

k = 0;
for i = 1:length(maskshape)
    %if alphashape exist and the volumn of that shape is larger than 1000
    if maskshape(i).alphaShape.Alpha ~= Inf & maskshape(i).nuc_props.Volume > volumecutoff
        k = k+1;
        ReshapeStruct(k) = maskshape(i);
    end
end
end
% do quality control for the dot, do dot in the L1 shape.
% restrict the number of dots in each nuclei no more than 4
function Dot = DotQualityControl(foci,TerriStruct)

l = 0;
for i = 1:length(TerriStruct)
    X = [];
    Y = [];
    Z = [];
   for j = 1:length(foci.x)
       TF = inShape(TerriStruct(i).alphaShape, foci.x(j),foci.y(j),foci.z(j));
           if TF == 1 %if the TAD dots are in the shape
              X = [X foci.x(j)];
              Y = [Y foci.y(j)];
              Z = [Z foci.z(j)];
           end
   end
   % need a quality control for dots number over 4 in each nucleus
   % empty those dots number over 4
   if length(X) <= 4
      for ii = 1:length(X)
          l = l+1;
          Dot.x(l) = X(ii);
          Dot.y(l) = Y(ii);
          Dot.z(l) = Z(ii);
      end
   end
end

end
% associating mRNA to Nuc
function [pha4Nuc, nonpha4Nuc] = classify_mRNA_nuc(mRNA, AllNuc,Set)
mRNAtoNuc = {};
for i = 1:length(AllNuc)
    D = [];
    for j = 1:length(mRNA.x) 
        tf = inShape(AllNuc(i).alphaShape, mRNA.x(j),mRNA.y(j),mRNA.z(j));
        if tf == 0
            I = nearestNeighbor(AllNuc(i).alphaShape,mRNA.x(j),mRNA.y(j),mRNA.z(j));
            near_dot_x = AllNuc(i).alphaShape.Points(I,1)*0.11;
            near_dot_y = AllNuc(i).alphaShape.Points(I,2)*0.11;
            near_dot_z = AllNuc(i).alphaShape.Points(I,3)*0.2;

            d = ((mRNA.x(j)*0.11 - near_dot_x).^2+ ...
                 (mRNA.y(j)*0.11 - near_dot_y).^2+ ... 
                 (mRNA.z(j)*0.2 - near_dot_z).^2).^0.5;
        end
        D = [D d]; 

    end
    mRNAtoNuc = [mRNAtoNuc D];
end        
mRNAtoNuc = mRNAtoNuc.';

m=0;
n=0;
for i = 1:length(mRNAtoNuc)
    A = find(mRNAtoNuc{i,1} < Set.mRNA_dis_thresh);
    Num = length(A);
    if Num >= Set.mRNA_Num_thresh
       m = m+1;
       pha4Nuc(m).label = AllNuc(i).label;
       pha4Nuc(m).struct = AllNuc(i);
    else
       n = n+1;
       nonpha4Nuc(n).label = AllNuc(i).label;
       nonpha4Nuc(n).struct = AllNuc(i);
    end
end

end
% associating nRNA ID to Nuc
function List_Nuc_Intron = classify_intron_nuc(nRNA, NucStruct, nRNA_dis_thresh)
k = 0;
for i = 1:length(NucStruct)
    IntronID = [];
    for j = 1:length(nRNA.x) 
        tf = inShape(NucStruct(i).alphaShape, nRNA.x(j),nRNA.y(j),nRNA.z(j));
        if tf == 0 % if the nRNA is not in nuc but with distance <0.5 um
            I = nearestNeighbor(NucStruct(i).alphaShape,nRNA.x(j),nRNA.y(j),nRNA.z(j));
            near_dot_x = NucStruct(i).alphaShape.Points(I,1)*0.11;
            near_dot_y = NucStruct(i).alphaShape.Points(I,2)*0.11;
            near_dot_z = NucStruct(i).alphaShape.Points(I,3)*0.2;

            d = ((nRNA.x(j)*0.11 - near_dot_x).^2+ ...
                 (nRNA.y(j)*0.11 - near_dot_y).^2+ ... 
                 (nRNA.z(j)*0.2 - near_dot_z).^2).^0.5;
            if d < nRNA_dis_thresh
               IntronID = [IntronID j];
            end
        elseif tf == 1 % if the nRNA is in nuc
               IntronID = [IntronID j];
        end
    end
    if ~isempty(IntronID)
    k = k+1;
    List_Nuc_Intron(k).NucID = NucStruct(i).label;
    List_Nuc_Intron(k).Num = length(IntronID);
    List_Nuc_Intron(k).Intron = IntronID;
    end
end        
end
% Association territory to the list of nuclei
function Table_Asso = NucTerrAssociation(NucStruct, Territory,filter_over4terri)
k=0;
for i = 1:length(NucStruct)
    Terri = [];
    for j = 1:length(Territory)
        tf = inShape(NucStruct(i).alphaShape, Territory(j).nuc_props.Centroid(1,1),...
        Territory(j).nuc_props.Centroid(1,2), Territory(j).nuc_props.Centroid(1,3));
        if tf == 1
        Terri = [Terri Territory(j)];
        end
    end
%     filter out the terrinumber over 4 if you want
if filter_over4terri
    if length(Terri)<=4
        k=k+1;
        Table_Asso(k).NucID = NucStruct(i).label;
        Table_Asso(k).Terri = Terri;  
    end
else
       Table_Asso(i).NucID = NucStruct(i).label;
       Table_Asso(i).Terri = Terri; 
end
end
end
% Associate pha-4foci dot to nuclei
function List_Nuc_Dot = NucDot(NucStruct, Dot, Set)

k = 0;
for i = 1:length(NucStruct)
   DotID = [];
   for j = 1:length(Dot.x)
       TF = inShape(NucStruct(i).alphaShape, Dot.x(j),Dot.y(j),Dot.z(j));
           if TF == 1 %if the TAD dots are in the nucleus
           DotID = [DotID j];
           else %if the TAD dots are outside of the nucleus
                I = nearestNeighbor(NucStruct(i).alphaShape,Dot.x(j),Dot.y(j),Dot.z(j));
                near_dot_x = NucStruct(i).alphaShape.Points(I,1)*0.11;
                near_dot_y = NucStruct(i).alphaShape.Points(I,2)*0.11;
                near_dot_z = NucStruct(i).alphaShape.Points(I,3)*0.2;

                d = ((Dot.x(j)*0.11 - near_dot_x).^2+ ...
                     (Dot.y(j)*0.11 - near_dot_y).^2+ ... 
                     (Dot.z(j)*0.2 - near_dot_z).^2).^0.5;
                %if the distance to the nuclear edge is smaller than 0.5um
                if d < Set.dot_Nuc_thresh
                   DotID = [DotID j];
                end
           end
   end
   k = k+1;
   List_Nuc_Dot(k).BlobID = NucStruct(i).label;
   List_Nuc_Dot(k).DotID = DotID;
   end
end
% classify the pha-4foci dot, Terri2, intron with pha and nonpha nuc
function List_Nuc_Dot_Terri = concatlist(pha4Nuc,List_Nuc_Dot,List_Nuc_Terr2,List_Nuc_Intron)
%associate nuc has dot and terri together
for i = 1:length(List_Nuc_Dot)
    for j = 1:length(List_Nuc_Terr2)
        if List_Nuc_Terr2(j).NucID == List_Nuc_Dot(i).BlobID
           List_Nuc_Dot(i).Terri = List_Nuc_Terr2(j).Terri;
        end
    end
%classify the nuc into the pha nuc and nonpha nuc
    List_pha_nuc = [pha4Nuc.label].';
    if ismember(List_Nuc_Dot(i).BlobID,List_pha_nuc)
       List_Nuc_Dot(i).BelongToPha = 1; %if the nuc is a pha nuc yes then assign 1
    elseif ~ismember(List_Nuc_Dot(i).BlobID,List_pha_nuc)
        List_Nuc_Dot(i).BelongToPha = 0;%if nonpha no then assign 0
    end
end

for i = 1:length(List_Nuc_Dot)
    for j = 1:length(List_Nuc_Intron)
        if List_Nuc_Dot(i).BlobID == List_Nuc_Intron(j).NucID
           List_Nuc_Dot(i).IntronNum = length(List_Nuc_Intron(j).Intron);
           List_Nuc_Dot(i).Intron = List_Nuc_Intron(j).Intron;
        end
    end
end
List_Nuc_Dot_Terri = List_Nuc_Dot;
end
% classify the nearDisTerr2 to pha, nonpha, 2ON, ON, 2OFF calculate the distance between dot28 and Terri2
% then calculate the distance between the edge of terr2 and the dot28. Find
% the shortest distance of all terr2 with dot 28 to be its dis.
function [pha_dis, nonpha_dis, BothON_dis, OneON_dis, BothOFF_dis] = calculateDis(List_Nuc_Dot_Terri, Dot, NucStruct)
% in each nuclei, culculate the dot to each territories, not consider dot
% number over 3 vs 1 territories
for i = 1:length(List_Nuc_Dot_Terri)
        NearDisTerr2=[];
        %fetch the dot in That nucleus
        if ~isempty(List_Nuc_Dot_Terri(i).DotID)
        %if the DotID is not empty, itterate for each dot in that nucleus 
        % and calculate all the nearest distance to all the territories in 
        % the nucleus
        for DotID = List_Nuc_Dot_Terri(i).DotID
            D = [];
            if ~isempty(List_Nuc_Dot_Terri(i).Terri)
                for j = 1:length(List_Nuc_Dot_Terri(i).Terri)
                    tf = inShape(List_Nuc_Dot_Terri(i).Terri(j).alphaShape,...
                        Dot.x(DotID),Dot.y(DotID),Dot.z(DotID));
                    if tf == 1
                       d = 0;
                    else
                    I = nearestNeighbor(List_Nuc_Dot_Terri(i).Terri(j).alphaShape,...
                        Dot.x(DotID),Dot.y(DotID),Dot.z(DotID));
                    near_dot_x = List_Nuc_Dot_Terri(i).Terri(j).alphaShape.Points(I,1)*0.11;
                    near_dot_y = List_Nuc_Dot_Terri(i).Terri(j).alphaShape.Points(I,2)*0.11;
                    near_dot_z = List_Nuc_Dot_Terri(i).Terri(j).alphaShape.Points(I,3)*0.2;
        
                    d = ((Dot.x(DotID)*0.11 - near_dot_x).^2+ ...
                         (Dot.y(DotID)*0.11 - near_dot_y).^2+ ... 
                         (Dot.z(DotID)*0.2 - near_dot_z).^2).^0.5;
                    end
                    D = [D d];
                end
                dmin = min(D);
                %set cutoff distance based on the radius
                ind = find([NucStruct.label]==List_Nuc_Dot_Terri(i).BlobID);
                V = NucStruct(ind).nuc_props.Volume;
                R = (0.75*V/3.1416)^(1/3);
                %set the threshold to radius and only save the distance
                %less than R
                if dmin <= R
                NearDisTerr2 =[NearDisTerr2 dmin];
                end
            end
         end
         end
   List_Nuc_Dot_Terri(i).NearDisTerr2 = NearDisTerr2;
end

nonpha_dis = [];
pha_dis = [];
BothON_dis = [];
OneON_dis = [];
BothOFF_dis = [];

for i = 1:length(List_Nuc_Dot_Terri)
    % nonpha distance
    if List_Nuc_Dot_Terri(i).BelongToPha == 0
       nonpha_dis = [nonpha_dis List_Nuc_Dot_Terri(i).NearDisTerr2];
    % pha distance
    elseif List_Nuc_Dot_Terri(i).BelongToPha == 1
        pha_dis = [pha_dis List_Nuc_Dot_Terri(i).NearDisTerr2];
       if  isfield(List_Nuc_Dot_Terri,'IntronNum') 
            if isempty(List_Nuc_Dot_Terri(i).IntronNum) % 2OFF distance
               BothOFF_dis = [BothOFF_dis List_Nuc_Dot_Terri(i).NearDisTerr2];
            end
       end   
    end
    if  isfield(List_Nuc_Dot_Terri,'IntronNum')
        if List_Nuc_Dot_Terri(i).IntronNum >= 2 %both ON distance
           BothON_dis = [BothON_dis List_Nuc_Dot_Terri(i).NearDisTerr2];
        % 1ON distance
        elseif List_Nuc_Dot_Terri(i).IntronNum == 1
           OneON_dis = [OneON_dis List_Nuc_Dot_Terri(i).NearDisTerr2];
        end
    elseif ~isfield(List_Nuc_Dot_Terri,'IntronNum')
        BothON_dis = [];
        OneON_dis = [];
    end
end
end
% plot the figure with pha-4+ nuclei, pha-4-nuclei, with region territories
function display(mRNA,nRNA, pha4Nuc, nonpha4Nuc,Dot,List_Nuc_Dot_Terri,SampleNum,FOV,FigurePath,List_Nuc_Intron)

figure
hold on
for i = 1:length(pha4Nuc)
    h{1} = plot(pha4Nuc(i).struct.alphaShape  ,'FaceColor','red','FaceAlpha',0.05,'LineStyle', 'none');
end
for i = 1:length(nonpha4Nuc)
    h{2} = plot(nonpha4Nuc(i).struct.alphaShape  ,'FaceColor','blue','FaceAlpha',0.05,'LineStyle', 'none');
end
for i = 1:length(List_Nuc_Dot_Terri)
    for j = 1:length(List_Nuc_Dot_Terri(i).Terri)
    h{3} = plot(List_Nuc_Dot_Terri(i).Terri(j).alphaShape  ,'FaceColor',[0.4940 0.1840 0.5560],'FaceAlpha',0.3,'LineStyle', 'none');
    end
end
h{4} = scatter3(mRNA.x(1,:), mRNA.y(1,:), mRNA.z(1,:) , '.','blue');
% only fetch the nRNA in nuclei or within certain distance cutoff
ID_nRNA_qualified = [List_Nuc_Intron.Intron];
h{5} = scatter3(nRNA.x(1,ID_nRNA_qualified), nRNA.y(1,ID_nRNA_qualified), nRNA.z(1,ID_nRNA_qualified) ,'og', 'MarkerFaceColor', 'g');
h{6} = scatter3(Dot.x(1,:),Dot.y(1,:),Dot.z(1,:),'or', 'MarkerFaceColor', 'r');
legend([h{1}, h{2}, h{3}, h{4}, h{5}, h{6}], {'pha-nuc', 'nonpha-nuc','A-compart','mRNA','nRNA','pha-4region' } );
% for i = 1:length(Dot.x)
% text(Dot.x(i)+2, Dot.y(i),Dot.z(i), num2str(i), 'Color', 'Black', ...
% 'FontWeight', 'Bold', 'FontSize',10);
% end
% for i = 1:length(pha4Nuc)
%     px = pha4Nuc(i).struct.nuc_props.Centroid(1);
%     py = pha4Nuc(i).struct.nuc_props.Centroid(2);
%     pz = pha4Nuc(i).struct.nuc_props.Centroid(3);    
% text(px+2,py,pz, num2str(pha4Nuc(i).struct.label), 'Color', 'red', ...
% 'FontWeight', 'Bold', 'FontSize',10);
% end
% for i = 1:length(nonpha4Nuc)
%     nx = nonpha4Nuc(i).struct.nuc_props.Centroid(1);
%     ny = nonpha4Nuc(i).struct.nuc_props.Centroid(2);
%     nz = nonpha4Nuc(i).struct.nuc_props.Centroid(3);
% text(nx+2, ny,nz, num2str(nonpha4Nuc(i).struct.label), 'Color', 'green', ...
% 'FontWeight', 'Bold', 'FontSize',10);
% end
set(gca,'Ydir','reverse')
title(['FOV' num2str(SampleNum) '_Emb_' num2str(FOV) 'Reconstruction'])
hold off
saveas(gcf,fullfile(FigurePath,['FOV-' num2str(SampleNum) '-Emb-' num2str(FOV)]));
end
