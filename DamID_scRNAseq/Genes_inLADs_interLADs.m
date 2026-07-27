close all
clear all
clc

% find genes in LADs vs inter LADs in early embryo.
T1 = readtable("EMR-1_DamID-med1p-bound-bed-1707061998.txt");
T2 = readtable("Annotation_Gene.txt");

T = removevars(T2,["Var2","Var3","Var6","Var7","Var8","Var9","Var11","Var12","Var13","Var14","Var15","Var17"]);
T = renamevars(T,"x__genebuild_versionWS285","chr");
T = renamevars(T,"Var4","start");
T = renamevars(T,"Var5","end");
T = renamevars(T,"Var10","WB");
T = renamevars(T,"Var16","kind");
T = renamevars(T,"Var18","gene");
T1 = renamevars(T1,"xEnd","end");
head(T,5)
clear T2

%remove non protein coding gene
id = [];
for i=1:height(T)
T.mid(i) = mean(T.start(i),T.end(i));
    if ~strcmp(T.kind(i),'protein_coding') 
       id = [id i];
    end
end
T(id,:) = [];
%remove mtDNA gene
id = [];
for i=1:height(T)
    if strcmp(T.chr(i),'MtDNA') 
       id = [id i];
    end
end
T(id,:) = [];

for i = 37:height(T)
    chr = string(T.chr(i));
    if ~isempty(chr{1,1})
        if chr == "I"
            for ii = 1:403
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end
        elseif chr == "II"
            for ii = 404:880
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end
        elseif chr == "III"
            for ii = 881:1236
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end
        elseif chr == "IV"
            for ii = 1237:1873
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end
        elseif chr == "V"
            for ii = 1874:2741
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end
        elseif chr == "X"
            for ii = 2742:3113
                if T1.start(ii)<=T.mid(i) & T1.end(ii)>=T.mid(i)
                    T.inLADs(i) = 1;
                %first test if gene include a 10kb region, thus at least 50% of
                % the gene is in the region.
                elseif T.start(i)<=T1.start(ii) & T.end(i)>=T1.end(ii)
                    T.inLADs(i) = 1;
                end
                T.DamID(i) = T1.average(ii);
            end        
        end
    end
end

% %see the expression level of each stage.
list = [T.inLADs];
lia1 = find(list == 1); % in LADs genes
lia2 = find(list == 0); % inter LADs genes

TinLADs = T;
TinterLADs = T;

%extract inLADs genes
id = [];
for i=1:height(TinLADs)
    if TinLADs.inLADs(i) == 0
       id = [id i];
    end
end
TinLADs(id,:) = [];
%save genes in LADs
writetable(TinLADs,'Genes_inLADs','Delimiter',',');

%extract inLADs genes
id = [];
for i=1:height(TinterLADs)
    if TinterLADs.inLADs(i) == 1
       id = [id i];
    end
end
TinterLADs(id,:) = [];
%save genes in LADs
writetable(TinterLADs,'Genes_interLADs','Delimiter',',');

names = ["genes inLADs";"genes outLADs"];
number = [6163; 13806];
tbl = table(number,names);
piechart(tbl,"number","names")
mypalette = validatecolor(["#D95319" "#0072BD"], "multiple");
set(gca, 'LineWidth', 2, 'FontSize', 12, 'FontName', 'Arial');
colororder(mypalette)
