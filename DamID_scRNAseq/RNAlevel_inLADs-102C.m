clear all

T1 = readtable("Annotation_Gene.txt");
T3 = readtable("EMR-1_DamID-med1p-bound-bed-1707061998.txt");
T2 = readtable("AlisonMSxppp_over25.xlsx");

AnoGenArray = [T1.Var18];
AnoWBGenArray = [T1.Var10];
for i = 1:height(T2)
    Gene = string(T2.Var3(i));
    index1 = find(ismember(AnoWBGenArray,Gene));
    if ~isempty(index1)
        T2.Chr(i) = append('chr', T1.x__genebuild_versionWS285(index1));
        T2.start(i) = T1.Var4(index1);
        T2.end(i) = T1.Var5(index1);
        T2.mid(i) = mean([T1.Var4(index1);T1.Var5(index1)]);
    end
end
%Generate the mid of 10kb regions
for i=1:height(T3)
T3.mid(i) = mean(T3.start(i),T3.xEnd(i));
end
%Test if at least 50% of genes locates in 10kb DamID region
for i = 1:height(T2)
chr = T2.Chr(i);
if ~isempty(chr{1,1})
    if chr == "chrI"
        for ii = 1:403
            %first test if mid-gene in a 10kb region, thus at least 50% of
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    elseif chr == "chrII"
        for ii = 404:880
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    elseif chr == "chrIII"
        for ii = 881:1236
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    elseif chr == "chrIV"
        for ii = 1237:1873
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    elseif chr == "chrV"
        for ii = 1874:2741
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    elseif chr == "chrX"
        for ii = 2742:3113
            % the gene is in the region.
            if T3.start(ii)<=T2.mid(i) & T3.xEnd(ii)>=T2.mid(i)
                T2.inLADs(i) = 1;
            %first test if gene include a 10kb region, thus at least 50% of
            % the gene is in the region.
            elseif T2.start(i)<=T3.start(ii) & T2.end(i)>=T3.xEnd(ii)
                T2.inLADs(i) = 1;
            end
            T2.DamID(i) = T3.mid(ii);
        end
    end
end
end
writetable(T2,['DamID_MSxppp'],'Delimiter',',');
