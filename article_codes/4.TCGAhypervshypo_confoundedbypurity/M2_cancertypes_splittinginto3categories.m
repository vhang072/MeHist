cancer_types = { 'CESC','LGG','COAD','CHOL','ESCA','HNSC','READ','UCEC','GBM','BRCA','PRAD','BLCA','LUAD','LUSC','LIHC','LAML','PAAD',...
   'SKCM','ACC','STAD', 'SARC','KIRC','KIRP','PCPG','KICH','THCA','THYM'};
              
figure(61);    
id1 = 3;
id2 = 4;
vals = [];
for x = 1:27%length(cancer_types)
    data0 = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\ABSOLUTE_MeanCGIprobe\_SigPur.',...
    char(cancer_types(x)),".CGIMean_ABSOLUTE.txt"),...
        'filetype','text','readvariablenames',false, 'delimiter','\t','headerlines',1,'readrownames',false,...
         'TreatAsEmpty','NA');
     numres = data0{:,2:6};
     
     numres = [numres, strcmp(data0{:,end},"Included") *2.5+3.0];
     numres0 = data0{strcmp(data0{:,end},"Included"),2:6};
     numres = numres(numres(:,5)> min(numres0(:,5)),:);

     subplot(3,9,x);%histogram2(numres(:,id2),numres(:,id1), 0:800:8000,0:3000:30000,'DisplayStyle','tile','ShowEmptyBins','off');
     h1 = histogram(numres(:,id2)*2+numres(:,id1),0:4000:40000,'Normalization', 'probability');
     vals = [vals,sum(numres(:,id2)*2+numres(:,id1) < 8000)/size(numres,1)];
     title(char(cancer_types(x)));
     colormap(jet);

end
out = [cell2table(cancer_types'), table(vals','VariableNames',{'val'})];

[~,ind] = sort(out{:,2},"ascend");
out_sort = out(ind,:);
figure(20);
y = out_sort{:,2};
x = categorical(table2array(out_sort(:,1)));
xx = reordercats(x,table2array(out_sort(:,1)));
bar(xx',y')
