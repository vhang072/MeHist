cancer_types = {'ACC', 'BLCA', 'BRCA','CESC', 'COAD', 'GBM','HNSC',   'KICH', 'KIRC', 'KIRP', ...
                'LGG','LIHC','LUAD',...
                  'LUSC',  'PRAD','READ','SKCM', 'THCA', 'UCEC'};
cancer_types = {'LUSC','CESC','UCEC','GBM','HNSC','SKCM','BLCA','LIHC','READ','LUAD','COAD','ACC',...
    'BRCA','LGG','THCA','KIRP','KICH','PRAD','KIRC'};
              
figure(4);           
id1 = 4;
id2 = 3;
spnum = [];
for x = 1:length(cancer_types)
    data = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\CPE_MeanCGIprobe\',...
    char(cancer_types(x)),".CGIMean_CPE.txt"),...
        'filetype','text','readvariablenames',true, 'delimiter','\t','headerlines',0,'readrownames',false,...
         'TreatAsEmpty','NA');
     numres = data{:,:};
     [~,idxxx] = sort(numres(5,:),'ascend');
     numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(2)];
     colormap(jet);
     subplot(3,7,x);scatter(numres(id1,:),numres(id2,:),6,numres(5,:),"filled","markeredgecolor","k",'LineWidth',0.1);title(char(cancer_types(x)));
     caxis([0.5 1]);
     hold on;
     x1 = [min(numres(id1,:)),max(numres(id1,:))];
     x2 = [min(numres(id2,:)),max(numres(id2,:))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
 
end

cancer_types = cancer_types(spnum>=200);
figure(42);    
id1 = 4;
id2 = 3;
spnum = [];
for x = 1:length(cancer_types)
    data = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\ABSOLUTE_MeanCGIprobe\',...
    char(cancer_types(x)),".CGIMean_ABSOLUTE.txt"),...
        'filetype','text','readvariablenames',true, 'delimiter','\t','headerlines',0,'readrownames',false,...
         'TreatAsEmpty','NA');
     numres = data{:,:};
     [~,idxxx] = sort(numres(5,:),'ascend');
     numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(2)];
     colormap(jet);
     subplot(3,7,x);scatter(numres(id1,:),numres(id2,:),6,numres(5,:),"filled","markeredgecolor","k",'LineWidth',0.1);title(char(cancer_types(x)));
     caxis([0.2 1]);
     hold on;
     x1 = [min(numres(id1,:)),max(numres(id1,:))];
     x2 = [min(numres(id2,:)),max(numres(id2,:))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
 
end