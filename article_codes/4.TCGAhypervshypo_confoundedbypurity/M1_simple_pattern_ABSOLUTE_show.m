cancer_types = {'ACC', 'BLCA', 'BRCA','CESC', 'CHOL','COAD','ESCA', 'GBM','HNSC',   'KICH', 'KIRC', 'KIRP', ...
                'LAML','LGG','LIHC','LUAD',...
                  'LUSC', 'PAAD','PCPG', 'PRAD','READ','SARC','SKCM','STAD', 'THCA','THYM', 'UCEC'};
cancer_types = { 'GBM','KICH','PCPG','SKCM','LUSC','CHOL','LGG','CESC','STAD','LIHC','UCEC','ESCA','READ','BLCA','HNSC',...
    'COAD','SARC','ACC','BRCA','LUAD','PRAD','LAML','KIRC','KIRP','THCA','PAAD','THYM'};

              
figure(61);    
id1 = 4;
id2 = 3;
spnum = [];
for x = 1:length(cancer_types)
    data = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\ABSOLUTE_MeanCGIprobe\',...
    char(cancer_types(x)),".CGIMean_ABSOLUTE.txt"),...
        'filetype','text','readvariablenames',true, 'delimiter','\t','headerlines',0,'readrownames',false,...
         'TreatAsEmpty','NA');
     numres = data{:,:};
     %[~,idxxx] = sort(numres(5,:),'ascend');
     %numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(2)];
     colormap(jet);
     subplot(3,9,x);scatter(numres(id1,:),numres(id2,:),4,numres(5,:),"filled","markeredgecolor","k",'LineWidth',0.03);title(char(cancer_types(x)));
     caxis([0.3 1]);
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
figure(62);    
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
     subplot(2,8,x);scatter(numres(id1,:),numres(id2,:),10,numres(5,:),"filled","markeredgecolor","k",'LineWidth',0.1);title(char(cancer_types(x)));
     caxis([0.3 1]);
     hold on;
     x1 = [min(numres(id1,:)),max(numres(id1,:))];
     x2 = [min(numres(id2,:)),max(numres(id2,:))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
 
end