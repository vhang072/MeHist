cancer_types = {'BLCA', 'BRCA','COAD', 'GBM','HNSC', 'KIRC','LUAD',...
                  'LUSC',  'UCEC'};
              
figure(7);           
id1 = 4;
id2 = 3;
spnum = [];
for x = 1:length(cancer_types)
    data = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\ESTIMATE_MeanCGIprobe\',...
    char(cancer_types(x)),".CGIMean_ESTIMATE.txt"),...
        'filetype','text','readvariablenames',true, 'delimiter','\t','headerlines',0,'readrownames',false,...
         'TreatAsEmpty','NA');
     numres = data{:,:};
     [~,idxxx] = sort(numres(5,:),'ascend');
     numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(2)];
     colormap(jet);
     subplot(2,5,x);scatter(numres(id1,:),numres(id2,:),6,numres(5,:),"filled","markeredgecolor","k",'LineWidth',0.1);title(char(cancer_types(x)));
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