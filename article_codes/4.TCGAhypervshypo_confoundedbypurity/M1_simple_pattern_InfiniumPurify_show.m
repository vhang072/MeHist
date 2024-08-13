cancer_types = {'ACC', 'BLCA', 'BRCA','CESC', 'CHOL','COAD','ESCA', 'GBM','HNSC',   'KICH', 'KIRC', 'KIRP', ...
                'LAML','LGG','LIHC','LUAD',...
                  'LUSC', 'PAAD','PCPG', 'PRAD','SARC','SKCM','STAD', 'THCA','THYM', 'UCEC'};
 cancer_types = { 'LUSC', 'SKCM',  'BLCA',  'CESC',  'STAD', 'LGG','PCPG','LIHC','GBM','HNSC',...  
     'CHOL','UCEC','ESCA','LUAD','COAD','SARC','BRCA','ACC','LAML','PAAD','PRAD','KICH','KIRP',...
     'KIRC','THCA','THYM'};
           
id1 = 4;
id2 = 3;
figure(81);
for x = 1:9%length(cancer_types)
    data0 = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\InfiniumPurify_MeanCGIprobe\_SigPur.',...
    char(cancer_types(x)),".CGIMean_InfiniumPurify.txt"),...
        'filetype','text','readvariablenames',false, 'delimiter','\t','headerlines',1,'readrownames',false,...
         'TreatAsEmpty','NA');
    
     numres = data0{:,2:6};
     
     %numres = [numres, strcmp(data0{:,end},"Included") *2.5+3.0];
     numres = [numres,ones(size(numres,1),1)*5.5];
     
     numres0 = data0{strcmp(data0{:,end},"Included"),2:6};
     numres = numres(numres(:,5)>= 0,:);%min(numres0(:,5))
     %numres = numres0;
     %numres = numres(numres(:,5)> min(numres0(:,5)),:);
     %numres = data{:,:};
     %[~,idxxx] = sort(numres(5,:),'ascend');
     %numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(1)];
     %colormap(parula);
     colormap(jet);
     subplot(1,9,x);scatter(numres(:,id1),numres(:,id2),numres(:,6),numres(:,5),"filled","markeredgecolor","k",'LineWidth',0.01);title(char(cancer_types(x)));
     caxis([0.3 1]);
     hold on;
     x1 = [min(numres(:,id1)),max(numres(:,id1))];
     x2 = [min(numres(:,id2)),max(numres(:,id2))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
     
    hold on; 
    p = polyfit(numres0(:,id1),numres0(:,id2),1);
    xx = [numres0(:,id1),ones(length(numres0(:,id1)),1)];
    ycal = xx*p';
    plot(xx(:,1),ycal,'r:','LineWidth',1); hold off
 
end

figure(82);
for x = 10:18%length(cancer_types)
    data0 = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\InfiniumPurify_MeanCGIprobe\_SigPur.',...
    char(cancer_types(x)),".CGIMean_InfiniumPurify.txt"),...
        'filetype','text','readvariablenames',false, 'delimiter','\t','headerlines',1,'readrownames',false,...
         'TreatAsEmpty','NA');
    
     numres = data0{:,2:6};
     
     %numres = [numres, strcmp(data0{:,end},"Included") *2.5+3.0];
     numres = [numres,ones(size(numres,1),1)*5.5];
     
     numres0 = data0{strcmp(data0{:,end},"Included"),2:6};
     numres = numres(numres(:,5)>= 0,:); %min(numres0(:,5))
     %numres = numres0;
     %numres = numres(numres(:,5)> min(numres0(:,5)),:);
     %numres = data{:,:};
     %[~,idxxx] = sort(numres(5,:),'ascend');
     %numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(1)];
     %colormap(parula);
     colormap(jet);
     subplot(1,9,x-9);scatter(numres(:,id1),numres(:,id2),numres(:,6),numres(:,5),"filled","markeredgecolor","k",'LineWidth',0.01);title(char(cancer_types(x)));
     caxis([0.3 1]);
     hold on;
     x1 = [min(numres(:,id1)),max(numres(:,id1))];
     x2 = [min(numres(:,id2)),max(numres(:,id2))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
     
    hold on; 
    p = polyfit(numres0(:,id1),numres0(:,id2),1);
    xx = [numres0(:,id1),ones(length(numres0(:,id1)),1)];
    ycal = xx*p';
    plot(xx(:,1),ycal,'r:','LineWidth',1); hold off
 
 
end

figure(83);
for x = 19:26%length(cancer_types)
    data0 = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\InfiniumPurify_MeanCGIprobe\_SigPur.',...
    char(cancer_types(x)),".CGIMean_InfiniumPurify.txt"),...
        'filetype','text','readvariablenames',false, 'delimiter','\t','headerlines',1,'readrownames',false,...
         'TreatAsEmpty','NA');
    
     numres = data0{:,2:6};
     
     %numres = [numres, strcmp(data0{:,end},"Included") *2.5+3.0];
     numres = [numres,ones(size(numres,1),1)*5.5];
     
     numres0 = data0{strcmp(data0{:,end},"Included"),2:6};
     numres = numres(numres(:,5)>= 0,:);% min(numres0(:,5))
     %numres = numres0;
     %numres = numres(numres(:,5)> min(numres0(:,5)),:);
     %numres = data{:,:};
     %[~,idxxx] = sort(numres(5,:),'ascend');
     %numres = numres(:,idxxx);
     s1 = size(numres);
     spnum = [spnum,s1(1)];
     %colormap(parula);
     colormap(jet);
     subplot(1,9,x-18);scatter(numres(:,id1),numres(:,id2),numres(:,6),numres(:,5),"filled","markeredgecolor","k",'LineWidth',0.01);title(char(cancer_types(x)));
     caxis([0.3 1]);
     hold on;
     x1 = [min(numres(:,id1)),max(numres(:,id1))];
     x2 = [min(numres(:,id2)),max(numres(:,id2))];
     %plot(x1,x2,"--");
     hold off;
     ext = (x1(2)-x1(1))/10;
     ext2 = (x2(2)-x2(1))/10;
     xlim([x1(1)-ext,x1(2)+ext]);ylim([x2(1)-ext2,x2(2)+ext2])
     
    hold on; 
    p = polyfit(numres0(:,id1),numres0(:,id2),1);
    xx = [numres0(:,id1),ones(length(numres0(:,id1)),1)];
    ycal = xx*p';
    plot(xx(:,1),ycal,'r:','LineWidth',1); hold off
 
 
end