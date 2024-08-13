cancer_types = {'ACC', 'BLCA', 'BRCA','CESC','CHOL', 'COAD','DLBC', 'ESCA', 'GBM','HNSC',   'KICH', 'KIRC', 'KIRP', ...
                'LAML','LGG','LIHC','LUAD',...
                  'LUSC','PAAD', 'PRAD', 'PCPG','READ','SARC','SKCM','STAD', 'THCA','THYM', 'UCEC'};
              
figure(1);           
for x = 1:length(cancer_types)
    data = readtable(strcat('Y:\4.basic_data\TCGA_PancanAtlas\methylation_cgiprobe\Diff_TCGA\_DiffCGIprobeMean_',char(cancer_types(x)),".txt"),...
        'filetype','text','readvariablenames',true, 'delimiter','\t','headerlines',0,'readrownames',true,...
         'TreatAsEmpty','NA');
     numres = data{:,:};
     subplot(5,6,x);scatter(numres(2,:),numres(1,:),10,"filled","markeredgecolor","k");title(char(cancer_types(x)));
     xlim([min(numres(2,:))-0.1,max(numres(2,:))+0.1]);ylim([min(numres(1,:))-0.1,max(numres(1,:))+0.1])
 
end