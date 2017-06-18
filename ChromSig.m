files=dir('*.tab.gz'); for i=1:length(files), tsv2mat(files(i).name); end;

% Load tracks for c14
load('ZLD_1hr_peaks.2000.H3K27ac-c14.mat'); XX(:,:,1)=X;
load('ZLD_1hr_peaks.2000.H3K4me1-c14.mat'); XX(:,:,2)=X;
load('ZLD_1hr_peaks.2000.H3K27me3-c14.mat');XX(:,:,3)=X;
load('ZLD_1hr_peaks.2000.H3K18ac-c14.mat'); XX(:,:,4)=X;
load('ZLD_1hr_peaks.2000.H3K4me3-c14.mat'); XX(:,:,5)=X;

[I,LL,RM,T,mu,sigma] = EM_flip_3D(XX);
% Average improvement in Likeihood (per sequence/timepoint)
2.^(LL(end)-LL(1))

N=length(I);
for i=1:N, if I(i), s{i}='+'; else s{i}='-'; end; end;
fid=fopen('ZLD_1hr_peaks.2000_oriented.bed','w');
for i=1:N,
    fprintf(fid,'%s\t%d\t%d\t%s\t0\t%s\n',chr{i},from(i),to(i),name{i},s{i});
end;
fclose(fid);

% expression
X = dlmread('ZLD_1hr_peaks.2000_oriented.expression','\t',0,1); % X(end,:)=[];
for i=1:3, I=X(:,1)==i;, T(i,:)=mean(X(I,2:end)); end;
clf; 
bar(T',2); axis([0.5 8.5 0 275]); set(gca,'XTick',[1:8],'XTickLabel',{'10','11','12','13','14A','14B','14C','14D'});
legend({'Cluster 1','Cluster 2','Cluster 3'},'location','nw'); ylabel('Average Expression');
set(gcf,'PaperSize',[4.2 3],'PaperPosition',[0 0 4.2 3]); print(gcf,'-dpdf','-r300',sprintf('Exp.pdf'));
