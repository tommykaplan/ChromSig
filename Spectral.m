files=dir('ZLD_1hr_peaks.2000_oriented.H*.tab.gz'); for i=1:length(files), tsv2mat(files(i).name); end;

[chr,from,to,name,~,str] = textread('ZLD_1hr_peaks.2000_oriented.bed','%s%d%d%s%s%s%*[^\n]','headerlines',0,'delimiter','\t');
% Load tracks for c14
A=load('ZLD_1hr_peaks.2000_oriented.H3K27ac-c14.mat'); XX(:,:,1)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K4me1-c14.mat'); XX(:,:,2)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K27me3-c14.mat');XX(:,:,3)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K18ac-c14.mat'); XX(:,:,4)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K4me3-c14.mat'); XX(:,:,5)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
% Load tracks for c13
A=load('ZLD_1hr_peaks.2000_oriented.H3K27ac-c13.mat'); XX(:,:,6)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K4me1-c13.mat'); XX(:,:,7)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K27me3-c13.mat');XX(:,:,8)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K18ac-c13.mat'); XX(:,:,9)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
A=load('ZLD_1hr_peaks.2000_oriented.H3K4me3-c13.mat'); XX(:,:,10)=A.X; if ~strcmp(name{1232}, A.name{1232}), error(''); end;
XX=single(XX);
hdrs={'H3K27ac-c14','H3K4me1-c14','H3K27me3-c14','H3K18ac-c14','H3K4me3-c14',...
     'H3K27ac-c13','H3K4me1-c13','H3K27me3-c13','H3K18ac-c13','H3K4me3-c13'};
% Load (wrong order!)
clear A;

K = 3;
N = size(XX,1);
len = size(XX,2);
d = size(XX,3);

% distance to adjacency
for i=1:d,
    dist(:,:,i)=squareform(pdist(XX(:,:,i),'euclidean'))/sqrt(len);
    sigma(i)=prctile(reshape(dist(:,:,i),1,N*N),10);
    tmp = exp(-dist(:,:,i).^2/(2*sigma(i).^2)); tmp(find(eye(N)))=0; wei(:,:,i)=tmp;
    if 1,
	clf; xx = 0:.1:100; yy = exp(-xx.^2/(2*sigma(i).^2));
	zz = histc(reshape(dist(:,:,i),1,N*N),xx); zz = zz / max(zz); xx2 = xx(xx<=sigma(i)); zz2 = zz(xx<=sigma(i));
	plot(xx,yy,'-','Color',[37,106,225]/255,'LineWidth',2); hold on;
	plot(xx,zz,'-','Color',[141,56,199]/255,'LineWidth',2); bar(xx2,zz2,2,'FaceColor',[141,56,199]/255,'EdgeColor','none');
	plot([sigma(i) sigma(i) 0],[0 exp(-0.5) exp(-0.5)],'r-','LineWidth',3); hold off;
	axis([0.1 30 0 1]); xlabel('Distance'); ylabel('Weight'); title(hdrs{i})
	set(gcf,'PaperSize',[3 2],'PaperPosition',[0 0 3 2]); print(gcf,'-dpdf','-r300',sprintf('Dist-Weight.%d.pdf',i));

    end
end
A = sum(wei,3); clear dist wei xx yy zz xx2 zz2;

D=diag(sum(A,2));
isD=inv(sqrt(D));
L = double(eye(N) - isD * A * isD);

[t1,dd]=eigs(L,K,'sm'); dd=diag(dd)';
% normalize to the unit sphere
n1=sqrt(sum(t1.^2,2)); t2 = t1./repmat(n1,1,K); clear n1;
[idx,~,sumd]=kmeans(t2,K,'replicates',100,'Start','cluster'); clear sumd;
scatter3(t2(:,1),t2(:,2),t2(:,3),25,idx,'filled');

subplot(2,3,4); plot(permute(mean(XX(idx==1,:,1:5),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);
subplot(2,3,5); plot(permute(mean(XX(idx==2,:,1:5),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);
subplot(2,3,6); plot(permute(mean(XX(idx==3,:,1:5),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);
subplot(2,3,1); plot(permute(mean(XX(idx==1,:,6:10),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);
subplot(2,3,2); plot(permute(mean(XX(idx==2,:,6:10),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);
subplot(2,3,3); plot(permute(mean(XX(idx==3,:,6:10),1),[3,2,1])','LineWidth',2); set(gca,'YLim',[0 80]);

clear t1 L isD fid j I i;
save('spectral.mat');

% dump bed
fid=fopen('ZLD_1hr_peaks.2000_oriented.cl.bed','w');
for j=1:K
    I=find(idx==j)';
    for i=[I],
	fprintf(fid,'%s\t%d\t%d\t%s\t1000\t%s\n',chr{i},from(i),to(i),name{i},str{i});
    end;
    fprintf(fid,'#cl%d\n',j);
end
fclose(fid);
clear i j I fid

