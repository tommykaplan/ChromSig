load('spectral.mat','XX','idx','K','d','hdrs');

% prepare signatures
hdrs = hdrs(1:5); d=5; XX=XX(:,:,1:5);
for k=1:K,
    heights = permute(mean(XX(idx==k,:,:),2),[1,3,2]);
    for i=1:d,
        % 20 bin histogram
	D=mean(XX(idx==k,:,i),2); bin=(max(D)-min(D))/20;
	x=[max(D)-[20:-1:0]*bin]; x(1)=0; x(end)=Inf; h=histc(D,x)'; h(end)=[];
	h=100*h./sum(h);
	HistC{i}{k} = h; % counts
	HistE{i}{k} = x; % edges
    end
    sig(:,:,k) = double(permute(mean(XX(idx==k,:,:),1),[3,2,1]));
    sig(:,:,k) = sig(:,end:-1:1,k);
end
clear XX;

% Load tracks for c14
A=load('../mat/H3K27ac-c14.mat');  G{1}=A.chip;
A=load('../mat/H3K4me1-c14.mat');  G{2}=A.chip;
A=load('../mat/H3K27me3-c14.mat'); G{3}=A.chip;
A=load('../mat/H3K18ac-c14.mat');  G{4}=A.chip;
A=load('../mat/H3K4me3-c14.mat');  G{5}=A.chip;
% A=load('../mat/H3K27ac-c13.mat');  G{6}=A.chip;
% A=load('../mat/H3K4me1-c13.mat');  G{7}=A.chip;
% A=load('../mat/H3K27me3-c13.mat'); G{8}=A.chip;
% A=load('../mat/H3K18ac-c13.mat');  G{9}=A.chip;
% A=load('../mat/H3K4me3-c13.mat');  G{10}=A.chip;
clear A;

chrs = fieldnames(G{1});

step=10;
for c=1:length(chrs),
    X = []; for i=1:d, X(i,:)=double(G{i}.(chrs{c})); end
    len = floor((size(X,2)-length(sig))/step);
    % correlations
    cr=NaN*ones(d,len,K,'single');

    for k=1:len,
	for i=1:d,
	    for class=1:K,
	        % current window
		XX = 1+step*(k-1)+[0:999];
		% correlation with signature of class i
		[cc,pp]=corr(X(i,XX)', sig(i,:,class)');
		% prob
		mPP = max([0,HistC{i}{class}(find(HistE{i}{class}<mean(X(i,XX)),1,'last'))]);
		cr(i,k,class) = cc * mPP;
	    end
	end
	if mod(k,10)==0, fprintf('\r %s %d/%d                    ',chrs{c},step*k,size(X,2)); end;
    end
    clear k i mPP class
    totcorr.(chrs{c})=permute(sum(cr,1),[3,2,1]);
end
fprintf('\r                                \n');

% dump to bedGraph
for class=1:K, fid(class)=fopen(sprintf('Scan.cl%d.-.bedGraph',class),'w'); end
for c=1:length(chrs),
    len=size(totcorr.(chrs{c}),2);
    % zero NaN values
    I=isnan(totcorr.(chrs{c})(class,:)); totcorr.(chrs{c})(class,I)=0;
    for k=1:len
	% find middle point
	pos = 10*step*(k-1) + 5000 - 10*step/2;
	for class=1:K, fprintf(fid(class), '%s\t%d\t%d\t%.1f\n', chrs{c}, pos, pos+10*step, totcorr.(chrs{c})(class,k)); end
    end
end
for class=1:K, fclose(fid(class)); end;

save('Scan.-.mat');
