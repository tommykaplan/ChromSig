function [I,LL,RM,T,mu,sigma] = EM_flip_3D(X)
[N,K,D]=size(X);
for d=1:D, [mu(d,:),sigma(d,:)]=normfit(X(:,:,d)); end;
LL0=mean(mean(sum(log2(1e-10+normpdf(zscore(X))),2))),;
RM0=mean(mean(sqrt(mean((X-repmat(mean(X,1),N,1)).^2,2)),3));

% init randomly
I=rand(N,1)<0.5;

for i=1:20,
    % flip
    for d=1:D, T(:,:,d)=X(:,:,d); T(I==0,:,d)=X(I==0,end:-1:1,d); end;

    % Likelihood
    LL(i)=mean(mean(sum(log2(1e-10+normpdf(zscore(T))),2)));
    RM(i)=mean(mean(sqrt(mean((T-repmat(mean(T,1),N,1)).^2,2)),3));

    % M-step
    for d=1:D, [mu(d,:),sigma(d,:)]=normfit(T(:,:,d)); sigma(d,:)=repmat(mean(sigma(d,:)),1,K); end;
    plot(1:K,mu); shg; pause(0.1);

    % E-step
    for d=1:D, t1(:,:,d) = repmat(mu(d,:),N,1); t2(:,:,d) = repmat(sigma(d,:),N,1); end;
    LLf = mean(sum(log2(1e-10+normpdf((X - t1) ./ t2)),3),2);
    LLr = mean(sum(log2(1e-10+normpdf((X(:,end:-1:1,:) - t1) ./ t2)),3),2);

    % Well, actually MM
    I = LLf>LLr;
end
% rebuild
for d=1:D, T(:,:,d)=X(:,:,d); T(I==0,:,d)=X(I==0,end:-1:1,d); end;
LL = [LL0 LL mean(mean(sum(log2(1e-10+normpdf(zscore(T))),2)))];
RM = [RM0 RM mean(mean(sqrt(mean((T-repmat(mean(T,1),N,1)).^2,2)),3))];

for d=1:D, [mu(d,:),sigma(d,:)]=normfit(T(:,:,d)); sigma(d,:)=repmat(mean(sigma(d,:)),1,K); end;
