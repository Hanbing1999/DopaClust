function [TPpDAs,pDAs] = test_stable(table_raw,range,k,youfun)
%youfun: K-means, gm,....
load ('D:\graduation thesis\results\table_raw.mat');
% TPpDAs_total = [];
% pDAs_total = [];
% for t = 1:10
TPpDAs = [];
pDAs = [];
%for q = 19:50:289
for q  = range
    table = table_raw{2,1};
    random = randperm(489,q);
    table_train = table (random,:);
    [lower,upper] = bandwidth(table_train);
    for m = 1: upper
    p = find(~isnan(table_train(:,m)));
    [Z,mu,sigma] = zscore(table_train(p,m));
    for n = 1:length(p)
    table_ztrain (p(n),m) = Z(n);
    end 
    end 
    table_ztrain = [table_ztrain table_train(:,upper+1)] ;
    x = table_ztrain (:,1:end-1);
    y = table_ztrain (:,end);
    
    m = 38;
    [U,Z] = pca(x,'NumComponents',m); 
    xTrain = Z;
    yTrain = y;
    
    % k-means
    if youfun == 1
    name = 'k-means';
    idx = kmeans (xTrain,k);
    idx2 = yTrain(find (idx == 2));
    length_idx2 = length(find(idx2==1));
    idx1 = yTrain(find (idx == 1));
    legnth_idx1 = length(find(idx1==1));
    if length_idx2 >legnth_idx1 % idx = 2 is DA 
    else 
    % idx = 1 is DA 
    idx(find(idx ==2)) = 3;
    idx(find(idx ==1)) = 4;
    idx(find(idx ==3)) = 1;
    idx(find(idx ==4)) = 2;
    end 
    end 
    
    if youfun == 4
        rng default  % For reproducibility
        name = 'Gaussian Mixture Model';
        GMModel = fitgmdist(xTrain,2,'CovType','Diagonal');%'SharedCovariance',true 
        idx = cluster(GMModel,xTrain);
        idx(find(idx==1))= 0;
        idx(find(idx==2)) = 3;
        idx(find(idx==3)) = 1;
        idx(find(idx==0)) = 2;
    end
    
    if  youfun == 2
        
    % 1.  kmedoids
    name = 'k-medoids';
    k = 2;
    idx = kmedoids(xTrain,k);
    idx2 = yTrain(find (idx == 2));
    length_idx2 = length(find(idx2==2));
    idx1 = yTrain(find (idx == 1));
    legnth_idx1 = length(find(idx1==1));
    if length_idx2 >legnth_idx1 % idx = 2 is DA 
    else 
    % idx = 1 is DA 
    idx(find(idx ==2)) = 3;
    idx(find(idx ==1)) = 4;
    idx(find(idx ==3)) = 1;
    idx(find(idx ==4)) = 2;
    end 
    end 

     if  youfun == 3
    
    idx = dbscan(xTrain,3,12);
    idx2 = yTrain(find (idx == 2));
    length_idx2 = length(find(idx2==2));
    idx1 = yTrain(find (idx == -1));
    legnth_idx1 = length(find(idx1==-1));
    if length_idx2 >legnth_idx1 % idx = 2 is DA 
    else 
    % idx = -1 is DA 
    idx(find(idx ==2)) = 3;
    idx(find(idx ==-1)) = 4;
    idx(find(idx ==3)) = 1;
    idx(find(idx ==4)) = 2;
    end 
     end 

    if  youfun == 5
    rng(1); % For reproducibility
    name = 'OCSVM';
    yTrain_DA = ones(size(yTrain,1),1);
    SVMModel = fitcsvm(xTrain,yTrain_DA,'KernelFunction','RBF','KernelScale','auto','Standardize',false,'OutlierFraction',0.05);
    svInd =double(SVMModel.IsSupportVector);
    svInd(find(svInd==0))= 2;% 2 is DA
    svInd(find(svInd==1)) = 1; 
    idx = svInd;
    end 

    [TPpDA,pDA] = Clustering_plotting_2(yTrain,idx);
    TPpDAs = [TPpDAs TPpDA];
    pDAs = [pDAs pDA];
% end 
%  TPpDAs_total = [ TPpDAs_total;TPpDAs];
%   pDAs_total = [ pDAs_total; pDAs];
end
end 
