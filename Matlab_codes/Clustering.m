% This main script is for unsupervised clustering to conduct putative DA neuron identification
%% 0. Import path to save clustering results
filename = ['D:\graduation thesis\results\PCA\58 feature\K-means\'];
%% 0. Z-score normalization and divide training and test sample
% z-score the training group 
num_train = 320;
for l =1:2 % l =1, 38 features; l=2,59 features
table_ztrain = [];
random = randperm(489,num_train); % form a random group from raw dataset
table_random = table_raw{l,1}; % Import our feature table
table_train = table_random (random,:); % Get the training feature table
[lower,upper] = bandwidth(table_train); 
    for m = 1: upper
    p = find(~isnan(table_train(:,m)));
    [Z,mu,sigma] = zscore(table_train(p,m)); % Z-score normalization and disregard the nan
        for n = 1:length(p)
        table_ztrain (p(n),m) = Z(n);
        end 
    end 
table_ztrain = [table_ztrain table_train(:,upper+1)] ;
table_raw{l,2} = table_ztrain; % save the z-scored training feature table
table_raw{l,3} = random;  % save the training sample unit number
testing = setdiff(1:489,random);
end
%% Save table_raw
save ('D:\graduation thesis\results\table_raw.mat','table_raw');
%% 1. Load table_raw.mat
load ('D:\graduation thesis\results\table_raw.mat');
%% 1. PCA for dimentaional reduction
l = 2;  % l =1, 38 features; l=2,59 features
m = 38; % define the number of PCs (can be changed based on different dataset): to reach 99% variance
table_ztrain = table_raw{l,2};
x = table_ztrain (:,1:end-1);
y = table_ztrain (:,end);
[U,Z] = pca(x,'NumComponents',m); 
xTrain = Z;
yTrain = y;
table_ztrain_pca = Z;
random = table_raw{l,3}; 
table_random = table_raw{l,1};
testing = setdiff(1:489,random);
table_train = table_random (random,:);
table_testing = table_random(testing,:);
%% 1.2 PCA for dimentaional reduction and t-SNE for visualization
m = 38;
Y = tsne(x,'Algorithm','exact','NumPCAComponents',m);
%% 1.3 Test different number of PCs on its performance
figure;
subplot(2,2,1)
C = cov(x);
e = eig(C);
e = sort (eig(C),'descend');
plot (cumsum(e)/sum(e),'color','k')
ylim([0 1.1])
hold on
%plot (xlim,[0.99,0.99],':','linewidth',1,'color','k')
plot ([38,38],ylim,':','linewidth',1,'color','k')% The red line identifies the significant dimensions that explain 99.9% variance.
xlabel('Principal Component'); ylabel('Variance Explained')
text (40,0.9,'(38,0.993)')
subplot (2,2,3)
%scatter3 (Z(:,1),Z(:,2),Z(:,3),25, y,'filled');
scatter3 (xTrain(find(yTrain == 1),1),xTrain(find(yTrain == 1),2),xTrain(find(yTrain == 1),3),10,yTrain(find(yTrain == 1)),'filled', 'MarkerFaceColor','#FFC61E')
hold on 
scatter3(xTrain(find(yTrain == -1),1),xTrain(find(yTrain == -1),2),xTrain(find(yTrain == -1),3),10,yTrain(find(yTrain == -1)),'filled', 'MarkerFaceColor','#009ADE')
xlabel('PC1'); ylabel('PC2') ;zlabel('PC3');
h = legend('DA phototagged','Non-phototagged');
set(h,'Box','off','position',[0.188,0.3256,0.1128,0.0578]);
view(-47,28)
title('PCA')
subplot (2,2,[2 4])
rng default % for reproducibility
scatter (Y(find(yTrain == 1),1),Y(find(yTrain == 1),2),15,yTrain(find(yTrain == 1)),'filled', 'MarkerFaceColor','#FFC61E')
hold on 
scatter(Y(find(yTrain == -1),1),Y(find(yTrain == -1),2),15,yTrain(find(yTrain == -1)),'filled', 'MarkerFaceColor','#009ADE')
xlabel('t-SNE Dim.1'); ylabel('t-SNE Dim.2');
axis( [-22 25 -22 30])
title('PCA-t-SNE')
% saveas(gca,['D:\graduation thesis\results\PCA\58 feature\''PCA'],'svg')
% set (gcf,'PaperPosition',[-1,10,25,15],'PaperSize',[23 27])
% print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\' 'PCA TSNE with title']);
%% 2. K-means to identify pDA among all units
figure;
%0. k-means 
name = 'k-means'
k = 2; % the number of clusters
idx = kmeans (xTrain,k);
[idx,C] = kmeans (xTrain,k);
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
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,1)
table_raw{l,4} = idx; % save the idx (putative DA: 2, non-putative:1)
%% 2.1 Optional: Different clustering methods 
figure;
%0. k-means 
name = 'k-means'
k = 2;
idx = kmeans (xTrain,k);
[idx,C] = kmeans (xTrain,k);
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
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,1)
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
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,2)
%2. Density-based spatial clustering of applications with noise (DBSCAN)
figure;
name = 'DBSCAN';
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
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,1)
%3. Gaussian Mixture Model 
rng(1); % For reproducibility
name = 'Gaussian Mixture Model';
GMModel = fitgmdist(xTrain,2);%,'SharedCovariance',true
% figure;
% h = gscatter(xTrain(:,1),xTrain(:,2),yTrain);
% hold on
% gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
% g = gca;
% fcontour(gmPDF,[g.XLim g.YLim])
% 
% gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:) x3(:) x4(:) x5(:) x6(:) x7(:) x8(:) x9(:) x10(:) x11(:) x12(:) x13(:) x14(:) x15(:) x16(:) x17(:) x18(:) x19(:) x20(:) x21(:) x22(:) x23(:) x24(:) x25(:) x26(:) x27(:) x28(:) x29(:) x30(:) x31(:) x32(:) x33(:) x34(:) x35(:) x36(:) x37(:) x38(:)]),size(x1)); % Plot the contour of the pdf of the GMM.
% g= gca;
% fcontour(gmPDF,[g.XLim g.YLim]);
% 
% gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(GMModel,[x0 y0]),x,y);
% fcontour(gmPDF,[-6 8 -4 6])
% form clusters
% creat a 2D-form
% gm = gmdistribution(GMModel.mu(:,1:2), GMModel.Sigma(1:2,1:2,:));
% gmPDF = @(x1,x2)reshape(pdf(gm,[x1(:) x2(:)]),size(x1));
% g = gca;
% fcontour(gmPDF,[g.XLim g.YLim])
idx = cluster(GMModel,xTrain);
idx(find(idx==1))= 0;
idx(find(idx==2)) = 3;
idx(find(idx==3)) = 1;
idx(find(idx==0)) = 2;
P = posterior(GMModel,xTrain);
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,2)
% 4. one-class SVM (OCSVM)
figure;
rng(1); % For reproducibility
name = 'OCSVM';
yTrain_DA = ones(size(yTrain,1),1);
SVMModel = fitcsvm(xTrain,yTrain_DA,'KernelFunction','RBF','KernelScale','auto','Standardize',false,'OutlierFraction',0.05);
svInd =double(SVMModel.IsSupportVector);
svInd(find(svInd==0))= 2;% 2 is DA
svInd(find(svInd==1)) = 1; 
idx = svInd;
[TPpDA,pDA] = Clustering_plotting(x,m,yTrain,idx,name,1)
%% 3. Evaluate different clustering methods
rng('default');  % For reproducibility
figure;
subplot(1,2,1)
name = 'k-means';
eva = evalclusters(xTrain,'kmeans','CalinskiHarabasz','KList',[1:6])
plot(eva)
ylabel ('Calinski-Harabasz value');xlabel('Clusters')
title(name)
ylim([10 40])
subplot(1,2,2)
name = 'Gaussian Mixture Model';
eva = evalclusters(xTrain,'gmdistribution','CalinskiHarabasz','KList',[1:6])
plot(eva)
ylabel ('Calinski-Harabasz value');xlabel('Clusters')
ylim([10 40])
title(name)
set (gcf,'PaperPosition',[-1,10,20,10],'PaperSize',[20 10])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' 'compare']);
%% 3.1 Go back to waveforms
tDA_idx = find(yTrain ==1);
pDA_idx = find(idx ==2);
npDA_idx = find(idx ==1);
tDA_idx = random(tDA_idx);
pDA_idx = random(pDA_idx);
npDA_idx = random(npDA_idx);
figure();
tmin = -20;
tmax = 10;
set(gcf,'position',[150 150 300 900]);
subplot (3,1,1)
set(gca,'position',[0.3 0.65 0.4 0.2])
for n = 1: length(tDA_idx)
tDA_wf = FIGURE{tDA_idx(n),4};
plot(tDA_wf,'color','#3C93C2');
hold on
end 
ylim([-1500 1500]);
xlim ([0 33]);
text(tmin,tmax,'tDA','color','#3C93C2')
set(gca,'XColor','white') 
set(gca,'YColor','white')
hold off
subplot (3,1,2)
set(gca,'position',[0.3 0.45 0.4 0.2])
for n = 1: length(pDA_idx )
pDA_wf = FIGURE{pDA_idx (n),4};
plot(pDA_wf,'color','#40AD5A');
hold on
end 
set(gca,'XColor','white') 
set(gca,'YColor','white')
ylim([-1500 1500]);
text(tmin,tmax,'pDA','color','#40AD5A')
hold off
subplot (3,1,3)
set(gca,'position',[0.3 0.25 0.4 0.2])
for n = 1: length(npDA_idx)
npDA_wf = FIGURE{npDA_idx(n),4};
plot(npDA_wf,'color','#FC4E2A');
hold on
end 
set(gca,'XColor','white') 
set(gca,'YColor','white')
ylim([-1500 1500]);
text(tmin,tmax,'Non-pDA','color','#FC4E2A')
hold on
x = 24:30;
y = -700 *ones(1,length(x));
plot(x,y,'k')
y = -700:-200;
x = 30 *ones(1,length(y));
plot(x,y,'color','k')
text(30.5,-500,'500 v','FontSize', 6)
text (27,-800,'0.25 ms','FontSize', 6)
set (gcf,'PaperPosition',[-1,10,15,35],'PaperSize',[25 15])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' 'DA sorting3']);
npDA_m = [];
npDA_mu = [];
npDA_std = [];
for k = 1:32 
for n = 1: length(npDA_idx)
npDA_wf = FIGURE{npDA_idx(n),4};
npDA_m = [npDA_m npDA_wf(k)];
plot(npDA_wf,'color','r');
hold on
end 
npDA_mu = [npDA_mu mean(npDA_m)];
npDA_std = [npDA_std std(npDA_m)];
end 
figure;
x = 1:32;
std_dev = 1;
curve1 = npDA_mu + std_dev;
curve2 = npDA_mu - std_dev;
fill([x fliplr(x)], [curve1 fliplr(curve2)], [.9 .9 .9], 'linestyle', 'none')
hold all
plot(x,npDA_mu)
ylim([-1500 1500]);
%% 3.2 ISIs 
tDA_idx = find(yTrain ==1);
pDA_idx = find(idx ==2);
npDA_idx = find(idx ==1);
tDA_idx = random(tDA_idx);
pDA_idx = random(pDA_idx);
npDA_idx = random(npDA_idx);
tDA_ISI = [];
pDA_ISI =[];
npDA_ISI = [];
for n = 1:length(tDA_idx)
   % n=30
n=31;
last = max(find(FIGURE{tDA_idx(n),11}<300000)); % first 5 min
spiketime = FIGURE{tDA_idx(n),11}(1:last,:);
spiketime = spiketime';
Spikes = spiketime/1000;
ISI_plotting (Spikes,'#0D4A70','#3C93C2','#9EC9E2'); % burst tonic pasue
end 
for n = 1:length(pDA_idx)
    n=102
last = max(find(FIGURE{pDA_idx(n),11}<300000)); % first 5 min
spiketime = FIGURE{pDA_idx(n),11}(1:last,:);
spiketime = spiketime';
Spikes = spiketime/1000;
ISI_plotting (Spikes,'#06592A','#40AD5A','#6CBA7D'); % burst tonic pasue
end 
for n = 1:length(npDA_idx)
    %n=3;
    n=24
last = max(find(FIGURE{npDA_idx(n),11}<300000)); % first 5 min
spiketime = FIGURE{npDA_idx(n),11}(1:last,:);
spiketime = spiketime';
Spikes = spiketime/1000;
ISI_plotting (Spikes,'#B10026','#FC4E2A','#FEB24C'); % burst tonic pasue
end 
%% 3.3 stable test of different clustering methods
TPpDAs_5model = [];
pDAs_5model = [];
% range = 19:50:489;
range = 50:50:300; % CORRECTED: Was a scalar, now a vector for plotting
k = 2;
%youfun = 2; % 1:kmeans; 2:GMM
for youfun = 4
    TPpDAs_total = [];
    pDAs_total = [];
for n  = 1:100
[TPpDAs,pDAs] = test_stable(table_raw,range,k,youfun);
TPpDAs_total = [TPpDAs_total; TPpDAs];
pDAs_total = [pDAs_total; pDAs];
end
TPpDAs_5model = [TPpDAs_5model TPpDAs_total];
pDAs_5model = [pDAs_5model pDAs_total];
end 
figure ();
subplot (1,2,1)
boxplot(TPpDAs_5model,'Notch','on','Widths',0.7,'Labels',{'k-means','k-medoids','DBSCAN','GMM','OCSVM'})
ylabel('TP(%)')
subplot (1,2,2)
boxplot(pDAs_5model,'Notch','on','Widths',0.7,'Labels',{'k-means','k-medoids','DBSCAN','GMM','OCSVM'})
ylabel('pDA(%)')
saveas(gca,['D:\graduation thesis\results\PCA\58 feature\K-means\''five models'],'svg')
set (gcf,'PaperPosition',[-1,10,20 10],'PaperSize',[20 10])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\''five models']);
std(TPpDAs_5model)
std(pDAs_5model)
figure();
subplot(2,1,1)
errorbar(range, mean(TPpDAs_total), std(TPpDAs_total))
ylim([-10 130]);
xlabel('Number of cells')
title('TPpDAs(%)')
subplot(2,1,2)
errorbar(range, mean(pDAs_total), std(pDAs_total))
ylim([-10 130]);
xlabel('Number of cells')
title('pDAs(%)')
if youfun ==1
    name = 'k-means';
else
     name = 'Gaussian Mixture Model';
end
suptitle(name)
set (gcf,'PaperPosition',[-1,10,15,15],'PaperSize',[20 20])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' name]);
%% 3.4 GMM test for more k
figure
gscatter(Y(:,1),Y(:,2),idx)
hold on 
plot(Y(find(yTrain == 1),1),Y(find(yTrain == 1),2),'ko','MarkerSize',10,'Color', 'g');
hold on 
%plot(Y(find(yTrain == -1),1),Y(find(yTrain == -1),2),'ko','MarkerSize',10,'Color', 'b');
legend ('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','tDA')
xlabel('t-SNE Dim.1'); ylabel('t-SNE Dim.2');
axis( [-20 20 -20 30] ) 
title ('GMM with 5 clusters')
%% 4 Optional: Distribution of repolarization time and ttp 
data_ttp = [];
data_rt = [];
premarks = [];
truemarks = [];
tds = [];
hws = [];
bws = [];
pws = [];
frs = [];
for n = 1: num_train
ptt =  FIGURE{random(n),13};
rt = FIGURE{random(n),7};
td = FIGURE{random(n),9};
hw = FIGURE{random(n),10}; % ms
last = max(find(FIGURE{random(n),11}<300000)); % first 5 min
fr = last/300;% Firing rate(sp/s)
premark = idx(n);
truemark = FIGURE{random(n),5};
data_ttp = [data_ttp ptt];
data_rt=[data_rt rt];
hws = [hws hw];
frs = [frs fr];
tds = [tds td];
premarks = [premarks premark]; % 2 is DA
truemarks = [truemarks truemark]; % 1 is DA
end 
data_ttp = data_ttp';
data_rt = data_rt';
hws = hws';
frs = frs';
% keySet = 1:2; 
% valSet = {'Non-pDA','pDA'}; 
% map = containers.Map(keySet,valSet); 
% premarks = map.values(num2cell(premarks));
% premarks = premarks';
tds = tds';
premarks = premarks';
truemarks = truemarks';
table = [data_ttp data_rt hws frs tds premarks truemarks];
T = array2table(table,'VariableNames',{'Peak-to-trough duration','Repolarization time','Spike width','Firing rate','tds', 'premarks','truemarks'});
figure();
subplot (2,2,1)
s = scatterhistogram(T,'Spike width','Firing rate','GroupVariable','premarks','NumBins',15,'LineWidth',1.2,'ScatterPlotLocation','NorthEast','MarkerSize', 5);
ylabel('Firing rate (sp/s)');xlabel('Spike width (1/3k s)');
s.Color = {'#40AD5A','#FC4E2A'};
s.LegendVisible = 'off';
subplot (2,2,2)
s = scatterhistogram(T,'Peak-to-trough duration','Firing rate','GroupVariable','premarks','NumBins',15,'LineWidth',1.2,'ScatterPlotLocation','NorthEast','MarkerSize', 5);
ylabel('Firing rate (sp/s)');xlabel('Peak-to-trough duration (1/3k s)');
s.Color = {'#40AD5A','#FC4E2A'};
s.LegendVisible = 'off';
subplot (2,2,3)
s = scatterhistogram(T,'Repolarization time','Firing rate','GroupVariable','premarks','NumBins',15,'LineWidth',1.2,'ScatterPlotLocation','NorthEast','MarkerSize', 5);
ylabel('Firing rate (sp/s)');xlabel('Repolarization time (1/3k s)');
s.Color = {'#40AD5A','#FC4E2A'};
s.LegendVisible = 'off';
subplot (2,2,4)
s = scatterhistogram(T,'tds','Firing rate','GroupVariable','premarks','NumBins',15,'LineWidth',1.2,'ScatterPlotLocation','NorthEast','MarkerSize', 5);
ylabel('Firing rate (sp/s)');xlabel('Total duration 2 (1/3k s)');
s.Color = {'#40AD5A','#FC4E2A'};
s.LegendVisible = 'off';
set (gcf,'PaperPosition',[-1,10,25 20],'PaperSize',[25 20])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' 'FR vs waveform']);
% set (gcf,'PaperPosition',[-1,10,25 10],'PaperSize',[25 10])
% print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' 'FR vs waveform']);
% %%
% csvwrite('D:\graduation thesis\results\PCA\58 feature\table_ztrain_pca.csv',table_ztrain_pca)
