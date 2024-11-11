function [TPpDA,pDA]= Clustering_plotting(x,m,yTrain,idx,name,n)
rng default % for reproducibility
Y = tsne(x,'Algorithm','barneshut','NumPCAComponents',m,'NumDimensions',3);
scatter3(Y(find(yTrain == 1),1),Y(find(yTrain == 1),2),Y(find(yTrain == 1),3),10,yTrain(find(yTrain == 1)),'filled');
hold on 
scatter3(Y(find(yTrain == -1),1),Y(find(yTrain == -1),2),Y(find(yTrain == -1),3),10,yTrain(find(yTrain == -1)),'filled');
hold on 
plot3(Y(find(idx==1),1),Y(find(idx== 1),2),Y(find(idx== 1),3),'ko','MarkerSize',5,'Color', 'r')
plot3(Y(find(idx==2),1),Y(find(idx==2),2),Y(find(idx== 2),3),'ko','MarkerSize',5,'Color', 'g')
h1 = legend('DA phototagged','Non-phototagged','Non-putative DA','Putative DA');
set(h1,'Box','off');
xlabel('t-SNE Dim.1'); ylabel('t-SNE Dim.2');zlabel('t-SNE Dim.3')
axis( [-20 20 -20 30 -10 30] )   
% title(['t-SNE after ' name 'clustering (with PCA)']) 
title(name,'FontSize',13) 
view(-47,28)
% the percentage of true DA neurons assigned to pDA group
TPpDA = length(intersect(find(yTrain == 1),find(idx==2)))/length(find(yTrain == 1))*100; 
% total propotion of pDA neurons among all neurons
pDA = length(find(idx==2))/length(yTrain) *100;
text(-20,20,35,['TP:' num2str(roundn(TPpDA,-2)) '%'])
text(-20,20,30,['pDA:' num2str(roundn(pDA,-2)) '%'])
% set (gcf,'PaperPosition',[-1,10,15,20],'PaperSize',[27 23])
% print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' name]);
% set (gca,'PaperPosition',[-1,10,15,20],'PaperSize',[27 23])
% saveas(gca,['D:\graduation thesis\results\PCA\58 feature\K-means\' name],'tiff')
end

