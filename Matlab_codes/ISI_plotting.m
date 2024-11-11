function ISI_plotting (Spikes,burstc,tonicc,pausec)
Steps = -3:0.005:1.5;
p = 0.05;
alpha = 0.05;
ISIs = Spikes(2:end)-Spikes(1:end-1); %Offset spikes by 1 and subtract
%for ISI's
LogISIs = log10(ISIs); %Take the log10 of ISI's
N = length(LogISIs); %N = number of ISIs
Q = max([20,floor(0.2*N)]); %Set window length as max of 20 and 20% of N
NLISITrain = zeros(1,length(LogISIs))'; %Initialize Normalized Log ISI Train
%Central Location of first 2¡ÁQ+1 ISIs
CentralLoc1 = ComputeCL(LogISIs(1:(2*Q+1)), Steps, p);
%Subtract Central Location (Normalize)
NLISITrain(1:Q) = LogISIs(1:Q) - CentralLoc1;
%Central Location of last 2¡ÁQ+1 ISIs
CentralLoc2 = ComputeCL(LogISIs((N-2*Q):N), Steps, p);
%Subtract Central Location (Normalize)
NLISITrain(N-Q+1:end) = LogISIs(N-Q+1:end) - CentralLoc2;
for i = Q+1:N-Q
%Compute central location for portion of index +/? Q and subtract
NLISITrain(i) = LogISIs(i) - ComputeCL(LogISIs(i-Q:i+Q), Steps, p);
end
%Get statistics of the NLISI train
med = median(NLISITrain);
pool_MAD = mad(NLISITrain);
CentralDistBounds = [med - pool_MAD*2.58 med + pool_MAD*2.58];
mu = median(NLISITrain);
sigma = mad(NLISITrain);

ISIidx = ones (1,length(NLISITrain));
Pause_ISIidx =find (NLISITrain>CentralDistBounds(2));
Burst_ISIidx =find (NLISITrain<CentralDistBounds(1));
ISIidx(Pause_ISIidx) = 3;
ISIidx(Burst_ISIidx) = 2;

% % plot histograms 
% figure(); 
% % NLISITrain_p = NLISITrain./sum(NLISITrain);
% if min(NLISITrain) == -inf
% hist(NLISITrain,[-3:0.025:max(NLISITrain)]);
% hold on; 
% else
% hist(NLISITrain,[min(NLISITrain):0.025:max(NLISITrain)]);
% hold on; 
% end 
% plot([CentralDistBounds(1) CentralDistBounds(1)], get(gca, 'YLim'), '-r');
% plot([CentralDistBounds(2) CentralDistBounds(2)], get(gca, 'YLim'), '-b');
% xlabel ('Normalized Log ISI');
% ylabel ('Number of spikes');
% legend ('Normalized Log ISI','burst threshold','pause threshold')
% plot histograms with different colors

% figure(); 
% if min(NLISITrain) == -inf
% histogram(NLISITrain(find(ISIidx == 3)),[-3:0.025:max(NLISITrain)],'FaceColor',pausec,'EdgeColor','none');
% hold on; 
% histogram(NLISITrain(find(ISIidx == 1)),[-3:0.025:max(NLISITrain)],'FaceColor',tonicc,'EdgeColor','none');
% hold on 
% histogram(NLISITrain(find(ISIidx == 2)),[-3:0.025:max(NLISITrain)],'FaceColor',burstc,'EdgeColor','none');
% hold on
% else
% histogram(NLISITrain(find(ISIidx == 3)),[min(NLISITrain):0.025:max(NLISITrain)],'FaceColor',pausec,'EdgeColor','none');
% hold on; 
% histogram(NLISITrain(find(ISIidx == 1)),[min(NLISITrain):0.025:max(NLISITrain)],'FaceColor',tonicc,'EdgeColor','none');
% hold on 
% histogram(NLISITrain(find(ISIidx == 2)),[min(NLISITrain):0.025:max(NLISITrain)],'FaceColor',burstc,'EdgeColor','none');
% hold on
% end 
% ylim([0 400])  
% plot([CentralDistBounds(1) CentralDistBounds(1)], get(gca, 'YLim'), '-k','LineWidth',1);
% plot([CentralDistBounds(2) CentralDistBounds(2)], get(gca, 'YLim'), '-k','LineWidth',1);
% xlabel ('Normalized Log ISI');
% ylabel ('Number of spikes');
% legend ('Pause','Burst','Tonic','Thresholds','location' ,'NorthWest'  )

% 
% % Plot the NLISI Distribution 
% figure
% hold on
% %Run a smoothing pdf kernel.
% NLISIpdf = pdf(fitdist(NLISITrain,'Kernel'), Steps);
% NLISIpdf = NLISIpdf./sum(NLISIpdf);
% plot(Steps, NLISIpdf,'g')
% %Plot treshold lines
% plot([CentralDistBounds(1) CentralDistBounds(1)], [0 max(NLISIpdf)], '-r')
% plot([CentralDistBounds(2) CentralDistBounds(2)], [0 max(NLISIpdf)], '-b')
% xlabel 'Normalized Log ISIs'
% ylabel 'Probability'
% title 'Normalized Log ISI Distribution'
minx = -2.5; 
maxx = 2.5;
figure(); 
if min(NLISITrain) == -inf
histogram(NLISITrain(find(ISIidx == 3)),[minx:0.025:maxx],'FaceColor',pausec,'EdgeColor','none');
hold on; 
histogram(NLISITrain(find(ISIidx == 1)),[minx:0.025:maxx],'FaceColor',tonicc,'EdgeColor','none');
hold on 
histogram(NLISITrain(find(ISIidx == 2)),[minx:0.025:maxx],'FaceColor',burstc,'EdgeColor','none');
hold on
else
histogram(NLISITrain(find(ISIidx == 3)),[minx:0.025:maxx],'FaceColor',pausec,'EdgeColor','none');
hold on; 
histogram(NLISITrain(find(ISIidx == 1)),[minx:0.025:maxx],'FaceColor',tonicc,'EdgeColor','none');
hold on 
histogram(NLISITrain(find(ISIidx == 2)),[minx:0.025:maxx],'FaceColor',burstc,'EdgeColor','none');
hold on
end 
ylim([0 350])  
%ylim([0 100])  
plot([CentralDistBounds(1) CentralDistBounds(1)], get(gca, 'YLim'), '-k','LineWidth',1);
plot([CentralDistBounds(2) CentralDistBounds(2)], get(gca, 'YLim'), '-k','LineWidth',1);
xlabel ('Normalized Log ISI');
ylabel ('Count');
legend ('Pause','Tonic','Burst','Thresholds','location' ,'NorthWest'  )
set (gcf,'PaperPosition',[-1,10,15,5],'PaperSize',[10 15])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\PCA\58 feature\K-means\' burstc]);
end 