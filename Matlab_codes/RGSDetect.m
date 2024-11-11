function [Bursts, Pauses] = RGSDetect(Spikes, N_min, Steps, p, alpha,plotting)
%% RGSDETECT
%
% Determines burst and pause interspike interval (ISI) thresholds and
% identifies burst and pause strings based on the Robust Gaussian Surprise
% (RGS) method.
%
% INPUTS:
% Spikes (s): Times of spikes in seconds.
% N_min: Minimum number of spikes to be considered a burst/pause
% string.
% Steps (log10(s)): Bin edges for histogram count (histc) of the
% log ISIs.
% p: Bottom and top p% used as outliers to calculate central
% location; keep p in range [0.05, 0.30] (default 0.05).
% alpha: Value used in Bonferroni correction; lower value of
% alpha to filter out false positives (default 0.05).
% NOTE: Requires MATLAB statistics toolbox.
%
% OUTPUTS:
% Bursts: Structure containing burst information.
% Bursts.BurstingSpikes (s): Column of times of all spike times
% included in a burst. before Bonferroni correction! 
% Bursts.IBF (Hz): Column of intraburst frequency (IBF) of each
% burst.
% Bursts.NumSpikes: Column of number of spikes in each burst.
% Bursts.Windows (s): 2 Columns of start and end times of each
% burst.
% Pauses: Structure containing pause information.
% Pauses.AllSpikes (s): Start times of all ISIs that satisfy
% pause threshold.
% Pauses.AllLengths (s): Lengths of all ISIs that satisfy pause
% threshold.
% Pauses.PausingSpikes (s): Column of all spike times
% included in a pause string.
% Pauses.IPF (Hz): Column of intrapause frequency (IPF) of each
% pause string.
% Pauses.NumSpikes: Column of number of spikes in each pause
% string.
% Pauses.Windows (s): 2 Columns of start and end times of each
% pause string.
% NOTE: Rows of structure elements correspond to the same burst or
% pause.
% NOTE: Normalized Log ISI Distribution (NLISI) plot is used confirm
% central distribution is centered on 0. If distribution is not
% centered on 0, change p until it is. Use steps to adjust the
% x-axis.
%
% EXAMPLE values
% Spikes = (Spike times go here);
% N_min = 2;
% Steps = -3:0.005:1.5;
% p = 0.05;
% alpha = 0.05;
% [Bursts, Pauses] = RGSDetect(Spikes, N_min, Steps, p, alpha);
%
% REFERENCE: Ko et al. (2012)
% Detection of bursts and pauses in spike trains
% J Neurosci Methods 211:145¨C158
%
%%% NORMALIZED LOG ISI DISTRIBUTION
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
%For the middle portion
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

%Plot the NLISI Distribution
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

% --------------- important plotting ----------
if plotting == 1
% plot histograms 
figure(); 
% NLISITrain_p = NLISITrain./sum(NLISITrain);
hist(NLISITrain,[min(NLISITrain):0.025:max(NLISITrain)]);
hold on; 
plot([CentralDistBounds(1) CentralDistBounds(1)], get(gca, 'YLim'), '-r');
plot([CentralDistBounds(2) CentralDistBounds(2)], get(gca, 'YLim'), '-b');
xlabel ('Normalized Log ISI');
ylabel ('Number of spikes');
legend ('Normalized Log ISI','burst threshold','pause threshold')
else 
end 
% -----------------------------------------

% figure(); 
% % NLISITrain_p = NLISITrain./sum(NLISITrain);
% bin_edges = min(NLISITrain):0.025:max(NLISITrain);
% mask = NLISITrain< CentralDistBounds(1);
% hist(NLISITrain(mask),bin_edges);
% hold on
% histogram(NLISITrain(~mask),bin_edges,'FaceColor','red');
% hold off
% 
% hist(NLISITrain,bin_edges);
% hold on; 
% plot([CentralDistBounds(1) CentralDistBounds(1)], get(gca, 'YLim'), '-r');
% plot([CentralDistBounds(2) CentralDistBounds(2)], get(gca, 'YLim'), '-b');
% xlabel ('Normalized Log ISI');
% ylabel ('Number of spikes');

%% BURST detection
if sigma == inf
   Bursts = nan;
else 
%Get index and ISI lengths of all ISIs that satisfy burst threshold
Burst_Thresh = CentralDistBounds(1);
BurstINDXS = 1:length(NLISITrain);
BurstINDXS(NLISITrain >= Burst_Thresh) = []; %Delete all indexes greater than
%the burst threshold
if ~isempty(BurstINDXS)
%Matrix of all potential burst ISIs and their indexes
BurstsM = [NLISITrain(NLISITrain < Burst_Thresh)';BurstINDXS];
[~,c] = size(BurstsM);
Burst_Seed = mat2cell(BurstsM,2,ones(1,c,1));
else
Burst_Seed = {};
Bursts = Burst_Seed;
end

%Loop through each potential burst ISI (Burst Seed)
for i = 1:length(Burst_Seed);
%Go forward and backward from the current burst until both conditions
%are unsatisfied
forward = 1;
backward = 1;
while forward || backward
currBurst = cell2mat(Burst_Seed(i));
%Go forward 1 ISI
if forward
%Set current ISI as end of the current burst
currSpike = currBurst(:,end);
if currSpike(2) ~= length(NLISITrain)
%q is number of spikes
[~,q] = size(currBurst);
%P1 is probability burst will occur assuming Gaussian
%distribution with mean, mu¡Áq, and std, sqrt(q)¡Ásigma
P1 = normcdf(sum(currBurst(1,:)), mu*q, sqrt(q).*sigma);
testBurst = [currBurst [NLISITrain(currSpike(2)+1);currSpike(2)+1]];
%P2 is the same probability with the next ISI added to the
%burst
P2 = normcdf(sum(testBurst(1,:)), mu*(q+1),sqrt(q+1).*sigma);
%If the next ISI increased the probability of the burst
%occurring
if P2 >= P1
%Stop going forward
forward = 0;
else
%Otherwise, set the current burst seed to the tested
%burst
Burst_Seed{i} = testBurst;
end
else
%Stop going forward if at the end of the ISI train
forward = 0;
end
end
currBurst = cell2mat(Burst_Seed(i));
%Go backward 1 ISI
if backward
%Set current ISI as end of the current burst
currSpike = currBurst(:,1);
if currSpike(2) ~= 1
%q is number of spikes
[~,q] = size(currBurst);
%P1 is probability burst will occur assuming Gaussian
%distribution with mean, mu¡Áq, and std, sqrt(q)¡Ásigma
P1 = normcdf(sum(currBurst(1,:)), mu*q, sqrt(q).*sigma);
testBurst = [[NLISITrain(currSpike(2)-1);currSpike(2)-1] currBurst];
%P2 is the same probability with the next ISI added to the
%burst
P2 = normcdf(sum(testBurst(1,:)), mu*(q+1), sqrt(q+1).*sigma);
%If the next ISI increased the probability of the burst
%occurring
if P2 >= P1
%Stop going backward
backward = 0;
else
%Otherwise, set the current burst seed to the tested
%burst
Burst_Seed{i} = testBurst;
end
else
%Stop going backward if at the end of the ISI train
backward = 0;
end
end
end
end

if ~isempty(Burst_Seed)
%Initialize BurstInfo
BurstInfo = zeros(length(Burst_Seed),3);
%Get start index of each burst
BurstInfo(:,1) = cellfun(@(x) x(2,1),Burst_Seed);
%Get end index of each burst
BurstInfo(:,2) = cellfun(@(x) x(2,end),Burst_Seed);
%Get P-value of each burst (probability of occurence assuming Gaussian 
%distribution)
BurstInfo(:,3) = cellfun(@(x) normcdf(sum(x(1,:)), mu*length(x),sqrt(length(x)).*sigma),Burst_Seed);
%Filter out bursts less than minimum number of spikes specified by N_min
BurstInfo(BurstInfo(:,2)-BurstInfo(:,1)+2 < N_min,:) = [];
%Filter out overlapping bursts
no_overlap = 0;
i=1;
if ~isempty(BurstInfo)
[r,~] = size(BurstInfo);
if r ~= 1
while ~no_overlap
%If the indexes of the burst ISIs don't intersect
if isempty(intersect(BurstInfo(i,1):BurstInfo(i,2),BurstInfo(i+1,1):BurstInfo(i+1,2)))
%move to the next burst
i = i+1;
else
%If the intersect, choose the burst with the lower P
%value
if BurstInfo(i,3) <= BurstInfo(i+1,3)
BurstInfo(i+1,:) = [];
else
BurstInfo(i,:) = [];
end
end
%When the end is reached, stop
[r,~] = size(BurstInfo);
if i == r
no_overlap = 1;
end
end
end
%r is the number of rows or the number of bursts
[r,~] = size(BurstInfo);
Bursts.BurstingSpikes = {};
% for each burst, append the burst spikes
% for i = 1:r
% Bursts.BurstingSpikes = [Bursts.BurstingSpikes ;Spikes(BurstInfo(i,1):BurstInfo(i,2)+1)]; % important
% end

for i = 1:r
Bursts.BurstingSpikes{i,1} = Spikes(BurstInfo(i,1):BurstInfo(i,2)+1); % important
end

end
%Bonferroni correction for false positives
KB = length(find(BurstInfo(:,3) < alpha));
BurstInfo(BurstInfo(:,3)*KB >= alpha,:) = [];
%Use the indexes in burst info to find the burst windows
Bursts.Windows1 = [];
Bursts.Windows2 = [];
Bursts.Windows1 = [Bursts.Windows1 Spikes(BurstInfo(:,1))];
Bursts.Windows2 = [Bursts.Windows2 Spikes(BurstInfo(:,2)+1)];
Bursts.Windows = [Bursts.Windows1;Bursts.Windows2]';
%Use the indexes to find the number of spikes in each burst
Bursts.NumSpikes = BurstInfo(:,2) - BurstInfo(:,1) + 2;
%Use the number of spikes and windows to calculate the IBF
Bursts.IBF = Bursts.NumSpikes./(Bursts.Windows(:,2) - Bursts.Windows(:,1));
end
end 
%% Pause
if sigma == inf
   Pauses = nan;
else 
%Get index and ISI lengths of all NLISIs that satisfy pause threshold
Pause_Thresh = CentralDistBounds(2);
PauseINDXS = 1:length(NLISITrain);
%Delete all indexes less than the pause threshold
PauseINDXS(NLISITrain <= Pause_Thresh) = [];

if ~isempty(PauseINDXS)
%Matrix of all potential pause string NLISIs and their indexes
PausesM = [NLISITrain(NLISITrain > Pause_Thresh)';PauseINDXS];
[~,c] = size(PausesM);
Pause_Seed = mat2cell(PausesM,2,ones(1,c,1));
else
Pause_Seed = {};
Pauses = [];
end

%Loop through each potential pause string NLISI (Pause Seed)
for i = 1:length(Pause_Seed);
%Go forward and backward from the current pause string until both conditions
%are unsatisfied
forward = 1;
backward = 1;
while forward || backward
currPause = cell2mat(Pause_Seed(i));
%Go forward 1 ISI
if forward
%Set current ISI as end of the current pause string
currPauseind = currPause(:,end);
if currPauseind(2) ~= length(NLISITrain)
[~,q] = size(currPause);
%P1 is probability pause string will occur assuming Gaussian
%distribution with mean, mu¡Áq, and std, sqrt(q)¡Ásigma
P1 = (1-normcdf(sum(currPause(1,:)), mu*q, sqrt(q).*sigma));
testPause = [currPause [NLISITrain(currPauseind(2)+1);currPauseind(2)+1]]; % important
%P2 is the same probability with the next ISI added to the
%pause string
P2 = (1-normcdf(sum(testPause(1,:)), mu*(q+1), sqrt(q+1).*sigma));
%If the next ISI increased the probability of the pause
%string occurring
if P2 >= P1
%Stop going forward
forward = 0;
else
%Otherwise, set the current pause seed to the tested
%pause string
Pause_Seed{i} = testPause;
end
else
forward = 0;
end
end
currPause = cell2mat(Pause_Seed(i));
%Go backward 1 ISI
if backward
currPauseind = currPause(:,1);
if currPauseind(2) ~= 1
[~,q] = size(currPause);
%P1 is probability pause string will occur assuming Gaussian
%distribution with mean, mu¡Áq, and std, sqrt(q)¡Ásigma
P1 = (1-normcdf(sum(currPause(1,:)), mu*q, sqrt(q).*sigma));
testPause = [[NLISITrain(currPauseind(2)-1);currPauseind(2)-1] currPause];
%P2 is the same probability with the next ISI added to the
%pause string
P2 = (1-normcdf(sum(currPause(1,:)), mu*(q+1), sqrt(q+1).*sigma));
%If the next ISI increased the probability of the pause
%string occurring
if P2 >= P1
%Stop going forward
backward = 0;
else
%Otherwise, set the current pause seed to the tested
%pause string
Pause_Seed{i} = testPause;
end
else
backward = 0;
end
end
end
end

if ~isempty(Pause_Seed)
%Initialize PauseInfo variable
PauseInfo = zeros(length(Pause_Seed),3);
%Starting indexes of pause strings
PauseInfo(:,1) = cellfun(@(x) x(2,1),Pause_Seed);
%Ending indexes of pause strings
PauseInfo(:,2) = cellfun(@(x) x(2,end),Pause_Seed);
%P-value of the pause strings
PauseInfo(:,3) = cellfun(@(x) (1-normcdf(sum(x(1,:)), mu*length(x), sqrt(length(x)).*sigma)),Pause_Seed);
%Minimum number of spikes filter
PauseInfo(PauseInfo(:,2)-PauseInfo(:,1) + 2 < N_min,:) = [];
%Filter out overlaps
no_overlap = 0;
i=1;
if ~isempty(PauseInfo)
%r is number of current pause strings
[r,~] = size(PauseInfo);
if r ~= 1
while ~no_overlap
%If the indexes of the burst ISIs don't intersect
if isempty(intersect(PauseInfo(i,1):PauseInfo(i,2),PauseInfo(i+1,1):PauseInfo(i+1,2)))
%Move to next pause string
i = i+1;
else
%Choose the pause string with lower P-value
if PauseInfo(i,3) <= PauseInfo(i+1,3)
PauseInfo(i+1,:) = [];
else
PauseInfo(i,:) = [];
end
end
%End if the last pause string is reached
[r,~] = size(PauseInfo);
if i == r
no_overlap = 1;
end
end
end
end
%Use indexes to find start and end times
Pauses.Windows1 = [];
Pauses.Windows2 = [];
Pauses.Windows1 = [Pauses.Windows1 Spikes(PauseInfo(:,1))];
Pauses.Windows2 = [Pauses.Windows2 Spikes(PauseInfo(:,2)+1)];
Pauses.Windows = [Pauses.Windows1;Pauses.Windows2]';
%Use indexes to determine number of spikes
Pauses.NumSpikes = PauseInfo(:,2) - PauseInfo(:,1) + 2;
%Use windows and numspikes to calculate IPF
Pauses.IPF = Pauses.NumSpikes./(Pauses.Windows(:,2) - Pauses.Windows(:,1));
Pauses.PausingSpikes = {};
%Add pausing spikes from each pause string to the pausingspikes element
[r,~] = size(PauseInfo);
for i = 1:r
Pauses.PausingSpikes{i,1} = Spikes(PauseInfo(i,1):PauseInfo(i,2)+1); % important
end
%Bonferroni correction
KP = length(find(PauseInfo(:,3) < alpha));
PauseInfo(PauseInfo(:,3)*KP >= alpha,:) = [];
%Get all pauses using the pause indexes
Pauses.AllSpikes = Spikes(PauseINDXS);
%Get the lengths of all the pauses that satisfy the threshold
Pauses.AllLengths = ISIs(PauseINDXS);
end

end 
end
