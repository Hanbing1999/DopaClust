function [CentralLocation] = ComputeCL(ISIs, Steps, p)
%% COMPUTECL help
%
% Subroutine required for MATLAB code for RGS burst and pause detection
% (Table 4). This subroutine computes the central location given an ISI train
% using robust measures
% of the median absolute difference (MAD), median, and central set.
%
% INPUTS:
% ISIs (s): Lengths of ISIs in seconds.
% Steps: Bin edges for histogram count (histc) of the ISIs.
% p: Bottom and top p% used as outliers to calculate central
% location; keep p in range [0.05, 0.30] (default 0.05).
% NOTE: RGSDetect inputs log scale ISIs and Steps.
% OUTPUTS:
% CentralLocation: Central location of the ISI distribution
%
% REFERENCE: Ko et al. (2012) Detection of bursts and
% pauses in spike trains. J Neurosci Methods 211:145¨C158
%
%% COMPUTECL
%
%Locates the bin centers of the steps on a linear scale
bincenters = [(Steps(1:end-1) + Steps(2:end))./2 Steps(end)];
%Histogram counts the ISIs using Steps
ISIhist = histc (ISIs,Steps)';
%Converts to probability distribution
normhist = ISIhist./sum(ISIhist);
%Creates cumulative probability distribution
cumprob = cumsum(normhist);
%Calculates thresholds for bottom and top p quantiles
[~,burstquantid] = min(abs(cumprob-(p)));
[~,pausequantid] = min(abs(cumprob-(1-p)));
burstquant = bincenters(burstquantid);
pausequant = bincenters(pausequantid);
%Caclulates E-Center as average of 2 thresholds
E_Center = (burstquant + pausequant)/2;
%Calculates central set using MAD
CentralSetBoundaries = [E_Center - 1.64*mad(ISIs,1) E_Center + 1.64*mad(ISIs,1)];

CentralSet = [bincenters(bincenters >=CentralSetBoundaries(1) &bincenters <=CentralSetBoundaries(2)); (normhist(bincenters >=CentralSetBoundaries(1) & bincenters <=CentralSetBoundaries(2)))'];

%Calculates median of central set
CentralDistCumProb = cumsum(CentralSet(2,:)./sum(CentralSet(2,:)));
[~,CentralLocationid] = min(abs(CentralDistCumProb - 0.5));
CentralLocation = CentralSet(1,CentralLocationid);
end