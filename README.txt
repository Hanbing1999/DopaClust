Unsupervised clustering for putative dopamine neuron identification from tetrode recordings	

***********************************************************************************
%Scripts and functions

King_2.m : Main script, call for other functions to conduct feature extractions on raw dataset

RGSDetect.m : A function to detect bursts and pauses 

ComputeCL.m : A function within RGSDetect.m

plot_ptt_new.m : A function to compute peak-to-trough and repolarization time 

plot_ptt_2.m: A function to plot peak-to-trough and repolarization time

plot_td_2.m: A function to compute the total duration

Rasterplot.m : A function to plot all spikes along time scale

Clustering.m: Main script, to normalize features, reduce dimensionality, and unsupervised clustering

Clustering_plotting.m: A function to plot the clustering results

feature_table.m: A function to form 38 feature-table

test_stable.m: test the stability of different clustering methods
***********************************************************************************
HOW TO USE THESE FUNCTIONS?
1. Open King_2.m to preprocess dataset and extract features from electrode recordings 
2. Open Clustering.m to test different clustering methods

**************************************************************************************
% FIGURE: main cell included all information
1. Cell names 
2. DAN.unref_avg_WF
3. DAN.best_unref_channel: name of the best channel 
4. signals from the best channel
5. photoidentified or no-photoidentified
6. the order of the cell
7. repolarizatime time.
8. total duration_1
9. total duration_2 (80% amplitude)
10. half-width
11. Input time scale in FIGURE
12.
13. ptt: peak to trough duration 
14. Bursts
15. Pauses
16.Bursts_features = [mean(Bursts_matrix);median(Bursts_matrix);std(Bursts_matrix);
    skewness(Bursts_matrix);kurtosis(Bursts_matrix)]; &&& matrix == [Bursts_windows(s) Bursts_IBFs (HZ) Bursts_NumSpikes];
17.per_Bursts_windows
18.Pauses_features = [mean(Pauses_matrix);median(Pauses_matrix);std(Pauses_matrix);
    skewness(Pauses_matrix);kurtosis(Pauses_matrix)];
19.per_Pauses_windows
20. Bursts 
21. Pauses
22. Bursts
23. Pauses
24. Tonics
25. [per_Bursts_windows Bursts_ISIfeature' Bursts_spikefeature']
26.[per_Pauses_windows Pauses_ISIfeature' Pauses_spikefeature']
27.  per_Tonics_windows
***************************************************************************************
*********************************************************
table_raw.m: A cell to save feature table, training sample, and testing sample
{1,1}:37 features
{2,1}:58 features
****************************************************************************************888
Feature table (58 features + 1 mark[1: photo-tagged,-1: non-identified ]): 
1. data_ttp % trough to peak duration (ms)
2. data_rt % repolarization time (ms)
3. td1s % total duration (ms)
4. td2s 80% amplitude total duration (ms)
5. hws  % half-width (ms)
6-31. Burst_features: 26
32-57. Pause_features: 26
58:Tonic_features: 1
59: Mark

********************
to do list 
2. IPF
