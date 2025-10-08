% This main script is for feature extraction to conduct putative DA neuron identification
%
% (c) 2020-2021 Hanbing Shen
%% 0. import path and save path
mfilePath = mfilename('fullpath');
% Check if mfilename returned a Live Editor temporary file path
if contains(mfilePath, 'LiveEditorEvaluationHelper')
    % If so, get the path of the currently active file in the editor
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
% Extract just the folder path
current_path = fileparts(mfilePath);
% Extract the parent folder 
path = fileparts(current_path);
path = [path '/'];
total_cell = 489;
save_path = [path 'Results_Hanbingrun/'];
mkdir(save_path);
disp(['Directory created (or already exists): ', save_path]);
%% 0.1. load raw data and check name (Only use when the cells need to be changed)
Data = load ([path 'DA_spikedata/DAN_list.mat']); % import the photo-tagged neuron mat file
DANs = Data.DANs;
FIGURE = {}; % form an empty cell to save information
m=0;
for n= 1:length(DANs)
    name = DANs{n,1};
    if isempty(name) % eliminate empty file
        m = m+1;
    end
    if ~isempty(name)
        DAN = load ([path 'SpikeDataOL_v2/' name '.mat'] ); % read the specific name of one cell
        FIGURE {n-m,1} = name; % save name
        FIGURE {n-m,2} = [DAN.unref_avg_WF];  % save unref_avg_WF
        FIGURE {n-m,3} = [DAN.best_unref_channel];
        FIGURE {n-m,4} = [DAN.unref_avg_WF(:,DAN.best_unref_channel)];
        FIGURE {n-m,6} = n;
    end
end
for n =1:length(FIGURE)
    FIGURE {n,5} = 1; % mark the photo-tagged DA units as 1
end
dir_info = dir( [path 'SpikeDataOL_v2']); % read the name of each cell
dir_info = dir_info(3:491); % eliminate useless lines
s = 122;
for n = 1:total_cell
    name_cell = dir_info(n).name;
    name_cell = name_cell (1:end-4);
    if isempty(find(strcmp(DANs, name_cell))) % if this cell is not in DAN_list, this cell is not photoidentified
        s = s+1;
        cell = load ([path 'SpikeDataOL_v2/' name_cell '.mat'] );
        FIGURE {s,1} = name_cell;
        FIGURE {s,2} = [cell.unref_avg_WF];
        FIGURE {s,3} = [cell.best_unref_channel];
        FIGURE {s,4} = [cell.unref_avg_WF(:,cell.best_unref_channel)];
        FIGURE {s,5} = -1; % mark unidentified units as -1
    end
end
% Input spikes (firing patterns) of each cell
for n = 1:total_cell
    name = FIGURE{n,1};
    d = h5read([path 'timestamps/' name '_public.h5'],'/spike');
    FIGURE {n,11} =d;
end
save ([save_path 'FIGURE.mat'],'FIGURE');
%% 1.1 load or save data for analysis
% load ('D:\graduation thesis\results\FIGURE.mat');
%% 1.2 Plot avg_wv for all units
for m = 1:20
    figure();
    for n = 8*(m-1)+1:8*m
        subplot (2,4,(n-8*(m-1)))
        plot (FIGURE{n,2})
        name = FIGURE{n,1};
        title(name,'fontsize',10);
        xlabel ('Times Step','fontsize',10);
        ylabel ('unref avg WF','fontsize',10);
        h= legend({'C1','C2','C3','C4'});
        set(h,'fontsize',9);
        h;
        text(5,-400,['Best Channel: C ', num2str(FIGURE{n,3})],'fontsize',9);
        ylim([-800 800]);
    end
    %     set (gcf,'Position',[1500,700,1500,800],'PaperSize',[30 20])
    %     F = getframe(gcf);
    %     filename = [save_path ,num2str(m),'.pdf'];
    %     saveas(gcf,filename,'pdf')
end
%% 2.1 Plot and save total duration and half-width
plotting = 1; %0 not plot
ts = [1:32];
n = total_cell; % number of cells
plot_infp ={}; % create an empty cell to save information for plotting
% Call the total duration and half-width function
for cell = 1:n
    [td_sp, td_sp_2, hw, cs, x_half,xStart_sp,xmax_sp,xStart_sp_2,xmax_sp_2] = plot_td_2(cell,ts,FIGURE);
    FIGURE {cell,8} = td_sp;
    FIGURE {cell,9} = td_sp_2;
    FIGURE {cell,10} = hw;
    plot_infp{cell,1} =cs;
    plot_infp{cell,2} =x_half;
    plot_infp{cell,3} =xStart_sp;
    plot_infp{cell,4} =xmax_sp;
    plot_infp{cell,5} =xStart_sp_2;
    plot_infp{cell,6} =xmax_sp_2;
end
%% 2.2 Plotting half width
for m = 1:61
    figure();
    num = 8*(m-1)+1:8*m;
    for h = 1:8
        cell = num(h);
        subplot(2,4,h)
        fnplt (plot_infp{cell,1})
        hold on
        ylim([-600 600])
        %plot (ts, data,'ro');
        plot  (repmat (plot_infp{cell,2},2,1),repmat([-600;600],1,length(plot_infp{cell,2})),':','linewidth',2,'color','#E32977');
        %title ('half width')
        if cell == 1
            h1 = legend('Data','Half width');
            set(h1,'Box','off');
        else
        end
    end
    sgtitle('Half width')
    xlabel('Time step (1/3k s)','position',[-50,-800,4],'fontsize',13)
    hAxis = axes('visible','off');
    h = text(-0.1,0.5,'Amplitude (v)');
    set(h,'fontsize',13,'rotation',90,'HorizontalAlignment','center')
    %     set (gcf,'Position',[1500,700,1500,800],'PaperSize',[30 20])
    %     print(gcf,'-dtiff','-r300',[save_path 'half-width' num2str(m)]);
end
%% 2.3 Plotting: Total duration and local extreme
for m = 1:61
    figure();
    num = 8*(m-1)+1:8*m;
    for h = 1:8
        cell = num(h);
        subplot(2,4,h)
        fnplt (plot_infp{cell,1});
        %plot(ts,data,'x','color','[.667 .667 1]'); % 'color','[.667 .667 1]'
        hold on
        %plot(x,y,'ro')
        hold on
        ylim([-600 600])
        plot([plot_infp{cell,3},plot_infp{cell,3}],ylim,':','linewidth',2,'Color','#40AD5A');
        plot([plot_infp{cell,4},plot_infp{cell,4}],ylim,':','linewidth',2,'Color','#40AD5A');
        plot([plot_infp{cell,5},plot_infp{cell,5}],ylim,':','linewidth',2,'Color','#E32977');
        plot([plot_infp{cell,6},plot_infp{cell,6}],ylim,':','linewidth',2,'Color','#E32977');
        %plot([1 32],[bl,bl],'color','#A0B1BA');
        % legend('data','model','local extreme','xStart','xmax','xStart_2','xmax_2','bl_1','Location','SouthEast');
        if cell == 1
            h1 = legend('Data','Start1','Peak1','Start2','Peak2','Location','SouthEast');
            set(h1,'Box','off');
        else
        end
    end
    suptitle('Total duration')
    xlabel('Time step (1/3k s)','position',[-50,-800,4],'fontsize',13)
    hAxis = axes('visible','off')
    h = text(-0.1,0.5,'Amplitude (v)')
    set(h,'fontsize',13,'rotation',90,'HorizontalAlignment','center')
    %set (gcf,'Position',[1500,700,1500,800],'PaperSize',[30 20])
    %print(gcf,'-dtiff','-r300',[save_path 'total_duration' num2str(m)]);
end
%% 3.1 Save peak-to-trough and repolarization time NEW
ip_infor = {};
for cell = 1:n
    [rt,ptt,ip_1,ip_1y] = plot_ptt_new (cell,ts,FIGURE);
    FIGURE{cell,13} = ptt;
    FIGURE{cell,7} = rt;
    ip_infor{cell,1}=ip_1;
    ip_infor{cell,2}=ip_1y;
    close all
end
%% 3.2 Plotting peak-to-trough and repolarization
range = 21; % 8 cells per range
plotting =2;
for n = 1: range
    figure();
    for cell = (n-1)*8+1:n*8
        plot_ptt_2 (range,ts,FIGURE,plotting,ip_infor,cell,n)
    end
    if plotting ==2
        suptitle('Peak-to-trough Duration and Repolarization Time')
        xlabel('Time step (1/3k s)','position',[-50,-800,4],'fontsize',13);
        hAxis = axes('visible','off');
        h = text(-0.1,0.5,'Amplitude (v)');
        set(h,'fontsize',13,'rotation',90,'HorizontalAlignment','center');
    end
end
%% 4.1 Bursts, Pauses, and tonics analysis (whole-time scale)
% 472 not common
for n = 1
    spiketime = FIGURE{n,11}';
    ISI = diff(FIGURE{n,11});
    Spikes = spiketime/1000;
    N_min = 2;
    Steps = -3:0.005:1.5;
    p = 0.05;
    alpha = 0.05;
    plotting =0;
    [Bursts, Pauses] = RGSDetect(Spikes, N_min, Steps, p, alpha, plotting)
    FIGURE{n,14} = Bursts;
    FIGURE{n,15} = Pauses;
    % hist(ISI,[min(ISI):40:max(ISI)]);
    % xlabel ('ISI');
    % ylabel ('Number of spikes');
end
times = spiketime'
times = spiketime(400:600)
times = times'
rasterplot(times,1,5000)
% pd = fitdist(window,'Exponential'); % window is the bursts intra ISIs
% mean(pd);
% median (pd);
% var (pd);
% distributionFitter(window)
for n = 1:total_cell
    if isempty(FIGURE{n,14})
        FIGURE{n,14} = nan;
    else
        FIGURE{n,14} = FIGURE{n,14};
    end
end
for n = 1:total_cell
    if isempty(FIGURE{n,15})
        FIGURE{n,15} = nan;
    else
        FIGURE{n,15} = FIGURE{n,15};
    end
end
for n = 1:total_cell
    if isempty(FIGURE{n,9})
        FIGURE{n,9} = nan;
    else
        FIGURE{n,9} = FIGURE{n,9};
    end
end
%% 4.2 new Bursts and pauses detection (first 5 min)
for n = 472
    n
    last = max(find(FIGURE{n,11}<   300000)); % Change time-period!!!
    spiketime = FIGURE{n,11}(1:last,:);
    spiketime = spiketime';
    %ISI = diff(FIGURE{n,11});
    Spikes = spiketime/1000;
    N_min = 2;
    Steps = -3:0.005:1.5;
    p = 0.05;
    alpha = 0.05;
    plotting =0;
    try
        [Bursts, Pauses] = RGSDetect(Spikes, N_min, Steps, p, alpha, plotting)
    catch
        disp(['Inconsistent data in iteration, skipped n = ', num2str(n)]) % Detect unit that cannot be analyzed by RGSdetect
    end
    FIGURE{n,20} = Bursts;
    FIGURE{n,21} = Pauses;
end
% save ('D:\graduation thesis\results\FIGURE.mat','FIGURE'); % 340,472 problem
%% 4.3 Rastor plotting of representative neuron
n = 104
last = max(find(FIGURE{n,11}<300000));
spiketime = FIGURE{n,11}(1:last,:);
spiketime = spiketime';
Spikes = spiketime;
rasterplot(Spikes,1,8000)
set (gcf,'Position',[1500,700,1500,300],'PaperSize',[10 2])
print(gcf,'-dtiff','-r300',['D:\graduation thesis\results\' 'Spikes']);
%
% hold on
% x = 1000:2000;
% y = 1.15 *ones(1,length(x));
% plot(x,y,'k')
Spikes_2 = Spikes(92:end);
rasterplot(Spikes_2,10,8000)
BWindows = FIGURE{n,20}.Windows;
%PWindows = FIGURE{n,21}.Windows;
BWindows  = BWindows .* 1000;
%PWindows = PWindows .*1000;
for m = 1:length(BWindows)
    hold on
    x = BWindows(m,1):BWindows(m,2);
    y = 1.15 *ones(1,length(x));
    plot(x,y,'g')
end
%% 4.4 Fill the empty part with Nan in Bursts and Pauses for future analysis
for n = 1:total_cell
    if isempty(FIGURE{n,20})
        FIGURE{n,20} = nan;
    else
    end
end
for n = 1:total_cell
    if isempty(FIGURE{n,21})
        FIGURE{n,21} = nan;
    else
    end
end
%% 4.5 Extract features of the bursts, pauses and tonics
% Get more information: totoal duration, ISIs, ISIs_ins, Win, IPF,
% NumSpikes of Bursts, Pauses and Tonics
% spike within spikes
for n = 1:total_cell
    Bursts = FIGURE{n,20};
    Pauses = FIGURE{n,21};
    last = max(find(FIGURE{n,11}<300000));
    spiketime = FIGURE{n,11}(1:last,:);
    spiketime = spiketime';
    Spikes = spiketime/1000;
    Burstings = [];
    Pausings = [];
    % Bursts.IBF Bursts.NumSpikes
    if isstruct(Bursts)
        Bursts.BISIs = [];
        Bursts.BISIs_ins = [];
        Bursts.td =[];
        Bursts.win = [];
        % get corrected burstings
        for m = 1:length (Bursts.BurstingSpikes)
            Bursting = Bursts.BurstingSpikes{m};
            Burstings = [Burstings Bursting(1)];
        end
        a = [];
        for m= 1: length (Burstings)
            if find (Bursts.Windows(:,1) == Burstings(m))
                a = [a m];
            else
            end
        end
        % Get all spikes within Bursts
        Burstings = [];
        for x = 1: length(Bursts.BurstingSpikes)
            if find (x == a)
                Bursting = Bursts.BurstingSpikes{x};
                Burstings = [Burstings Bursting];
            else
            end
        end
        % Save withinBursts ISIs
        for x = 1: length(Bursts.BurstingSpikes)
            if find (x == a)
                Bursting = Bursts.BurstingSpikes{x};
                BISI = diff(Bursting);
                BISI_ins = log10(1./diff(Bursting));
                Bursts.BISIs = [Bursts.BISIs BISI];
                Bursts.BISIs_ins =[Bursts.BISIs_ins BISI_ins];
            end
        end
        Bursts.td = sum (Bursts.Windows(:,2) - Bursts.Windows(:,1)); % Burst duration
        Bursts.win = Bursts.Windows(:,2) - Bursts.Windows(:,1); % duration
    else
        Bursts = {};
        Bursts.BISIs = [];
        Bursts.BISIs_ins = [];
        Bursts.td =[];
        Bursts.win = [];
        Bursts.IBF = [];
        Bursts.NumSpikes =[];
    end
    % Pauses.IPF Pauses.NumSpikes
    if isstruct(Pauses)
        Pauses.PISIs = [];
        Pauses.PISIs_ins = [];
        Pauses.td = [];
        Pauses.win = [];
        % Get all spikes within Pauses
        for m = 1:length (Pauses.PausingSpikes)
            Pausing = Pauses.PausingSpikes {m};
            Pausings = [Pausings Pausing];
        end
        % save withinpause ISIs
        for m = 1: length(Pauses.PausingSpikes)
            Pausing = Pauses.PausingSpikes{m};
            PISI = diff(Pausing); % s
            PISI_ins = log10(1./diff(Pausing));
            Pauses.PISIs = [Pauses.PISIs PISI];
            Pauses.PISIs_ins =[Pauses.PISIs_ins PISI_ins];
        end
        Pauses.td = sum (Pauses.AllLengths);
        Pauses.win = Pauses.Windows(:,2) - Pauses.Windows(:,1);
    else
        Pauses = {};
        Pauses.td = [];
        Pauses.PISIs = [];
        Pauses.PISIs_ins = [];
        Pauses.win = [];
        Pauses.IPF = [];
        Pauses.NumSpikes = [];
    end
    Tonics.allspikes = setdiff(Spikes,[Burstings Pausings]);
    TISIs = setdiff(diff(Spikes),[Pauses.PISIs Bursts.BISIs]);
    if ~isempty(Bursts.BISIs)
        Tonics.TISIs = TISIs(TISIs> max(Bursts.BISIs));
    else
    end
    if ~isempty(Pauses.td) && ~isempty(Bursts.td)
        Tonics.td = 300 -Pauses.td - Bursts.td;
    end
    if  isempty(Pauses.td) && ~isempty(Bursts.td)
        Tonics.td = 300 - Bursts.td;
    end
    if isempty(Pauses.td) && isempty(Bursts.td)
        Tonics.td = 300;
    end
    FIGURE{n,22} = Bursts;
    FIGURE{n,23} = Pauses;
    FIGURE{n,24} = Tonics;
end
%% 4.6 Get and save the distribution information of features in Bursts, Pauses and Tonics
% mean,median,std,skewness,and kurtosis
for n = 1:total_cell
    Bursts = FIGURE{n,22};
    % form matrix
    if length (Bursts.BISIs') >1
        Bursts_ISImatrix = [Bursts.BISIs' Bursts.BISIs_ins'];
        Bursts_ISIfeature = [mean(Bursts_ISImatrix);median(Bursts_ISImatrix);std(Bursts_ISImatrix) ;skewness(Bursts_ISImatrix);kurtosis(Bursts_ISImatrix)];
    else
        Bursts_ISIfeature = [mean(Bursts.BISIs') mean(Bursts.BISIs_ins');median(Bursts.BISIs')  median(Bursts.BISIs_ins');
            std(Bursts.BISIs') std(Bursts.BISIs_ins');
            skewness(Bursts.BISIs') skewness(Bursts.BISIs_ins');kurtosis(Bursts.BISIs') kurtosis(Bursts.BISIs_ins')];
    end
    
    if length (Bursts.win)>1
        Bursts_spikematrix = [Bursts.win Bursts.IBF Bursts.NumSpikes];
        % Bursts_feature matrix
        Bursts_spikefeature = [mean(Bursts_spikematrix);median(Bursts_spikematrix);std(Bursts_spikematrix);
            skewness(Bursts_spikematrix);kurtosis(Bursts_spikematrix)];
    else
        Bursts_spikefeature = [mean(Bursts.win) mean(Bursts.IBF) mean(Bursts.NumSpikes);
            median(Bursts.win) median(Bursts.IBF) median(Bursts.NumSpikes);
            std(Bursts.win) std(Bursts.IBF) std(Bursts.NumSpikes);
            skewness(Bursts.win)  skewness(Bursts.IBF) skewness(Bursts.NumSpikes);
            kurtosis(Bursts.win) kurtosis(Bursts.IBF) kurtosis(Bursts.NumSpikes)];
    end
    % The average percentage of time (ms) spent bursting (s) (%)
    if ~isempty(Bursts.td)
        per_Bursts_windows = (Bursts.td/300)*100;
    else
        per_Bursts_windows = nan;
    end
    Bursts_ISIfeature = Bursts_ISIfeature(:);
    Bursts_spikefeature = Bursts_spikefeature(:);
    FIGURE{n,25} = [per_Bursts_windows Bursts_ISIfeature' Bursts_spikefeature'];
    Pauses = FIGURE{n,23};
    % form matrix
    if length (Pauses.PISIs') >1
        Pauses_ISImatrix = [Pauses.PISIs' Pauses.PISIs_ins'];
        Pauses_ISIfeature = [mean(Pauses_ISImatrix);median(Pauses_ISImatrix) ;
            std(Pauses_ISImatrix); skewness(Pauses_ISImatrix);kurtosis(Pauses_ISImatrix)];
    else
        Pauses_ISIfeature = [mean(Pauses.PISIs') mean(Pauses.PISIs_ins');median(Pauses.PISIs') median(Pauses.PISIs_ins');
            std(Pauses.PISIs') std(Pauses.PISIs_ins');
            skewness(Pauses.PISIs') skewness(Pauses.PISIs_ins');kurtosis(Pauses.PISIs') kurtosis(Pauses.PISIs_ins')];
    end
    
    if length (Pauses.win)>1
        Pauses_spikematrix = [Pauses.win Pauses.IPF Pauses.NumSpikes];
        % Pauses_feature matrix
        Pauses_spikefeature = [mean(Pauses_spikematrix);median(Pauses_spikematrix);std(Pauses_spikematrix);
            skewness(Pauses_spikematrix);kurtosis(Pauses_spikematrix)];
    else
        Pauses_spikefeature = [mean(Pauses.win) mean(Pauses.IPF) mean(Pauses.NumSpikes);
            median(Pauses.win) median(Pauses.IPF) median(Pauses.NumSpikes);
            std(Pauses.win) std(Pauses.IPF) std(Pauses.NumSpikes);
            skewness(Pauses.win)  skewness(Pauses.IPF) skewness(Pauses.NumSpikes);
            kurtosis(Pauses.win) kurtosis(Pauses.IPF) kurtosis(Pauses.NumSpikes)];
    end
    % The average percentage of time (ms) spent bursting (s) (%)
    if ~isempty(Pauses.td)
        per_Pauses_windows = (Pauses.td/300)*100;
    else
        per_Pauses_windows = nan;
    end
    Pauses_ISIfeature = Pauses_ISIfeature(:);
    Pauses_spikefeature = Pauses_spikefeature(:);
    FIGURE{n,26} = [per_Pauses_windows Pauses_ISIfeature' Pauses_spikefeature'];
    Tonics = FIGURE{n,24};
    per_Tonics_windows = (Tonics.td/300)*100;
    FIGURE{n,27} = [per_Tonics_windows];
end
%% 5.1 Import all features into a new table
table_raw = {}; % Create a cell to save this new feature table
table_new = []; % row: different units; column: different features
for n= 1:total_cell
    ptt =  FIGURE{n,13};
    rt = FIGURE{n,7};
    td1 = FIGURE{n,8};
    td2= FIGURE {n,9};
    hw = FIGURE {n,10};
    Burst_features = FIGURE{n,25}; %
    Pause_features = FIGURE{n,26}; %
    Tonic_features = FIGURE{n,27}; %
    mark = FIGURE{n,5};
    cell = [ptt rt td1 td2 hw Burst_features Pause_features Tonic_features mark];
    table_new = [table_new; cell];
end
table_raw {2,1} = table_new; % save the new feature table (58 features) in a cell called table {2,1}
%save ('D:\graduation thesis\results\table_raw.mat','table_raw');
%% 5.2 Save all features of bursts and pauses
for n = 1:total_cell
    if isstruct(FIGURE{n,14})
        Bursts = FIGURE{n,14};
        % burst length (s) [1x(number of burst)]
        Bursts_windows = Bursts.Windows(:,2) - Bursts.Windows(:,1);
        % coloumn intraburst frequency (HZ) [1x(number of burst)]
        Bursts_IBF = Bursts.IBF;
        % spikes per burst [1x(number of burst)]
        Bursts_NumSpikes = Bursts.NumSpikes;
        % form matrix
        Bursts_matrix = [Bursts_windows Bursts_IBF Bursts_NumSpikes];
        % 16 features now
        % Bursts_feature matrix 5X3 = 15
        Bursts_feature = [mean(Bursts_matrix);median(Bursts_matrix);std(Bursts_matrix);
            skewness(Bursts_matrix);kurtosis(Bursts_matrix)];
        % The average percentage of time (ms) spent bursting (s) (%)  1
        per_Bursts_windows = (sum (Bursts_windows)*1000/max(FIGURE {n,11}))*100;
    else
        Bursts_feature = nan;
        per_Bursts_windows =nan;
    end
    FIGURE{n,16} = Bursts_feature;
    FIGURE{n,17} = per_Bursts_windows;
end
for n = 1:total_cell
    if isstruct(FIGURE{n,15})
        Pauses = FIGURE{n,15};
        % burst length (s) [1x(number of burst)]
        Pauses_windows = Pauses.Windows(:,2) - Pauses.Windows(:,1);
        % coloumn intraburst frequency (HZ) [1x(number of burst)]
        Pauses_IPF = Pauses.IPF;
        % spikes per burst [1x(number of burst)]
        Pauses_NumSpikes = Pauses.NumSpikes;
        % form matrix
        Pauses_matrix = [Pauses_windows Pauses_IPF Pauses_NumSpikes];
        % 16 features now
        % Pauses_feature matrix 5X3 = 15
        Pauses_feature = [mean(Pauses_matrix);median(Pauses_matrix);std(Pauses_matrix);
            skewness(Pauses_matrix);kurtosis(Pauses_matrix)];
        % The average percentage of time (ms) spent bursting (s) (%)  1
        per_Pauses_windows = (sum (Pauses_windows)*1000/max(FIGURE {n,11}))*100;
    else
        Pauses_feature = nan;
        per_Pauses_windows =nan;
    end
    FIGURE{n,18} = Pauses_feature;
    FIGURE{n,19} = per_Pauses_windows;
end
%% 5.3 Form a features table (38 features) from FIGURE
% table = [data_ttp' data_rt' td1s' td2s' hws' bws' pws' IBFs_1' IBFs_2'
%IBFs_3'  IBFs_4'  IBFs_5'  IBFs_6'  IBFs_7'  IBFs_8'  IBFs_9'  IBFs_10'
% IBFs_11'  IBFs_12'  IBFs_13'  IBFs_14'  IBFs_15' IPFs_1' IPFs_2'
%IPFs_3'  IPFs_4'  IPFs_5'  IPFs_6'  IPFs_7'  IPFs_8'  IPFs_9'  IPFs_10'
% IPFs_11'  IPFs_12'  IPFs_13'  IPFs_14'  IPFs_15'  marks'];
table = feature_table (FIGURE); % Extract features from FIGURE to form a table
table_raw {1,1} = table; % 38 features table
%save ([save_path 'table_raw.mat'],'table_raw');
% T = array2table(table,'VariableNames',{'Trough-to-peak duration','Repolarization time','mark'});
