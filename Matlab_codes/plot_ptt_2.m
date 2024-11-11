function plot_ptt_2 (range,ts,FIGURE,plotting,ip_infor,cell,n) % Plot and save trough-to-peak in dataset

data = FIGURE{cell,4}';
bl = mean ([data(1:5) data(27:32)]); %  calculate the baseline based on the first five and the last five data points
cs = csapi (ts,data);
x = fnzeros (fnder(cs));% Get local extreme value [fnder:Differentiate function; fnzeros: roots of spline]
x = unique(x(:)); % omit repeated numbers
y = fnval(cs,x); % Evaluate spline function, obtain the value
% Select the real peak (not smaller than the ...)
maksimum_ab = max(abs(fnval(cs,x))); % get the absolute number of local extreme
% Distinguish positive-going or negative-going
if isempty(find (y==maksimum_ab))
    % if cannot find (logical is 1), it means that the waveform is negative-going
    gmin_y = min(y); % global minimum
    gmin_x = x (find(y==gmin_y)); % find the minumum in a specific interval
    
    xmax_sp = x (find (y==-maksimum_ab));
    ptt_peak_x =   xmax_sp;
    maksimum = -maksimum_ab;
    ptt_peak_y =  maksimum; % negative-going, the minimum
    
    if find(y==gmin_y)+1 > length (y) % In case there is no following local extreme
        ptt = nan;
    else
        ptt_trough_x = x(((find(y==gmin_y)+1)));% find the following local maximum for the minimum
        ptt = ptt_trough_x - ptt_peak_x;
        ptt_trough_y = y(((find(y==gmin_y)+1)));
    end
    
else
    % if the waveform is positive-going
    ptt_trough_y = min(y); % global minimum
    ptt_trough_x = x (find(y==ptt_trough_y)); % find the minumum in a specific interval
    if find(y==ptt_trough_y)-1 <1 % In case there is no following local extreme
        ptt = nan;
    else
        maksimum = max(fnval(cs,x)); % get the maximum of local extreme
        ptt_peak_y =  maksimum;
        ptt_peak_x= x (find (y==maksimum));
        %ptt_peak_x = x(((find(y==ptt_trough_y)-1))); % find the local maximum before the global minimum
        ptt = ptt_trough_x - ptt_peak_x;  % Peak-to-trough duration
    end
end
% repolarization time
if plotting == 1
    db1 = fnder (cs,1); % the first derivative
    db2 = fnder(cs,2); % the second
    y_prime1=ppval(db1,ts);
    y_prime2= ppval(db2,y_prime1);
    x_prime2 = [1:32];
    %subplot (4,4,2*cell -1);
    fnplt (db2)
    refline(0,0)
    title ('The Second Derivatives: Piecewise Linear')
    infp = slmsolve(db2); % find the crossing points (linear equations)
    % subplot (4,4,2*cell);
    % fnplt (cs);
    len=length (infp);
    infp_2=infp(find (infp > ptt_trough_x )); % find the inflection points after the local peak defined previously
    if isempty(infp_2)
        ip_1  = nan;
        %y_gminx_lomax = y(find(y==ptt_trough_y)+1);
        rt = nan;
    else
        ip_1 = infp_2(1); % choose the first
        ip_1y = ppval(cs,ip_1);
        %y_gminx_lomax = y(find(y==ptt_trough_y)+1);
        rt =  ip_1 - ptt_trough_x;
    end
end

if plotting == 2
    % PLOTTING: PEAK TO TROUGH
    if isempty(ptt)
        subplot (2,4,cell-8*(n-1))
        fnplt (cs)
        hold on
        %plot (ts, data,'x','color','[.667 .667 1]');
        plot (ptt_trough_x,ptt_trough_y,'b.','markersize',20)
        legend ('Data','Trough')
        title ('Trough-to-peak duration without peak')
    else
        subplot (2,4,cell-8*(n-1))
        fnplt (cs)
        ylim([-600 600])
        hold on
        % plot (ts, data,'x','color','[.667 .667 1]');
        plot(ptt_peak_x, ptt_peak_y,'r.','markersize',25);
        plot (ptt_trough_x,ptt_trough_y,'b.','markersize',25)
    end
    plot (ip_infor{cell,1},ip_infor{cell,2},'.','markersize',25,'color','#22BB3B')
    if cell == 1
        h = legend ('Data','Peak','Trough','Inflection point');
        set(h,'Box','off');
    else
    end
end
end

