function [rt,ptt,ip_1,ip_1y] = plot_ptt_new (cell,ts,FIGURE) % Plot and save trough-to-peak in dataset
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
db1 = fnder (cs,1); % the first derivative
db2 = fnder(cs,2); % the second
y_prime1=ppval(db1,ts);
y_prime2= ppval(db2,y_prime1);
x_prime2 = [1:32];
figure();
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
ip_1y=nan;
else 
ip_1 = infp_2(1); % choose the first
ip_1y = ppval(cs,ip_1);
%y_gminx_lomax = y(find(y==ptt_trough_y)+1);
rt =  ip_1 - ptt_trough_x;
end         
end

