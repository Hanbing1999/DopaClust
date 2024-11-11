function [td_sp, td_sp_2, hw, cs, x_half,xStart_sp,xmax_sp,xStart_sp_2,xmax_sp_2] = plot_td_2(cell,ts,FIGURE)
data = FIGURE{cell,4}';
%% spline_cs 
% [ptt:peak-to-trough duration; ttp: Trough-to-peak duration;...
% td_sp: total duration; td_sp_2: total duration with 90% amplitude; amplitude; ]

bl = mean ([data(1:5) data(27:32)]); %  Calculate the baseline based on the first five and the last five data points 
cs = csapi (ts,data);
x = fnzeros (fnder(cs));% Get local extreme value [fnder:Differentiate function; fnzeros: roots of spline]
x = unique(x(:)); % omit repeated numbers
y = fnval(cs,x); % Evaluate spline function, obtain the value

% Calculate the roots for baselines by using SLM tool
fnmin (cs,ts); % find the minumum in a specific interval 
xStart_sp = slmsolve (cs,bl); % calculate the roots for baselines by using SLM tool
xStart_sp = xStart_sp(1); % Select the real start point

format long g 
fnval (cs,xStart_sp);
maksimum_ab = max(abs(fnval(cs,x))); % get the absolute number of local extreme

% Distinguish positive-going or negative-going
if isempty(find (y==maksimum_ab))  
    % if cannot find (logical is 1), it means that the waveform is negative-going 
    xmax_sp = x (find (y==-maksimum_ab));
    maksimum = -maksimum_ab; 
else 
maksimum = max(fnval(cs,x)); % get the maximum of local extreme
xmax_sp= x (find (y==maksimum));
end 

%FIGURE{cell,8} = [maksimum xmax_sp]; 
 
% To make sure the start point is before the max point
if xStart_sp < xmax_sp 
    xStart_sp = xStart_sp;
else 
    xStart_sp = slmsolve (cs,mean(data (1:5))); % only use the first five points as start points
    xStart_sp = xStart_sp(1);
end 

% Calculate two sets of total duration (1. from start to max. 2. 90% amplitude)
amplitude = maksimum - bl; 
pre_ap = bl + 0.1*amplitude;
tail_ap = maksimum - 0.1*amplitude;

xStart_sp_2 = slmsolve (cs,pre_ap,0);
xStart_sp_2 (xStart_sp_2 >xmax_sp) = []; % delete the solution of the second start point larger than the max
if length (xStart_sp_2) > 1 % if there are multiple solutions 
xStart_sp_2 = xStart_sp_2 (find (xStart_sp_2>xStart_sp));
xStart_sp_2 = xStart_sp_2(1);
else 
    xStart_sp_2 = xStart_sp_2;
end
xmax_sp_2 =  slmsolve (cs,tail_ap);
xmax_sp_2 (xmax_sp_2>xmax_sp) = []; % delete unreasonable solution

td_sp = xmax_sp - xStart_sp; % total duration
td_sp_2 = xmax_sp_2 - xStart_sp_2; % total duration_2

% Half-width
a_half = bl + amplitude*0.5; 
x_half = slmsolve (cs,a_half);
hw = x_half(2) - x_half (1);
if length (x_half) >3
   x_half = x_half (3:4);
end

end
