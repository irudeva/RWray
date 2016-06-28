%%%%%  put in 1 for pole; cutoff and sampr are in time domain units.
%%%%%  pole is the order of the butterworth filter

function Zl=zfltr(Zu,pole,cutoff,sampr)
% USE: Zl=zfltr(Zu,pole,cutoff,sampr)
% Run a time series through a lowpass butterworth filter.
% Zu can be an array of time series to be filtered.
%
% >> Zu     = unfiltered series: dimensiond MxN where N is the
%             number of series and M is the time dimension
% >> pole   = No. of filtere poles
% >> cutoff = cutoff period in units of basic period (days, months years etc.)
% >> sampr  = sampling rate in same units.
%
% Calculate coefficients:
[b1,a1]=butter(pole,2*sampr/cutoff);
%
[M,N] = size(Zu);             % determin how many time series are given
L = fix(cutoff/sampr);        % determine length of padding
for n=1:N
    pl=flipud(Zu(2:L+1,n));   % padding on left
    pr=flipud(Zu(M-L:M-1,n)); % padding on right
    d(1:L)=pl'; d(L+1:M+L)=Zu(:,n)'; d(M+L+1:M+2*L)=pr';
%    d=Zu(:,n)';
    mean_d=mean(d');
    d=d-mean_d*ones(1,length(d));
    y=filter(b1,a1,d);        % filter forward
    ry=fliplr(y);             % flip filtered series
    fry=filter(b1,a1,ry);     % filter backwords
    fd=fliplr(fry);           % flip again
    fd=fd+mean_d*ones(1,length(d));
    Zl(:,n)=fd(L+1:M+L)';
%    Zl(:,n)=fd';
end
