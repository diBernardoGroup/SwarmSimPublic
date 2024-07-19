function [] = smooth_hist(data,nbins)


[heights,centers] = histcounts(data,nbins);                             % Get the counts and bins of the histogram


n = length(centers);                                                    % Number of bins
w = centers(2)-centers(1);                                              % Size of the bins
t = linspace(centers(1)-w/2,centers(end)+w/2,n+1);                      % create a vector spanning all the bins

dt = diff(t);                                                           % Get the size of the bin
Fvals = heights;%cumsum([0,heights.*dt]);                                        % Cumulative density function
% F = spline(t, [0, Fvals, 0]);                                           % Interpolate the CDF using splines
F = spaps(t, [0, Fvals, 0],1000);

% DF = fnder(F);                                                          % computes the CDF derivative (PDF)
fnplt(F, 'r', 2)                                                       % plots the smoothed histogram

end