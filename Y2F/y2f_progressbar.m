function status = y2f_progressbar(status,i,N,width)
% Displays progress bar in command window

bar = repmat(sprintf('='),1,round(i/N*width));
msg = sprintf('[%-30s] %d/%d', bar, i, N);
fprintf([status, msg]);
status = repmat(sprintf('%c',8), 1, length(msg));
