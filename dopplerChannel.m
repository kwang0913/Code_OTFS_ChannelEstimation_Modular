function y = dopplerChannel(x, fs, chanParams)
%DOPPLERCHANNEL Apply multipath channel with delays, Doppler shifts, and gains.
    numPaths = length(chanParams.pathDelays);  % delays are in samples (integer taps)
    maxPathDelay = max(chanParams.pathDelays);
    txOutSize = length(x);

    y = zeros(txOutSize + maxPathDelay, 1);

    for k = 1:numPaths
        pathOut = zeros(txOutSize + maxPathDelay, 1);
        pathOut(1+chanParams.pathDelays(k):chanParams.pathDelays(k)+txOutSize) = x;
        pathOut = frequencyOffset(pathOut, fs, chanParams.pathDopplerFreqs(k));  % MATLAB Comm Toolbox: applies exp(j*2*pi*f*t)
        pathOut = pathOut * chanParams.pathGains(k);
        y = y + pathOut;
    end
end
