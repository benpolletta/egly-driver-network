function PowByThetaPlot(datafile, fs, freqs, no_cycles)
% datafile = name of data file
% fs = sampling freq. of data
% freqs = vector of frequencies for plot
% no_cycles = vector of cycle lengths for each frequency (same size as
% freqs)

data = load(datafile);
time = (1:length(data))/fs;

no_thetas = max(time)/.25;
length_theta = .25*fs;
theta_phase = 360*time(1:length_theta)/max(time(1:length_theta));

label = sprintf('%s_%.2gto%.2gHz_%.2gto%.2gcyc', datafile, min(freqs), max(freqs), min(no_cycles), max(no_cycles))

if exist([label, '.mat'], 'file')

    datamat = load([label, '.mat']);
    ws = datamat.ws;

else

    ws = wavelet_spectrogram(data, fs, freqs, no_cycles, 1, label);
    save(label, 'data', 'time', 'fs', 'freqs', 'no_cycles', 'ws', '-v7.3')

end

ws_folded = reshape(ws, [length_theta, no_thetas, length(freqs)]);
ws_tmean = squeeze(nanmean(abs(ws_folded), 2));

figure
imagesc(theta_phase, freqs, ws_tmean')
axis xy
set(gca, 'FontSize', 20)
saveas(gcf, [label, '.fig'])

ws_mean_mat = ones(size(ws_tmean))*diag(nanmean(ws_tmean));
ws_pct_tmean = (ws_tmean - ws_mean_mat)./ws_mean_mat;

figure
imagesc(theta_phase, freqs, ws_pct_tmean')
axis xy
set(gca, 'FontSize', 20)
saveas(gcf, [label, '_pct.fig'])