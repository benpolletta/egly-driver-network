function SelectedPowByThetaPlots
% datafile = name of data file
% fs = sampling freq. of data
% freqs = vector of frequencies for plot
% no_cycles = vector of cycle lengths for each frequency (same size as
% freqs)

fs = 10^5;

prefixes = {'LIP', 'FEFvm'};

frequencies = {9:60, 9:40};

suffix = '_LFP_V_RS';

all_cycles = linspace(4,12,52);

for p = 1:2
    
    freqs = frequencies{p};
    cycles = all_cycles(1:min(length(all_cycles),length(freqs)));
    
    datafile = [prefixes{p}, suffix];
    
    label = sprintf('%s_%.2gto%.2gHz_%.2gto%.2gcyc', datafile, min(freqs), max(freqs), min(cycles), max(cycles))
    
    if exist([label, '.mat'], 'file')
        
        datamat = load([label, '.mat']);
        ws = datamat.ws;
        
    else
        
        data = load(datafile);
        time = (1:length(data))/fs;
        ws = wavelet_spectrogram(data, fs, freqs, cycles, 1, label);
        save(label, 'data', 'time', 'fs', 'freqs', 'no_cycles', 'ws', '-v7.3')
        
    end
    
    time = (1:size(ws, 1))/fs;
    
    no_thetas = max(time)/.25;
    length_theta = .25*fs;
    theta_phase = 360*time(1:length_theta)/max(time(1:length_theta));
    
    ws_folded = reshape(ws, [length_theta, no_thetas, length(freqs)]);
    ws_tmean = squeeze(nanmean(abs(ws_folded), 2));
    
    ws_mean_mat = ones(size(ws_tmean))*diag(nanmean(ws_tmean));
    ws_pct_tmean = (ws_tmean - ws_mean_mat)./ws_mean_mat;
    
    subplot(2, 1, p)
    imagesc(theta_phase, freqs, ws_pct_tmean')
    axis xy
    set(gca, 'FontSize', 20)
    title(prefixes{p})
    
end

saveas(gcf, 'selected_PAC_pct.fig')
saveas(gcf, 'selected_PAC_pct', 'eps')

end