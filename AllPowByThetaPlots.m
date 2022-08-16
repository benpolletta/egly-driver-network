function AllPowByThetaPlots
% datafile = name of data file
% fs = sampling freq. of data
% freqs = vector of frequencies for plot
% no_cycles = vector of cycle lengths for each frequency (same size as
% freqs)

fs = 10^5;

prefixes = {'LIP', 'FEFvm'};

frequencies = {9:60, 9:40};

suffixes = {'_RS_spikes','_LFP_I_RS', '_LFP_V_RS'};

cycle_cell={4, linspace(4,12,52), 6, linspace(6,18,52), 8, linspace(8,24,52)};

for p = 1:2

    for s = 1:3

        for c = 1:6

            freqs = frequencies{p};
            no_cycles = cycle_cell{c}(1:min(length(cycle_cell{c}),length(freqs)));

            datafile = [prefixes{p}, suffixes{s}];

            label = sprintf('%s_%.2gto%.2gHz_%.2gto%.2gcyc', datafile, min(freqs), max(freqs), min(no_cycles), max(no_cycles))

            if exist([label, '.mat'], 'file')

                datamat = load([label, '.mat']);
                ws = datamat.ws;

            else

                data = load(datafile);
                time = (1:length(data))/fs;
                ws = wavelet_spectrogram(data, fs, freqs, no_cycles, 1, label);
                save(label, 'data', 'time', 'fs', 'freqs', 'no_cycles', 'ws', '-v7.3')

            end

            time = (1:size(ws, 1))/fs;

            no_thetas = max(time)/.25;
            length_theta = .25*fs;

            ws_folded = reshape(ws, [length_theta, no_thetas, length(freqs)]);
            ws_tmean = squeeze(nanmean(abs(ws_folded), 2));

            figure(2*p-1)
            subplot(3, 6, (s-1)*6+c)
            imagesc(time(1:length_theta), freqs, ws_tmean')
            axis xy
            set(gca, 'FontSize', 20)

            ws_mean_mat = ones(size(ws_tmean))*diag(nanmean(ws_tmean));
            ws_pct_tmean = (ws_tmean - ws_mean_mat)./ws_mean_mat;

            figure(2*p)
            subplot(3, 6, (s-1)*6+c)
            imagesc(time(1:length_theta), freqs, ws_pct_tmean')
            axis xy
            set(gca, 'FontSize', 20)

        end

    end
    
    figure(2*p-1)
    saveas(gcf, [prefixes{p}, '.fig'])
    figure(2*p)
    saveas(gcf, [prefixes{p}, '_pct.fig'])

end