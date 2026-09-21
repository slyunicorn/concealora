%% ========================================================================
%  COVERT (PARAMETER-HOPPING CSS) LoRa  vs  STANDARD LoRa
%  ------------------------------------------------------------------------
%  Changes vs original:
%   1. Single source of truth for noise_power / tx power (CONFIG block only).
%   2. Fs/BW integer-oversampling sanity check.
%   3. Narrow-band ENERGY DETECTOR (radiometer) instead of "peak FFT" as the
%      detectability metric  ->  ROC (Pd vs Pfa), Pd-vs-eavesdropper-distance,
%      and TX-power-vs-Pd, all at a fixed false-alarm rate.
%   4. Decode accuracy broken out BY SPREADING FACTOR.
%   5. Error bars (binomial std) on the SNR-vs-success Monte-Carlo curve.
%   6. Trimmed redundant figures: time-domain merged into one figure,
%      static-spectrum figure removed (spectrogram covers it), signal/noise
%      bar chart replaced by a printed line.
%
%  NOTE on covertness: both waveforms are constant-envelope, so a WIDEBAND
%  radiometer sees equal energy for both. The detector below is CHANNELIZED
%  to a standard 125 kHz LoRa band, which is where the hopping/spreading
%  actually buys low probability of detection.
%
%  Requires: Signal Processing Toolbox (spectrogram only). Detector uses fft.
%% ========================================================================
clc; clear; close all;

%% ===================== CONFIG (single source of truth) =====================
Fs            = 1e6;                 % sample rate [Hz]
message       = 'HELLO JAY KRISHNA';
tx_pos        = [175, 181];
rx_pos        = [180, 180];          % intended receiver
fc            = 868e6;               % carrier [Hz]
c             = 3e8;

noise_power   = 1e-6;                % AWGN power  (defined ONCE, used everywhere)
tx_power_scale= 0.01;                % transmit power scale (alpha)
tx_amp        = sqrt(tx_power_scale);

ref_BW        = 125e3;               % eavesdropper's channelized detector bandwidth

pathloss      = @(d) (4*pi*d*fc/c).^2;   % free-space path loss (power)

SF_list       = [7 8 9 10 11 12];
BW_list       = [125e3 250e3 500e3];

% --- (2) integer-oversampling sanity check --------------------------------
osr = Fs ./ BW_list;
assert(all(osr == round(osr)), ...
    'Fs/BW must be an integer for every bandwidth (got %s).', mat2str(osr));

%% ===================== PARAMETER HOPPING SEQUENCES =====================
ascii_vals  = double(message);
num_symbols = length(ascii_vals);

rng(42);
SF_seq          = SF_list(randi(length(SF_list),1,num_symbols));
BW_seq          = BW_list(randi(length(BW_list),1,num_symbols));
freq_offset_seq = (rand(1,num_symbols)-0.5)*1e4;   % +/-5 kHz covert offset

symbols = zeros(1,num_symbols);
for k = 1:num_symbols
    symbols(k) = mod(ascii_vals(k), 2^SF_seq(k));
end

%% ===================== COVERT TX =====================
tx_signal      = [];
symbol_lengths = zeros(1,num_symbols);

for k = 1:num_symbols
    SF = SF_seq(k);  BW = BW_seq(k);  freq_offset = freq_offset_seq(k);

    symbol_duration = (2^SF)/BW;
    t  = 0:1/Fs:symbol_duration-1/Fs;
    Ns = length(t);
    symbol_lengths(k) = Ns;

    f0      = -BW/2;
    k_chirp = BW/symbol_duration;

    chirp_sig = exp(1j*2*pi*(f0*t + (k_chirp/2)*t.^2));           % base up-chirp
    chirp_sig = chirp_sig .* exp(1j*2*pi*(symbols(k)*BW/2^SF)*t); % symbol shift
    chirp_sig = chirp_sig .* exp(1j*2*pi*freq_offset*t);         % covert offset

    tx_signal = [tx_signal chirp_sig]; %#ok<AGROW>
end

%% ===================== CHANNEL (intended link) =====================
distance  = norm(tx_pos - rx_pos);
path_loss = pathloss(distance);

rx_signal = tx_signal * tx_amp / sqrt(path_loss);
noise     = sqrt(noise_power/2)*(randn(size(rx_signal)) + 1j*randn(size(rx_signal)));
rx_signal = rx_signal + noise;

signal_power = mean(abs(rx_signal).^2);
measured_SNR = 10*log10(signal_power / noise_power);

fprintf('Distance (TX->RX): %.2f m\n', distance);
fprintf('Signal power = %.3e  |  Noise power = %.3e\n', signal_power, noise_power);
fprintf('Measured SNR: %.2f dB\n', measured_SNR);

%% ===================== COVERT RX (matched, knows the hop schedule) =====================
decoded_symbols = zeros(1,num_symbols);
ptr = 1;
for k = 1:num_symbols
    SF = SF_seq(k);  BW = BW_seq(k);  freq_offset = freq_offset_seq(k);
    Ns = symbol_lengths(k);
    t  = 0:1/Fs:(Ns-1)/Fs;

    segment = rx_signal(ptr:ptr+Ns-1);  ptr = ptr + Ns;

    f0      = -BW/2;
    k_chirp = BW/((2^SF)/BW);

    segment   = segment .* exp(-1j*2*pi*freq_offset*t);          % undo offset
    ref_chirp = exp(-1j*2*pi*(f0*t + (k_chirp/2)*t.^2));
    dechirped = segment .* ref_chirp;

    [~, idx] = max(abs(fft(dechirped)));
    decoded_symbols(k) = mod(idx-1, 2^SF);
end
fprintf('\nCovert TX -> Covert RX: %s\n', char(max(32,min(126,decoded_symbols))));

%% ===================== NORMAL (baseline) TX =====================
SF_norm = 12;  BW_norm = 125e3;
symbol_duration_norm = (2^SF_norm)/BW_norm;
t_norm  = 0:1/Fs:symbol_duration_norm-1/Fs;
Ns_norm = length(t_norm);

normal_tx = [];
for k = 1:num_symbols
    sym       = mod(double(message(k)), 2^SF_norm);
    f0        = -BW_norm/2;
    k_chirp   = BW_norm/symbol_duration_norm;
    chirp_sig = exp(1j*2*pi*(f0*t_norm + (k_chirp/2)*t_norm.^2));
    chirp_sig = chirp_sig .* exp(1j*2*pi*(sym*BW_norm/2^SF_norm)*t_norm);
    normal_tx = [normal_tx chirp_sig]; %#ok<AGROW>
end

normal_rx = normal_tx * tx_amp / sqrt(path_loss);
normal_rx = normal_rx + sqrt(noise_power/2)*(randn(size(normal_tx)) + 1j*randn(size(normal_tx)));

%% ===================== NORMAL RX =====================
decoded_normal = zeros(1,num_symbols);
for k = 1:num_symbols
    segment   = normal_rx((k-1)*Ns_norm+1 : k*Ns_norm);
    ref_chirp = exp(-1j*2*pi*(-BW_norm/2*t_norm + (BW_norm/symbol_duration_norm/2)*t_norm.^2));
    [~, idx]  = max(abs(fft(segment .* ref_chirp)));
    decoded_normal(k) = mod(idx-1, 2^SF_norm);
end
fprintf('Normal TX -> Normal RX: %s\n', char(max(32,min(126,decoded_normal))));

%% ===================== NORMAL RX ON COVERT (interception fails) =====================
decoded_fail = [];  ptr = 1;
for k = 1:num_symbols
    if ptr+Ns_norm-1 > length(rx_signal), break; end
    segment   = rx_signal(ptr:ptr+Ns_norm-1);  ptr = ptr + Ns_norm;
    ref_chirp = exp(-1j*2*pi*(-BW_norm/2*t_norm + (BW_norm/symbol_duration_norm/2)*t_norm.^2));
    [~, idx]  = max(abs(fft(segment .* ref_chirp)));
    decoded_fail = [decoded_fail mod(idx-1, 2^SF_norm)]; %#ok<AGROW>
end
fprintf('Covert TX -> Normal RX: %s   (garbage = not decodable)\n', ...
        char(max(32,min(126,decoded_fail))));

%% ===================== SPECTROGRAM COMPARISON (kept) =====================
figure;
clim = [-160 -100];
subplot(2,1,1);
spectrogram(normal_rx,256,200,256,Fs,'yaxis'); caxis(clim); colormap jet;
title('Normal LoRa Spectrogram (visible chirps)','FontWeight','bold');
subplot(2,1,2);
spectrogram(rx_signal,256,200,256,Fs,'yaxis'); caxis(clim); colormap jet;
title('Covert LoRa Spectrogram (noise-like)','FontWeight','bold');

%% ===================== TIME DOMAIN (merged into one figure) =====================
figure;
subplot(2,1,1);
plot(real(tx_signal(1:2000))); hold on; plot(real(rx_signal(1:2000)));
legend('TX signal','RX signal'); grid on;
title('Time domain: TX vs RX (raw amplitude)'); xlabel('Samples'); ylabel('Amplitude');

subplot(2,1,2);
plot(real(tx_signal(1:2000)/max(abs(tx_signal)))); hold on;
plot(real(rx_signal(1:2000)/max(abs(rx_signal))));
legend('TX (norm.)','RX (norm.)'); grid on;
title('Time domain: TX vs RX (normalized)'); xlabel('Samples'); ylabel('Amplitude');

%% ===================== PARAMETER HOPPING (kept) =====================
figure;
subplot(2,1,1); plot(SF_seq,'-o'); title('Spreading-factor hopping');
xlabel('Symbol index'); ylabel('SF'); grid on;
subplot(2,1,2); plot(BW_seq/1e3,'-o'); title('Bandwidth hopping');
xlabel('Symbol index'); ylabel('BW (kHz)'); grid on;

%% ===================== DISTANCE vs SNR (intended link budget) =====================
distances = 1:5:200;
snr_vals  = zeros(size(distances));
for i = 1:length(distances)
    rx_tmp     = tx_signal * tx_amp / sqrt(pathloss(distances(i)));
    snr_vals(i)= 10*log10(mean(abs(rx_tmp).^2) / noise_power);
end
figure; plot(distances, snr_vals,'LineWidth',2); grid on;
title('Intended link: distance vs SNR'); xlabel('Distance (m)'); ylabel('SNR (dB)');

%% ===================== SNR vs SUCCESS  (+ per-SF accuracy, + error bars) =====================
SNR_range  = -20:1:-4;
num_trials = 50;

success    = zeros(size(SNR_range));               % message-level success prob.
SF_present = unique(SF_seq);
corr_SF    = zeros(numel(SF_present), numel(SNR_range));  % per-SF correct symbols
tot_SF     = zeros(numel(SF_present), numel(SNR_range));  % per-SF total symbols

sig_pow_ref = mean(abs(rx_signal).^2);

for i = 1:length(SNR_range)
    correct_count = 0;
    noise_pow     = sig_pow_ref / (10^(SNR_range(i)/10));

    for trial = 1:num_trials
        test_signal = rx_signal + sqrt(noise_pow/2)*(randn(size(rx_signal)) + 1j*randn(size(rx_signal)));

        decoded_test = zeros(1,num_symbols);  ptr = 1;
        for k = 1:num_symbols
            SF = SF_seq(k);  BW = BW_seq(k);  freq_offset = freq_offset_seq(k);
            Ns = symbol_lengths(k);
            t  = 0:1/Fs:(Ns-1)/Fs;

            segment   = test_signal(ptr:ptr+Ns-1);  ptr = ptr + Ns;
            segment   = segment .* exp(-1j*2*pi*freq_offset*t);
            ref_chirp = exp(-1j*2*pi*(-BW/2*t + (BW/((2^SF)/BW)/2)*t.^2));
            [~, idx]  = max(abs(fft(segment .* ref_chirp)));
            sym_est   = mod(idx-1, 2^SF);
            decoded_test(k) = sym_est;

            % ---- (4) per-SF symbol accuracy bookkeeping ----
            sfi = find(SF_present == SF, 1);
            tot_SF(sfi,i)  = tot_SF(sfi,i) + 1;
            corr_SF(sfi,i) = corr_SF(sfi,i) + (sym_est == symbols(k));
        end

        if strcmp(char(max(32,min(126,decoded_test))), message)
            correct_count = correct_count + 1;
        end
    end
    success(i) = correct_count / num_trials;
end

% (5) message-level success with binomial standard-error bars
se = sqrt(success.*(1-success)/num_trials);
figure;
errorbar(SNR_range, success, se,'-o','LineWidth',2); grid on;
title('SNR vs decoding success (message-level, \pm1 s.e.)');
xlabel('SNR (dB)'); ylabel('P(message correct)'); ylim([-0.05 1.05]);

% (4) per-SF symbol accuracy vs SNR
figure; hold on;
for s = 1:numel(SF_present)
    plot(SNR_range, corr_SF(s,:)./tot_SF(s,:), '-o','LineWidth',1.5);
end
grid on; ylim([-0.05 1.05]);
title('Per-symbol decode accuracy by spreading factor');
xlabel('SNR (dB)'); ylabel('Symbol accuracy');
legend(arrayfun(@(x) sprintf('SF %d',x), SF_present,'uni',0),'Location','SE');

%% ===================== EAVESDROPPER: CHANNELIZED ENERGY DETECTOR =====================
% Radiometer tuned to a standard 125 kHz LoRa channel at baseband.
Nwin      = 8192;      % detector integration window (samples)
NtrialsD  = 1500;      % Monte-Carlo trials for ROC
Pfa_target= 0.01;      % operating false-alarm rate for Pd sweeps

% H0 (noise only) statistic + operating threshold
T0     = collect_stats(zeros(1,2*Nwin), noise_power, Fs, ref_BW, Nwin, NtrialsD);
T0s    = sort(T0);
thr_op = T0s(max(1, ceil((1-Pfa_target)*numel(T0s))));

% Pick a demo eavesdropper distance where NORMAL is moderately detectable,
% so the ROC lands in an informative regime (auto-calibrated).
d_scan = 2:2:60;  pdn = zeros(size(d_scan));
for i = 1:numel(d_scan)
    sn     = normal_tx * tx_amp / sqrt(pathloss(d_scan(i)));
    pdn(i) = mean(collect_stats(sn, noise_power, Fs, ref_BW, Nwin, 300) > thr_op);
end
[~, bi] = min(abs(pdn - 0.85));
d_eve   = d_scan(bi);
fprintf('\nEavesdropper demo distance (auto-picked): %.1f m\n', d_eve);

sig_cov_eve  = tx_signal * tx_amp / sqrt(pathloss(d_eve));
sig_norm_eve = normal_tx * tx_amp / sqrt(pathloss(d_eve));

T1_cov  = collect_stats(sig_cov_eve,  noise_power, Fs, ref_BW, Nwin, NtrialsD);
T1_norm = collect_stats(sig_norm_eve, noise_power, Fs, ref_BW, Nwin, NtrialsD);

% --- ROC: Pd vs Pfa ---
thr  = linspace(min([T0 T1_cov T1_norm]), max([T0 T1_cov T1_norm]), 400);
Pfa  = arrayfun(@(x) mean(T0     > x), thr);
Pd_c = arrayfun(@(x) mean(T1_cov > x), thr);
Pd_n = arrayfun(@(x) mean(T1_norm> x), thr);

figure; hold on;
plot(Pfa, Pd_n,'LineWidth',2);
plot(Pfa, Pd_c,'LineWidth',2);
plot([0 1],[0 1],'k--');
grid on; axis([0 1 0 1]);
legend('Normal LoRa','Covert LoRa','Chance','Location','SE');
xlabel('P_{fa}'); ylabel('P_d');
title(sprintf('Channelized energy-detector ROC (eavesdropper @ %.0f m)', d_eve));

%% ===================== Pd vs EAVESDROPPER DISTANCE (fixed Pfa) =====================
d_list   = 1:2:40;  Ntr = 300;
Pd_dist_c = zeros(size(d_list));  Pd_dist_n = zeros(size(d_list));
for i = 1:numel(d_list)
    sc = tx_signal * tx_amp / sqrt(pathloss(d_list(i)));
    sn = normal_tx * tx_amp / sqrt(pathloss(d_list(i)));
    Pd_dist_c(i) = mean(collect_stats(sc, noise_power, Fs, ref_BW, Nwin, Ntr) > thr_op);
    Pd_dist_n(i) = mean(collect_stats(sn, noise_power, Fs, ref_BW, Nwin, Ntr) > thr_op);
end
figure; hold on;
plot(d_list, Pd_dist_n,'-o','LineWidth',2);
plot(d_list, Pd_dist_c,'-o','LineWidth',2);
grid on; ylim([-0.05 1.05]);
legend('Normal LoRa','Covert LoRa','Location','NE');
xlabel('Eavesdropper distance (m)'); ylabel('P_d');
title(sprintf('Detection probability vs eavesdropper distance (P_{fa}=%.0f%%)', 100*Pfa_target));

%% ===================== TX POWER vs Pd (replaces peak-FFT metric) =====================
tx_scale_sweep = [1e-5 1e-4 1e-3 1e-2 1e-1 1];
Pd_pow_c = zeros(size(tx_scale_sweep));  Pd_pow_n = zeros(size(tx_scale_sweep));
for i = 1:numel(tx_scale_sweep)
    a  = sqrt(tx_scale_sweep(i));
    sc = tx_signal * a / sqrt(pathloss(d_eve));
    sn = normal_tx * a / sqrt(pathloss(d_eve));
    Pd_pow_c(i) = mean(collect_stats(sc, noise_power, Fs, ref_BW, Nwin, Ntr) > thr_op);
    Pd_pow_n(i) = mean(collect_stats(sn, noise_power, Fs, ref_BW, Nwin, Ntr) > thr_op);
end
figure; hold on;
semilogx(tx_scale_sweep, Pd_pow_n,'-o','LineWidth',2);
semilogx(tx_scale_sweep, Pd_pow_c,'-o','LineWidth',2);
set(gca,'XScale','log'); grid on; ylim([-0.05 1.05]);
legend('Normal LoRa','Covert LoRa','Location','SE');
xlabel('TX power scale (\alpha)'); ylabel('P_d');
title(sprintf('TX power vs detection probability (P_{fa}=%.0f%%, eve @ %.0f m)', ...
      100*Pfa_target, d_eve));

%% ======================= LOCAL FUNCTIONS =======================
function T = band_energy(x, Fs, ref_BW)
% Energy of x inside a +/- ref_BW/2 baseband channel (ideal brick-wall).
    N = numel(x);
    X = fftshift(fft(x));
    f = linspace(-Fs/2, Fs/2, N);
    mask = abs(f) <= ref_BW/2;
    T = sum(abs(X(mask)).^2) / N;
end

function T = collect_stats(sig, noise_power, Fs, ref_BW, Nwin, Ntrials)
% Monte-Carlo detector statistics: random window of sig + fresh AWGN, then
% channelized band energy. Pass a zero vector for the noise-only (H0) case.
    if numel(sig) < Nwin
        sig = [sig zeros(1, Nwin - numel(sig) + 1)];
    end
    L = numel(sig);
    T = zeros(1, Ntrials);
    for tr = 1:Ntrials
        s   = randi(max(1, L - Nwin + 1));
        seg = sig(s:s+Nwin-1);
        n   = sqrt(noise_power/2)*(randn(1,Nwin) + 1j*randn(1,Nwin));
        T(tr) = band_energy(seg + n, Fs, ref_BW);
    end
end