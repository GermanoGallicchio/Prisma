function [spectra, pr_cfg] = pr_STFFT(signal, pr_cfg)
% time-frequency decomposition via Short-Time Fast Fourier Transform (STFFT)
% 
% USAGE:
%   [spectra, pr_cfg] = pr_STFFT(signal, pr_cfg)
%
% INPUTS:
%   signal          : 1D array of signal values (e.g., time series)
%   pr_cfg          : configuration structure with following fields:
%
%   pr_cfg.srate                       : sampling rate (Hz)
%   pr_cfg.STFFT_windowLength_pnt      : window length in points (i.e., samples) (default = 250)
%   pr_cfg.STFFT_windowStep_pnt        : spacing between window centers in points (i.e., samples) (default = 50)
%   pr_cfg.STFFT_zeroPadding           : zero-padding multiplier (default = 1, no padding)
%   pr_cfg.physicalAxis_vec            : vector of physical quantity values (e.g., time, position) matching input signal length (default = [1:nSamples])
%   pr_cfg.physicalAxis_units          : [OPTIONAL] string describing units of physicalAxis_vec (default = 'a.u.')
%
%
%   pr_cfg.sanityCheck_text            : text output checks (boolean), default = true
%   pr_cfg.sanityCheck_fig             : visual checks with plots (boolean), default = true
%
%
% OUTPUTS:
%   spectra           : 3D array (frequencies x time windows x 2) where dimension 3 contains [real, imaginary] parts
%                       Note: you can use complex() to convert to complex 2D: (frequencies x time windows)
%   pr_cfg            : updated configuration structure same as before plus additional field(s)
%   pr_cfg.decomp_method               : 'STFFT'
%
% NOTES:
%   - Input signal is centered (DC removed) within each window
%   - Full Hann window is applied before FFT for spectral leakage reduction
%   - Normalization: amplitude is divided by window length (excluding zero padding)
%   - Output contains complex FFT coefficients; use abs() for amplitude or abs().^2 for power
%   - Keep in mind amplitude is spread between positive and negative frequencies

%% inpute validation and defaults

% convert signal to row vector (for consistent windowing logic)
signal = signal(:)';
nSamples = length(signal);

% configuration must be a structure
if ~isstruct(pr_cfg)
    error('pr_cfg must be a structure. To start, provide pr_cfg = struct()');
end

% sampling rate is required
if ~isfield(pr_cfg, 'srate')
    error('pr_cfg.srate (sampling rate in Hz) is required. Check your signal');
end
srate = pr_cfg.srate;

% window length
% longer = better frequency resolution
if ~isfield(pr_cfg, 'STFFT_windowLength_pnt')
    pr_cfg.STFFT_windowLength_pnt = 250;  % Window length in points (samples) 
    warning('pr_cfg.STFFT_windowLength_pnt was not entered by the user. I will use %d as default', pr_cfg.STFFT_windowLength_pnt)
end

% window step
if ~isfield(pr_cfg, 'STFFT_windowStep_pnt')
    pr_cfg.STFFT_windowStep_pnt = 50;  % Step size in points (samples)
    warning('pr_cfg.STFFT_windowStep_pnt was not entered by the user. I will use %d as default', pr_cfg.STFFT_windowStep_pnt)
end

% zero padding
if ~isfield(pr_cfg, 'STFFT_zeroPadding')
    pr_cfg.STFFT_zeroPadding = 1;  % Zero-padding factor (1 = no padding, 8x = higher frequency resolution)
    warning('pr_cfg.STFFT_zeroPadding was not entered by the user. I will use %d as default', pr_cfg.STFFT_zeroPadding)
end

% validate zero padding
if pr_cfg.STFFT_zeroPadding < 1
    error('pr_cfg.STFFT_zeroPadding must be >= 1. Received: %d', pr_cfg.STFFT_zeroPadding);
end


% physical axis is required (e.g., time, position)
% if not provided, create a sequence of points [1, 2, ..., nSamples]
if ~isfield(pr_cfg, 'physicalAxis_vec') || isempty(pr_cfg.physicalAxis_vec)
    pr_cfg.physicalAxis_vec = 1:nSamples;  % Default: sequence of point indices
    warning('pr_cfg.physicalAxis_vec was not entered by the user. I will use [1:nSamples] as default')
end

% align physical axis to the signal
physicalAxis_input = pr_cfg.physicalAxis_vec(:)'; % enforce row vector
if length(physicalAxis_input) ~= nSamples
    error('physicalAxis_vec length (%d) must match signal length (%d)', length(physicalAxis_input), nSamples);
end

% physical axis units
if ~isfield(pr_cfg, 'physicalAxis_units')
    pr_cfg.physicalAxis_units = 'a.u.';  % Default units: arbitrary units
    warning('pr_cfg.physicalAxis_units was not entered by the user. I will use %s as default', pr_cfg.physicalAxis_units)
end

% sanity checks
if ~isfield(pr_cfg, 'sanityCheck_text')
    pr_cfg.sanityCheck_text = true;
    warning('pr_cfg.sanityCheck_text was not entered by the user. I will use %s as default', pr_cfg.sanityCheck_text)
end
if ~isfield(pr_cfg, 'sanityCheck_fig')
    pr_cfg.sanityCheck_fig = true;
    warning('pr_cfg.sanityCheck_fig was not entered by the user. I will use %s as default', pr_cfg.sanityCheck_fig)
end

%% shortcuts

physicalAxis_vec = pr_cfg.physicalAxis_vec;

STFFT_windowLength_pnt = pr_cfg.STFFT_windowLength_pnt;
STFFT_windowStep_pnt   = pr_cfg.STFFT_windowStep_pnt;
STFFT_zeroPadding      = pr_cfg.STFFT_zeroPadding;

STFFT_range_pnt = [1, nSamples];

%% compute window parameters

% windows' centers, spaced every STFFT_windowStep_pnt samples
STFFT_windowCentre_pnt = STFFT_range_pnt(1) : STFFT_windowStep_pnt : STFFT_range_pnt(end);
STFFT_windowCentre_units = physicalAxis_vec(STFFT_windowCentre_pnt);
STFFT_noWindows = numel(STFFT_windowCentre_pnt);

% windows beginnings and ends
STFFT_windowStart_pnt = STFFT_windowCentre_pnt - STFFT_windowLength_pnt/2;
STFFT_windowEnd_pnt   = STFFT_windowCentre_pnt + STFFT_windowLength_pnt/2;
STFFT_windowEnd_pnt   = STFFT_windowEnd_pnt - 1; % to match STFFT_windowLength_pnt


% sanity check window length
winIdx = 1; % use this window as example
if STFFT_windowLength_pnt~=length(STFFT_windowStart_pnt(winIdx):STFFT_windowEnd_pnt(winIdx))
    error('likely a bug: the window length does not match what the user requested')
end

% combine windows' beginnings, centers, and ends (each row = a window)
STFFT_windowStartCentreEnd = [STFFT_windowStart_pnt' STFFT_windowCentre_pnt' STFFT_windowEnd_pnt'];

% compute overlap between contiguous windows
STFFT_overlapPnt = mean(STFFT_windowStartCentreEnd(1:end-1,3) - STFFT_windowStartCentreEnd(2:end,1));
STFFT_overlapPerc = STFFT_overlapPnt / STFFT_windowLength_pnt * 100;

%% compute frequency parameters
% frequency-bin spacing (precision) after zero-padding
% e.g., 0.5 s window -> 2 Hz bins; 8x zero-padding -> 0.25 Hz bins.

STFFT_HzPrecision = srate / (STFFT_windowLength_pnt * STFFT_zeroPadding);

% meaningful frequencies for real signals span DC to Nyquist.
% Zero-padding refines the grid but does not add new spectral content.
STFFT_nyquist = srate / 2;
STFFT_frex = linspace(0, STFFT_nyquist, floor(STFFT_windowLength_pnt*STFFT_zeroPadding/2)+1);

% sanity check: frequency precision
if abs(mean(diff(STFFT_frex)) - STFFT_HzPrecision) > eps*100
    warning('Frequency precision mismatch: computed %.4f Hz, expected %.4f Hz', ...
        mean(diff(STFFT_frex)), STFFT_HzPrecision);
end

%% sanity checks: on screen, print parameter summary

if pr_cfg.sanityCheck_text
    fprintf('\n=== STFFT: Parameters to be used ===\n');
    fprintf('Sampling rate: %.1f Hz\n', srate);
    fprintf('Window length: %d points (%.3f units at %.1f Hz)\n', STFFT_windowLength_pnt, STFFT_windowLength_pnt/srate, srate);
    fprintf('Window step: %d points (%.3f units at %.1f Hz)\n', STFFT_windowStep_pnt, STFFT_windowStep_pnt/srate, srate);
    fprintf('Number of windows: %d\n', STFFT_noWindows);
    fprintf('Window overlap: %.1f%%\n', STFFT_overlapPerc);
    fprintf('Zero padding multiplier: %d\n', STFFT_zeroPadding);
    fprintf('Frequency precision: %.4f Hz\n', STFFT_HzPrecision);
    fprintf('Frequency range: [%.1f %.1f] Hz\n', STFFT_frex(1), STFFT_frex(end));
    fprintf('=====================================\n\n');
end

%% Compute STFFT on each window

% initialize complex spectrum: DC..Nyquist for each window
STFFT_Data_complex = nan(length(STFFT_frex), STFFT_noWindows, 'like', 1i);

if pr_cfg.sanityCheck_fig
    figHandle = figure('Name', 'STFFT sanity checks', 'NumberTitle', 'off');
else
    figHandle = [];
end


for winIdx = 1:STFFT_noWindows
        %% identify window

        % shortcut: beginning and end of each window
        win_start = round(STFFT_windowStart_pnt(winIdx));
        win_end   = round(STFFT_windowEnd_pnt(winIdx));
        
        % don't compute FFT if not enough data
        % TO DO: in the future, I might reflect the signal but keep track
        % of boundaries similarly to COI for wavelet
        if win_start < 1 || win_end > nSamples
            continue;
        end

        % extract windowed segment
        STFFT_dataPre = signal(win_start:win_end);
        
        % mena center to remove DC to reduce leakage and center amplitude
        STFFT_dataPre = STFFT_dataPre - mean(STFFT_dataPre);
        STFFT_dataPre_raw = STFFT_dataPre; % store a copy before windowing

        %% apply windowing function

        % Hann window to reduce sidelobes; transpose to match row
        STFFT_dataPre = STFFT_dataPre .* hann(length(STFFT_dataPre))';
        
        %% visualize windowed data

        % TO DO: let the user choose whether they want to see all windows
        % or a random subset of them (e.g., a percentage)

        if ~isempty(figHandle)
            subplot(2, 3, 1:2)

            t = linspace(win_start,win_end,STFFT_windowLength_pnt);
            plot(t, STFFT_dataPre_raw, 'b-', 'DisplayName', 'Raw signal'); hold on;
            plot(t, STFFT_dataPre, 'r-', 'DisplayName', 'Hann windowed'); hold off
            set(gca,'XLim',t([1 end]))
            ylabel('Amplitude');
            xlabel('Time (s)');
            title(sprintf(['Window ' num2str(winIdx) ': pre and post Hann window']));
            legend('Location', 'northeast');
            grid on;
            
            subplot(2, 3, 3)
            t_win = (1:length(STFFT_dataPre_raw)) / srate;
            plot(t_win, hann(length(STFFT_dataPre_raw))', 'k-', 'LineWidth', 2);
            ylabel('Window Value');
            xlabel('Time (s)');
            title('Hann window function');
            xlim([t_win(1) t_win(end)]);
            grid on;

            
        end

        %% FFT on this window

        % FFT with zero-padding for finer frequency resolution
        STFFT_dataPost = fft(STFFT_dataPre, length(STFFT_dataPre) * STFFT_zeroPadding);

        % normalize by window length (not the padded length);
        % to keep amplitudes comparable regardless of padding factor
        STFFT_dataPost = STFFT_dataPost / length(STFFT_dataPre);

        if ~isempty(figHandle)
            subplot(2, 3, 4:5)
            plot(STFFT_frex, abs(STFFT_dataPost(1:length(STFFT_frex))), 'b-', 'LineWidth', 1.5);
            xlabel('Frequency (Hz)');
            ylabel('|FFT| Amplitude');
            title('FFT Amplitude Spectrum');
            xlim([0 min(50, STFFT_frex(end))]);
            grid on;
        end
        
        pause(0.05)

        %% store values

        STFFT_Data_complex(:, winIdx) = STFFT_dataPost(1:length(STFFT_frex));

end

%% convert complex output to 3D real/imaginary representation
% to reduce memory needed for this variable

spectra = nan(length(STFFT_frex), STFFT_noWindows, 2);
spectra(:, :, 1) = real(STFFT_Data_complex);  % real part
spectra(:, :, 2) = imag(STFFT_Data_complex);  % imaginary part

%% update pr_cfg with computed parameters (for reference and reproducibility)

pr_cfg.decomp_method = 'STFFT';
pr_cfg.STFFT_freq_Hz = STFFT_frex;              % Frequency axis for results
pr_cfg.physicalAxis_output = STFFT_windowCentre_units;
pr_cfg.STFFT_overlapPerc = STFFT_overlapPerc;
pr_cfg.STFFT_overlapPnt = STFFT_overlapPnt;

end