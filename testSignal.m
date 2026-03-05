%% validation of time-frequency analysis with Prisma
% - create a test signal with known properties
% - send it to time-frequency decomposition via Prisma, 
% - inspect the results to check they match what expected

%% add toolbox folder to path: Prisma

% relative path to toolbox
targetFolder = fullfile('work','functions','Prisma');
found = false;

if ispc % --- Windows, loop over drive letters A: to Z: ---
    letterVector = 65:90;   % ASCII codes for A–Z
    for diskLbl = string(char(letterVector) + ":")'
        addOnPath = fullfile(diskLbl, targetFolder);
        if exist(addOnPath,'dir')
            addpath(addOnPath)
            disp(['added to path: ' addOnPath])
            found = true;
            break
        end
    end

elseif isunix % --- Linux, look under /media/<user>
    mountFolder = fullfile('/media', getenv('USER'));
    % check one level down for drives, then into work/functions/PhysioExplorer
    matches = dir(fullfile(mountFolder, '*', targetFolder));
    if ~isempty(matches)
        addOnPath = matches(1).folder;
        addpath(addOnPath)
        disp(['added to path: ' addOnPath])
        found = true;
    end
        
end

if ~found
    warning([ targetFolder ' not added to path: folder not found'])
end

%% select signal type

signalType = 1;  % simple chirp

%% time and frequency setup
% common to all signal types

% sampling rate
fs = 2000; % Hz
dt = 1/fs; % delta t

% time info
timeRange = [0 5]; % seconds
timeVec = timeRange(1):dt:timeRange(end);
timeVec(end) = [];
N = length(timeVec);

% frequency info
df = fs/(N-1);
f_Fourier = (0:N-1)*df;
f_Fourier_shifted = (-(N/2) : (N/2-1)) * df;

%% generate test signal

switch signalType
    
    case 1
        %% signal 1: simple chirp with variable amplitude
        fprintf('signal type 1: simple logarithmic chirp\n');
        freqRange = [1 15];
        amplRange = [3 1];
        a = linspace(amplRange(1), amplRange(end), N);
        X1 = a' .* chirp((0:N-1)*dt, freqRange(1), (N-1)*dt, freqRange(end), 'logarithmic')';
        groundTruth_Lbl = 'Logarithmic chirp: 1-15 Hz, decreasing amplitude';
end

% frequency domain FFT (for all signal types)
X_fft = fft(X1);
X_abs = abs(X_fft)/N;
X_pow = X_abs.^2;

%% figure: see the test signal in time and frequency (FFT) domains

figure(1); clf
set(gcf, 'Position', [100 100 1200 500]);

% time domain
subplot(2,1,1)
plot(timeVec, X1, 'LineWidth', 1.5);
%xlabel(['time [' pr_cfg.physicalAxis_units ']'], 'FontSize', 12);
xlabel(['time [s]'], 'FontSize', 12);
ylabel('x_t [a.u.]', 'FontSize', 12);
title('time domain');
grid on;
xlim(timeRange);

% frequency domain (Fourier)
subplot(2,1,2)
plot(f_Fourier, X_abs, 'LineWidth', 1.5);
xlabel('frequency [Hz]', 'FontSize', 12);
ylabel('|X_f| [a.u.]', 'FontSize', 12);
title('frequency domain (FFT)');
grid on;
xlim([0 50]); 

% Add ground truth as suptitle
sgtitle(sprintf('signal type %d: %s', signalType, groundTruth_Lbl), 'FontSize', 12, 'FontWeight', 'normal');


%% 1- STFFT

% Configure STFFT parameters
pr_cfg = struct();
pr_cfg.srate = fs;                       % Sampling rate
pr_cfg.physicalAxis_vec = timeVec;       % Time axis for output
pr_cfg.physicalAxis_units = 's';         % Time units
pr_cfg.STFFT_windowLength_pnt = 2000;    % Window length (e.g., 500 points = 250 ms at 2000 Hz)
pr_cfg.STFFT_windowStep_pnt = 100;       % Window step (100 points = 50 ms at 2000 Hz)
pr_cfg.STFFT_zeroPadding = 8;            % 8x zero padding for frequency resolution


pr_cfg.sanityCheck_fig = false;

% run STFFT
[spectra, pr_cfg] = pr_STFFT(X1, pr_cfg);

% extract amplitude from complex spectra
complex_spectra = complex(spectra(:,:,1), spectra(:,:,2));
amplitude = abs(complex_spectra);
amplitude = amplitude * 2;
power = amplitude.^2;

% get time-frequency axes from output
freq_axis = pr_cfg.STFFT_freq_Hz;
time_axis_stfft = pr_cfg.physicalAxis_output;

%% Visualize STFFT results

% TO DO: convert this to its own function (e.g., pr_view)


figure(2); clf
set(gcf, 'Position', [150 50 1200 600]);

% Time-frequency spectrogram (raw amplitude in original units)
imagesc(time_axis_stfft, freq_axis, amplitude);
axis xy; 
colorbar;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Frequency (Hz)', 'FontSize', 12);
title('STFFT');
ylim([0 min(30, max(freq_axis))]);
grid on;

colormap(flipud(hot))
