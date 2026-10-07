%% Paired whole-cell patch-clamp analysis: spikes, IPSPs and cross-cell alignment
%
% Analyses a dual whole-cell current-clamp recording (.abf) of two
% simultaneously recorded neurons. For each cell the script:
%   1. removes 50 Hz mains noise with a zero-phase Butterworth band-stop filter
%   2. detects spikes and IPSPs
%   3. aligns events in one cell to events in the other cell
%   4. tests spike-amplitude and ISI distributions for unimodality
%      (Hartigan's dip test)
%   5. saves figures (.svg) and per-cell results (.mat)
%
% USAGE
%   Run from a folder containing the .abf recording. The membrane potential
%   of Cell 1 and Cell 2 is read from columns 1 and 3 of the ABF data; other
%   columns are not used.
%
% OUTPUTS
%   Cell1Data.mat (struct S1) and Cell2Data.mat (struct S2), saved in the
%   working folder and copied to ./figures. See "Save per-cell results"
%   below for field definitions.
%   ./figures/*.svg
%     Cell1IPSPstat, Cell2IPSPstat   IPSP detection summaries
%     SpikeAlignment                 spike-triggered spike alignment
%     AllIPSPAlignment               IPSP-triggered IPSP alignment
%     SpikestoIPSPalign              IPSP times relative to spikes
%     IPSPtospikesalign              spike times relative to IPSPs
%     Cell1SpikeStat, Cell2SpikeStat spike amplitude / ISI distributions
%
% REQUIREMENTS
%   MATLAB R2020b or later (swarmchart)
%   Signal Processing Toolbox (designfilt, filtfilt, findpeaks)
%   Statistics and Machine Learning Toolbox (boxplot, in IPSPdetect)
%   Curve Fitting Toolbox (smooth, in complspikeoverlay)
%   abfload       - ABF file reader (H. Hentschke, MATLAB File Exchange)
%   dipTest       - Hartigan's dip test of unimodality
%                   (Hartigan & Hartigan, 1985, Ann. Stat. 13:70-84)
%   spikecounter5_0, IPSPdetect, complspikeoverlay, darken_hex_color
%                 - analysis functions included with this code
%
% Author: Annie Park, Waddell Lab, Centre for Neural Circuits and Behaviour,
%         University of Oxford

clear
close all

%% Parameters
fs            = 10000;   % sampling rate (Hz)
windowsz      = 5000;    % event window passed to detection/alignment functions (samples)
diptestmin    = 30;      % run dip tests only if either cell has more spikes than this
vmStartSample = 1000;    % first sample used when averaging filtVm

% Plot colours (dark greens for MBONs on a white background)
color1 = "#175340";      % Cell 1
color2 = "#0e9168";      % Cell 2

% Output folder
figDir = fullfile(pwd, 'figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% Load recordings and remove 50 Hz mains noise
% 2nd-order Butterworth band-stop (49-51 Hz), applied with filtfilt (zero phase).
notchFilt = designfilt('bandstopiir', 'FilterOrder', 2, ...
    'HalfPowerFrequency1', 49, 'HalfPowerFrequency2', 51, ...
    'DesignMethod', 'butter', 'SampleRate', fs);

filelist = dir('*.abf');
filename = {filelist.name};
nFiles   = numel(filename);

% NOTE: every .abf file in the folder is loaded and filtered, but only the
% first file (filename{1}) is analysed below.
rawData  = cell(1, nFiles);
filtData = cell(1, nFiles);
for iFile = 1:nFiles
    rawData{iFile} = abfload(filename{iFile});
    filtData{iFile}(:,1) = filtfilt(notchFilt, rawData{iFile}(:,1));   % Cell 1 Vm
    filtData{iFile}(:,3) = filtfilt(notchFilt, rawData{iFile}(:,3));   % Cell 2 Vm
end

vm1 = filtData{1}(:,1);   % Cell 1, notch-filtered membrane potential
vm2 = filtData{1}(:,3);   % Cell 2, notch-filtered membrane potential

%% Spike detection
% Each call detects spikes in one cell. The other cell's trace is also
% passed so that the complementary (other-cell) segments around each spike
% are returned. ISIs are recomputed below.
[spikeidx1, spikeamp1, filtVm1, thresh_peak1, overlayspike1, spss1, binaryspk1, ...
    complcellsp2, binaryspkoverlay1, stval1, spnum1, subtrnum1, ~, answer_celltype1] = ...
    spikecounter5_0(vm1, fs, spikeisimin, vm2, windowsz);

[spikeidx2, spikeamp2, filtVm2, thresh_peak2, overlayspike2, spss2, binaryspk2, ...
    complcellsp1, binaryspkoverlay2, stval2, spnum2, subtrnum2, ~, answer_celltype2] = ...
    spikecounter5_0(vm2, fs, spikeisimin, vm1, windowsz);

spikeidx1 = spikeidx1.';
spikeidx2 = spikeidx2.';
spikeisi1 = diff(spikeidx1) ./ fs;   % inter-spike intervals (s)
spikeisi2 = diff(spikeidx2) ./ fs;

%% IPSP detection
% As for spikes, the other cell's trace is passed so that its segments
% around each IPSP are returned (complcellipsp).
[ipspamp1, ipspidx1, overlayipsp1, binaryipsp1, complcellipsp2, thresh_valipsp1] = ...
    IPSPdetect(vm1, fs, filtVm1, vm2, windowsz);
saveas(gcf, fullfile(figDir, 'Cell1IPSPstat.svg'));

[ipspamp2, ipspidx2, overlayipsp2, binaryipsp2, complcellipsp1, thresh_valipsp2] = ...
    IPSPdetect(vm2, fs, filtVm2, vm1, windowsz);
saveas(gcf, fullfile(figDir, 'Cell2IPSPstat.svg'));

ipspipspint1 = diff(ipspidx1);   % inter-IPSP intervals (samples)
ipspipspint2 = diff(ipspidx2);

% Binary vectors with a 1 at each IPSP sample index (all IPSPs)
binaryipspall1 = zeros(size(ipspidx1));
binaryipspall1(ipspidx1) = 1;
binaryipspall2 = zeros(size(ipspidx2));
binaryipspall2(ipspidx2) = 1;

%% Cross-cell event alignment
% Spike-triggered: spike times relative to each spike
[complbinspkoverlay1, complbinspkoverlay2] = complspikeoverlay( ...
    overlayspike1, overlayspike2, complcellsp1, complcellsp2, ...
    spikeidx1, spikeidx2, binaryspk1, binaryspk2, windowsz, color1, color2, fs);
saveas(gcf, fullfile(figDir, 'SpikeAlignment.svg'));

% IPSP-triggered: IPSP times relative to each IPSP
[complipspoverlay1, complipspoverlay2] = complspikeoverlay( ...
    overlayipsp1, overlayipsp2, complcellipsp1, complcellipsp2, ...
    ipspidx1, ipspidx2, binaryipsp1, binaryipsp2, windowsz, color1, color2, fs);
saveas(gcf, fullfile(figDir, 'AllIPSPAlignment.svg'));

% Spike-triggered: IPSP times (all IPSPs) relative to each spike
[complbinspkipspoverlay1, complbinspkipspoverlay2] = complspikeoverlay( ...
    overlayspike1, overlayspike2, complcellsp1, complcellsp2, ...
    spikeidx1, spikeidx2, binaryipspall1, binaryipspall2, windowsz, color1, color2, fs);
saveas(gcf, fullfile(figDir, 'SpikestoIPSPalign.svg'));

% IPSP-triggered: spike times relative to each IPSP
[complbinipspspkoverlay1, complbinipspspkoverlay2] = complspikeoverlay( ...
    overlayipsp1, overlayipsp2, complcellipsp1, complcellipsp2, ...
    ipspidx1, ipspidx2, binaryspk1, binaryspk2, windowsz, color1, color2, fs);
saveas(gcf, fullfile(figDir, 'IPSPtospikesalign.svg'));

%% Unimodality of spike amplitude and ISI distributions (Hartigan's dip test)
if max([numel(spikeamp1), numel(spikeamp2)]) > diptestmin
    [pspamp1, dipspamp1, psisi1, dipisi1] = plotSpikeStats(spikeamp1, spikeisi1, 'Cell1');
    saveas(gcf, fullfile(figDir, 'Cell1SpikeStat.svg'));

    [pspamp2, dipspamp2, psisi2, dipisi2] = plotSpikeStats(spikeamp2, spikeisi2, 'Cell2');
    saveas(gcf, fullfile(figDir, 'Cell2SpikeStat.svg'));
end

%% Save per-cell results
% Both cells come from the same recording, so both structs share a file name.
S1.name          = filename{1};
S1.spikeidx      = spikeidx1;                                    % spike times (samples)
S1.spiketimesec  = spikeidx1 ./ fs;                              % spike times (s)
S1.spikeamp      = spikeamp1;                                    % spike amplitudes (mV): notch-filtered trace minus filtVm
S1.spikeisi      = spikeisi1;                                    % inter-spike intervals (s)
S1.ipspidx       = ipspidx1;                                     % IPSP times (samples)
S1.ipsptimesec   = ipspidx1 ./ fs;                               % IPSP times (s)
S1.ipspamp       = ipspamp1;                                     % IPSP amplitudes (mV): filtVm minus notch-filtered trace (not low-passed)
S1.ipspisi       = ipspipspint1 ./ fs;                           % inter-IPSP intervals (s)
S1.threshpeak    = thresh_peak1;                                 % spike detection threshold
S1.threshipsp    = thresh_valipsp1;                              % IPSP detection threshold
S1.filtVm        = mean(filtVm1(vmStartSample:end), 'omitnan');  % mean baseline (1 s median filter of notch-filtered trace)
S1.Vm            = mean(vm1);                                    % mean of notch-filtered trace (not low-passed)
save('Cell1Data.mat', 'S1');
copyfile('Cell1Data.mat', figDir);

S2.name          = filename{1};
S2.spikeidx      = spikeidx2;
S2.spiketimesec  = spikeidx2 ./ fs;
S2.spikeamp      = spikeamp2;
S2.spikeisi      = spikeisi2;
S2.ipspidx       = ipspidx2;
S2.ipsptimesec   = ipspidx2 ./ fs;
S2.ipspamp       = ipspamp2;
S2.ipspisi       = ipspipspint2 ./ fs;
S2.threshpeak    = thresh_peak2;
S2.threshipsp    = thresh_valipsp2;
S2.filtVm        = mean(filtVm2(vmStartSample:end), 'omitnan');
S2.Vm            = mean(vm2);
save('Cell2Data.mat', 'S2');
copyfile('Cell2Data.mat', figDir);


%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function [pAmp, dipAmp, pISI, dipISI] = plotSpikeStats(spikeamp, spikeisi, cellLabel)
% PLOTSPIKESTATS  Hartigan's dip test and swarm/histogram plots of spike
% amplitude and inter-spike interval for one cell.

[pAmp, dipAmp, ~, ~] = dipTest(spikeamp);
[pISI, dipISI, ~, ~] = dipTest(spikeisi);

x = ones(1, numel(spikeamp));

figure
subplot(1,4,1)
swarmchart(x, spikeamp, 'k')
ylabel('Spike amplitude (mV)')
set(gcf, 'Visible', 'on')
title(sprintf('Dip test p = %.2f, dip = %.2f', pAmp, dipAmp));

subplot(1,4,2)
histogram(spikeamp, 'FaceColor', 'k');
camroll(90)

subplot(1,4,3)
swarmchart(x(2:end), spikeisi, 'k');
ylabel('Inter-spike interval (s)');
title(sprintf('Dip test p = %.2f, dip = %.2f', pISI, dipISI));

subplot(1,4,4)
histogram(spikeisi, 'FaceColor', 'k');
camroll(90)
title(cellLabel);
end
