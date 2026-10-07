function analyze_bursts_dual_gui_isi(file1, file2, varargin)
% ANALYZE_BURSTS_DUAL_GUI_ISI
% Interactive GUI to pick an ISI threshold (slider + numeric entry) and
% detect bursts in two spike .mat files, overlay bursts on raw traces
% (reads .abf if present via abfload).
%
% Usage:
%   analyze_bursts_dual_gui_isi('Cell1Data.mat','Cell2Data.mat');
% Options:
%   'UseABF'       (true/false) default true
%   'ABFfile1'     explicit ABF filename for file1 folder (optional)
%   'ABFfile2'     explicit ABF filename for file2 folder (optional)
%   'ABFchannel1'  channel number for file1 ABF (default 1)
%   'ABFchannel2'  channel number for file2 ABF (default 3)

% ---------------- Parse inputs ----------------
p = inputParser;
addRequired(p,'file1',@ischar);
addRequired(p,'file2',@ischar);
addParameter(p,'UseABF',true,@islogical);
addParameter(p,'ABFfile1','',@ischar);
addParameter(p,'ABFfile2','',@ischar);
addParameter(p,'ABFchannel1',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'ABFchannel2',3,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'StrictBoth',false,@islogical);
params.StrictBoth = false;


parse(p,file1,file2,varargin{:});

useABF = p.Results.UseABF;
ABFfile1 = p.Results.ABFfile1;
ABFfile2 = p.Results.ABFfile2;
abfChan1 = p.Results.ABFchannel1;
abfChan2 = p.Results.ABFchannel2;

% ---------------- Load spikes and raw ----------------
[S1, times1, spikeidx1] = load_struct_and_spikes(file1,'S1');
[S2, times2, spikeidx2] = load_struct_and_spikes(file2,'S2');
times1 = unique(sort(times1));
times2 = unique(sort(times2));

[raw1, Fs1, rawName1] = try_load_abf_or_mat_raw(file1, ABFfile1, abfChan1, useABF, S1);
[raw2, Fs2, rawName2] = try_load_abf_or_mat_raw(file2, ABFfile2, abfChan2, useABF, S2);

% ---------------- Initial parameters ----------------
params.BurstISI = 0.08;    % seconds default
params.MinSpikes = 3;
params.CoincWin = 0.05;
params.StrictBoth = false;

% window: 10s around center of data if possible
if ~isempty(times1) || ~isempty(times2)
    allT = [times1(:); times2(:)];
    if isempty(allT)
        params.PlotStart = 0;
    else
        mid = (min(allT) + max(allT))/2;
        params.PlotStart = max(0, mid - 5);
    end
else
    params.PlotStart = 0;
end
params.PlotDuration = 10;

% slider bounds
smin = 0.001; smax = 0.5;

% ---------------- compute ISIs ----------------
[ISI1, ISI1_t, ISI1_idx] = compute_isi_values(times1);
[ISI2, ISI2_t, ISI2_idx] = compute_isi_values(times2);

% initial detection over full trains
[bursts1, spikeMask1] = detect_bursts_threshold(times1, params.BurstISI, params.MinSpikes);
[bursts2, spikeMask2] = detect_bursts_threshold(times2, params.BurstISI, params.MinSpikes);
coincSummary = match_coincidence(bursts1, bursts2, params.CoincWin);

% ---------------- Build GUI layout ----------------
fig = figure('Name','Burst GUI — ISI slider','Units','normalized','Position',[0.05 0.05 0.88 0.84],...
    'Color',[1 1 1],'NumberTitle','off','Resize','on');

% Left plotting area (reduced width to give controls room)
plotLeft = 0.03; plotWidth = 0.58;
axRaw1 = axes('Parent',fig,'Units','normalized','Position',[plotLeft 0.69 plotWidth 0.26]);
axRaw2 = axes('Parent',fig,'Units','normalized','Position',[plotLeft 0.41 plotWidth 0.24]);
axRaster = axes('Parent',fig,'Units','normalized','Position',[plotLeft 0.09 plotWidth 0.28]);

% Right-side control panel
ctrlLeft = plotLeft + plotWidth + 0.03; ctrlWidth = 0.34;
ctrlBottom = 0.05; ctrlTop = 0.94; ctrlHeight = ctrlTop - ctrlBottom;
ctrlPanel = uipanel('Parent',fig,'Units','normalized','Position',[ctrlLeft ctrlBottom ctrlWidth ctrlHeight],'Title','Controls','FontWeight','bold');

% ISI axes inside panel (top)
panelPad = 0.04;
axISI_pos = [panelPad 0.58 1-2*panelPad 0.36];
axISI = axes('Parent',ctrlPanel,'Units','normalized','Position',axISI_pos);
title(axISI,'Inter-spike interval (ISI)');

% Slider row under ISI axes: slider + numeric edit (no big label)
sliderH = 0.065;
sliderY = axISI_pos(2) - sliderH - 0.02;
sldThresh = uicontrol('Parent',ctrlPanel,'Style','slider','Units','normalized', ...
    'Position',[panelPad sliderY 0.92 sliderH], 'Min',smin,'Max',smax,'Value',params.BurstISI,'Callback',@cb_onThreshSlider);
editW = 0.28;
editThresh = uicontrol('Parent',ctrlPanel,'Style','edit','Units','normalized',...
    'Position',[panelPad+0.64 sliderY editW sliderH], 'String',num2str(params.BurstISI,'%.4f'),'Callback',@cb_onThreshEdit,'KeyPressFcn',@cb_onEditKey,'BackgroundColor','white');

% Controls stacked below slider
lineH = 0.055; gap = 0.03;
y = sliderY - (lineH + 0.02);

uicontrol('Parent',ctrlPanel,'Style','text','Units','normalized','Position',[panelPad y+lineH*0 0.24 lineH],...
    'String','ISI plot:','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
popupISI = uicontrol('Parent',ctrlPanel,'Style','popupmenu','Units','normalized','Position',[panelPad+0.24 y+lineH*0 0.72 lineH],...
    'String',{'Trace 1','Trace 2'},'Callback',@cb_onChangeISItrace);

y = y - (lineH + gap);
chkBoth = uicontrol('Parent',ctrlPanel,'Style','checkbox','Units','normalized','Position',[panelPad y 0.95 lineH],...
    'String','(compat) Require prev OR next ISI < thresh (default)','Value',1,'BackgroundColor',[1 1 1],'Enable','off'); %#ok<NASGU>

y = y - (lineH + gap);
uicontrol('Parent',ctrlPanel,'Style','text','Units','normalized','Position',[panelPad y 0.45 lineH],...
    'String','Min spikes/burst','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
editMinSp = uicontrol('Parent',ctrlPanel,'Style','edit','Units','normalized','Position',[panelPad+0.47 y 0.53 lineH],...
    'String',num2str(params.MinSpikes),'Callback',@cb_onParamEdit,'BackgroundColor','white');

y = y - (lineH + gap);
uicontrol('Parent',ctrlPanel,'Style','text','Units','normalized','Position',[panelPad y 0.45 lineH],...
    'String','Coinc window (s)','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
editCoinc = uicontrol('Parent',ctrlPanel,'Style','edit','Units','normalized','Position',[panelPad+0.47 y 0.53 lineH],...
    'String',num2str(params.CoincWin),'Callback',@cb_onParamEdit,'BackgroundColor','white');

% Buttons row (big) - Apply, Analyze & Save (Full)
y = y - (lineH + 1.3*gap);
btnH = lineH + 0.01;
btnApply = uicontrol('Parent',ctrlPanel,'Style','pushbutton','Units','normalized','Position',[panelPad y 0.46 btnH],...
    'String','Apply & Update (window)','FontWeight','bold','Callback',@cb_onApply);
btnFull = uicontrol('Parent',ctrlPanel,'Style','pushbutton','Units','normalized','Position',[panelPad+0.49 y 0.46 btnH],...
    'String','Analyze & Save (Full)','FontWeight','bold','BackgroundColor',[0.85 0.95 1],'Callback',@cb_onAnalyzeFull);

% Save current results small
y = y - (btnH + gap);
btnSave = uicontrol('Parent',ctrlPanel,'Style','pushbutton','Units','normalized','Position',[panelPad y 0.95 btnH],...
    'String','Save Current Results (.mat)','Callback',@cb_onSave);

% Stats listbox: compute height safely (prevent negative)
y = y - (btnH + 1.2*gap);
minListboxHeight = 0.06;
listboxHeight = max(minListboxHeight, y - 0.03);
lstStats = uicontrol('Parent',ctrlPanel,'Style','listbox','Units','normalized','Position',[panelPad 0.03 0.94 listboxHeight],...
    'String',{'Initializing...'},'BackgroundColor',[1 1 1],'FontName','FixedWidth');

% store handles & data in appdata
gui.h.fig = fig;
gui.h.axRaw1 = axRaw1; gui.h.axRaw2 = axRaw2; gui.h.axRaster = axRaster; gui.h.axISI = axISI;
gui.h.sldThresh = sldThresh; gui.h.editThresh = editThresh; gui.h.popupISI = popupISI;
gui.h.editMinSp = editMinSp; gui.h.editCoinc = editCoinc;
gui.h.btnApply = btnApply; gui.h.btnFull = btnFull; gui.h.btnSave = btnSave;
gui.h.lstStats = lstStats;
gui.file1 = file1; gui.file2 = file2;
gui.times1 = times1; gui.times2 = times2;
gui.ISI1 = ISI1; gui.ISI1_t = ISI1_t; gui.ISI1_idx = ISI1_idx;
gui.ISI2 = ISI2; gui.ISI2_t = ISI2_t; gui.ISI2_idx = ISI2_idx;
gui.raw1 = raw1; gui.raw2 = raw2; gui.Fs1 = Fs1; gui.Fs2 = Fs2;
gui.rawName1 = rawName1; gui.rawName2 = rawName2;
gui.params = params;
gui.bursts1 = bursts1; gui.bursts2 = bursts2; gui.coincSummary = coincSummary;
gui.spikeMask1 = spikeMask1; gui.spikeMask2 = spikeMask2;
gui.activeISI = 1;
setappdata(fig,'gui',gui);

% initial draw
draw_all(fig);

% ----------------- NESTED CALLBACKS -----------------
% These are nested so they can access GUI via getappdata/setappdata.
    function cb_onThreshSlider(src, ~) %#ok<DEFNU>
        figLocal = ancestor(src,'figure');
        guiLocal = getappdata(figLocal,'gui');
        val = get(guiLocal.h.sldThresh,'Value');
        set(guiLocal.h.editThresh,'String',num2str(val,'%.4f'));
        guiLocal.params.BurstISI = val;
        setappdata(figLocal,'gui',guiLocal);
        apply_threshold_and_update(figLocal);
        draw_all(figLocal);
    end

    function cb_onThreshEdit(src, ~) %#ok<DEFNU>
        figLocal = ancestor(src,'figure');
        guiLocal = getappdata(figLocal,'gui');
        val = str2double(get(guiLocal.h.editThresh,'String'));
        if isnan(val)
            errordlg('Enter numeric ISI threshold (seconds).','Invalid input');
            set(guiLocal.h.editThresh,'String',num2str(guiLocal.params.BurstISI,'%.4f'));
            return;
        end
        val = max(min(val,smax),smin);
        set(guiLocal.h.sldThresh,'Value',val);
        set(guiLocal.h.editThresh,'String',num2str(val,'%.4f'));
        guiLocal.params.BurstISI = val;
        setappdata(figLocal,'gui',guiLocal);
        apply_threshold_and_update(figLocal);
        draw_all(figLocal);
    end

    function cb_onEditKey(src, ev) %#ok<INUSL,DEFNU>
        if isequal(ev.Key,'return') || isequal(ev.Key,'enter')
            cb_onThreshEdit(src,[]);
        end
    end

    function cb_onChangeISItrace(src, ~) %#ok<DEFNU>
        figLocal = ancestor(src,'figure'); guiLocal = getappdata(figLocal,'gui');
        guiLocal.activeISI = get(guiLocal.h.popupISI,'Value');
        setappdata(figLocal,'gui',guiLocal);
        draw_ISI(figLocal, guiLocal);
    end

    function cb_onParamEdit(~, ~) %#ok<DEFNU>
        figLocal = gcbf; guiLocal = getappdata(figLocal,'gui');
        v = str2double(get(guiLocal.h.editMinSp,'String'));
        if ~isnan(v) && v>=1, guiLocal.params.MinSpikes = round(v); end
        v2 = str2double(get(guiLocal.h.editCoinc,'String'));
        if ~isnan(v2) && v2>=0, guiLocal.params.CoincWin = v2; end
        setappdata(figLocal,'gui',guiLocal);
        apply_threshold_and_update(figLocal);
        draw_all(figLocal);
    end

    function cb_onApply(~, ~) %#ok<DEFNU>
        figLocal = gcbf;
        apply_threshold_and_update(figLocal);
        draw_all(figLocal);
    end

    function cb_onSave(~, ~) %#ok<DEFNU>
        figLocal = gcbf; guiLocal = getappdata(figLocal,'gui');
        results = struct();
        results.file1 = guiLocal.file1; results.file2 = guiLocal.file2;
        results.params = guiLocal.params;
        results.bursts1 = guiLocal.bursts1; results.bursts2 = guiLocal.bursts2;
        results.spikeMask1 = guiLocal.spikeMask1; results.spikeMask2 = guiLocal.spikeMask2;
        results.coincidence = guiLocal.coincSummary;
        results.rawName1 = guiLocal.rawName1; results.rawName2 = guiLocal.rawName2;
        results.times1 = guiLocal.times1; results.times2 = guiLocal.times2;
        save('burst_current_results.mat','results','-v7.3');
        msgbox('Current results saved to burst_current_results.mat','Saved');
    end

    function cb_onAnalyzeFull(~, ~) %#ok<DEFNU>
        % Robust full-analysis callback with try/catch and error logging.
        figLocal = gcbf;
        guiLocal = getappdata(figLocal,'gui');

        try
            p = guiLocal.params;

            % Run detection on whole spike trains
            [fullBursts1, fullMask1] = detect_bursts_threshold(guiLocal.times1, p.BurstISI, p.MinSpikes);
            [fullBursts2, fullMask2] = detect_bursts_threshold(guiLocal.times2, p.BurstISI, p.MinSpikes);

            % Match coincidence between detected bursts
            fullCoinc = match_coincidence(fullBursts1, fullBursts2, p.CoincWin);

            % Compute per-burst stats
            stats1 = compute_burst_stats(fullBursts1);
            stats2 = compute_burst_stats(fullBursts2);

            % Use file-local helper for recording duration
            recDur1 = recording_duration(guiLocal.times1, guiLocal.raw1, guiLocal.Fs1);
            recDur2 = recording_duration(guiLocal.times2, guiLocal.raw2, guiLocal.Fs2);

            % --- Build summaries for cell1 ---
            summary1 = struct();
            summary1.nBursts = numel(fullBursts1);
            if isempty(stats1) || ~isfield(stats1,'spikesPerBurst'), stats1.spikesPerBurst = []; end
            if isempty(stats1.spikesPerBurst)
                summary1.meanSpikesPerBurst = NaN; summary1.medianSpikesPerBurst = NaN; summary1.stdSpikesPerBurst = NaN;
            else
                summary1.meanSpikesPerBurst = mean(stats1.spikesPerBurst);
                summary1.medianSpikesPerBurst = median(stats1.spikesPerBurst);
                summary1.stdSpikesPerBurst  = std(stats1.spikesPerBurst);
            end
            if isempty(stats1) || ~isfield(stats1,'durations'), stats1.durations = []; end
            if isempty(stats1.durations)
                summary1.meanBurstDuration = NaN; summary1.medianBurstDuration = NaN; summary1.stdBurstDuration = NaN;
            else
                summary1.meanBurstDuration = mean(stats1.durations);
                summary1.medianBurstDuration = median(stats1.durations);
                summary1.stdBurstDuration  = std(stats1.durations);
            end
            if isempty(stats1) || ~isfield(stats1,'periods'), stats1.periods = []; end
            if isempty(stats1.periods)
                summary1.meanPeriod = NaN; summary1.medianPeriod = NaN; summary1.stdPeriod = NaN;
            else
                summary1.meanPeriod = mean(stats1.periods);
                summary1.medianPeriod = median(stats1.periods);
                summary1.stdPeriod = std(stats1.periods);
            end
            summary1.meanInterBurstInterval = summary1.meanPeriod;

            if ~isnan(recDur1) && recDur1>0
                summary1.burstsPerSecond = summary1.nBursts / recDur1;
            else
                summary1.burstsPerSecond = NaN;
            end

            totalSpikes1 = numel(guiLocal.times1);
            if isempty(fullMask1), fullMask1 = false(totalSpikes1,1); end
            nBurstSpikes1 = sum(fullMask1(:));
            nTonicSpikes1 = totalSpikes1 - nBurstSpikes1;
            summary1.nTonicSpikes = nTonicSpikes1;
            if ~isnan(recDur1) && recDur1>0
                summary1.tonicSpikeRate_Hz = nTonicSpikes1 / recDur1;
            else
                summary1.tonicSpikeRate_Hz = NaN;
            end
            if isfield(fullCoinc,'nPairs'), summary1.nCoincidentPairs = fullCoinc.nPairs; else summary1.nCoincidentPairs = 0; end
            if isfield(fullCoinc,'pctBursts1Coincident'), summary1.pctBurstsCoincident = fullCoinc.pctBursts1Coincident; else summary1.pctBurstsCoincident = 0; end
            summary1.recordingDuration_s = recDur1;

            % Units
            summary1.units = struct('nBursts','count','spikesPerBurst','spikes','meanBurstDuration','s',...
                                   'period','s','burstsPerSecond','Hz','nTonicSpikes','spikes',...
                                   'tonicSpikeRate_Hz','Hz','recordingDuration_s','s','pctBurstsCoincident','percent');

            % --- Build summaries for cell2 (same as cell1) ---
            summary2 = struct();
            summary2.nBursts = numel(fullBursts2);
            if isempty(stats2) || ~isfield(stats2,'spikesPerBurst'), stats2.spikesPerBurst = []; end
            if isempty(stats2.spikesPerBurst)
                summary2.meanSpikesPerBurst = NaN; summary2.medianSpikesPerBurst = NaN; summary2.stdSpikesPerBurst = NaN;
            else
                summary2.meanSpikesPerBurst = mean(stats2.spikesPerBurst);
                summary2.medianSpikesPerBurst = median(stats2.spikesPerBurst);
                summary2.stdSpikesPerBurst  = std(stats2.spikesPerBurst);
            end
            if isempty(stats2) || ~isfield(stats2,'durations'), stats2.durations = []; end
            if isempty(stats2.durations)
                summary2.meanBurstDuration = NaN; summary2.medianBurstDuration = NaN; summary2.stdBurstDuration = NaN;
            else
                summary2.meanBurstDuration = mean(stats2.durations);
                summary2.medianBurstDuration = median(stats2.durations);
                summary2.stdBurstDuration  = std(stats2.durations);
            end
            if isempty(stats2) || ~isfield(stats2,'periods'), stats2.periods = []; end
            if isempty(stats2.periods)
                summary2.meanPeriod = NaN; summary2.medianPeriod = NaN; summary2.stdPeriod = NaN;
            else
                summary2.meanPeriod = mean(stats2.periods);
                summary2.medianPeriod = median(stats2.periods);
                summary2.stdPeriod = std(stats2.periods);
            end
            summary2.meanInterBurstInterval = summary2.meanPeriod;

            if ~isnan(recDur2) && recDur2>0
                summary2.burstsPerSecond = summary2.nBursts / recDur2;
            else
                summary2.burstsPerSecond = NaN;
            end

            totalSpikes2 = numel(guiLocal.times2);
            if isempty(fullMask2), fullMask2 = false(totalSpikes2,1); end
            nBurstSpikes2 = sum(fullMask2(:));
            nTonicSpikes2 = totalSpikes2 - nBurstSpikes2;
            summary2.nTonicSpikes = nTonicSpikes2;
            if ~isnan(recDur2) && recDur2>0
                summary2.tonicSpikeRate_Hz = nTonicSpikes2 / recDur2;
            else
                summary2.tonicSpikeRate_Hz = NaN;
            end
            if isfield(fullCoinc,'nPairs'), summary2.nCoincidentPairs = fullCoinc.nPairs; else summary2.nCoincidentPairs = 0; end
            if isfield(fullCoinc,'pctBursts2Coincident'), summary2.pctBurstsCoincident = fullCoinc.pctBursts2Coincident; else summary2.pctBurstsCoincident = 0; end
            summary2.recordingDuration_s = recDur2;
            summary2.units = summary1.units;

            % ------------------ Build results struct ------------------
            results = struct();
            results.file1 = guiLocal.file1;
            results.file2 = guiLocal.file2;
            results.params = p;
            results.bursts1 = fullBursts1;
            results.bursts2 = fullBursts2;
            results.spikeMask1 = fullMask1;
            results.spikeMask2 = fullMask2;
            results.coincidence = fullCoinc;
            results.stats1 = stats1;
            results.stats2 = stats2;
            results.summary.cell1 = summary1;
            results.summary.cell2 = summary2;
            results.summary.units = summary1.units;

            % -------- threshold and provenance metadata --------
            results.threshold_used = p.BurstISI;
            results.threshold_strict = logical(p.StrictBoth);
            results.threshold_minSpikes = p.MinSpikes;
            results.threshold_coincWindow = p.CoincWin;

            % detect if slider/edit used (best-effort)
            th_source = 'gui';
            try
                if isfield(guiLocal.h,'sldThresh') && ishandle(guiLocal.h.sldThresh)
                    sval = get(guiLocal.h.sldThresh,'Value');
                    if abs(sval - p.BurstISI) < max(1e-9,1e-6*abs(p.BurstISI)), th_source = 'slider'; end
                end
                if isfield(guiLocal.h,'editThresh') && ishandle(guiLocal.h.editThresh)
                    estr = get(guiLocal.h.editThresh,'String');
                    if ~isempty(estr) && abs(str2double(estr) - p.BurstISI) < max(1e-9,1e-6*abs(p.BurstISI)), th_source = 'edit'; end
                end
            catch
                % ignore
            end
            results.threshold_source = th_source;
            results.threshold_timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');

            results.rawName1 = guiLocal.rawName1;
            results.rawName2 = guiLocal.rawName2;
            results.times1 = guiLocal.times1;
            results.times2 = guiLocal.times2;
            results.generatedOn = datestr(now,'yyyy-mm-dd HH:MM:SS');

            % Save results atomically: write to a temp file then move/rename
            tmpname = ['tmp_' char(java.util.UUID.randomUUID) '.mat']; %#ok<CHAR>
            save(tmpname,'results','-v7.3');
            % move to final name (overwrite)
            try
                if exist('burst_full_results.mat','file'), delete('burst_full_results.mat'); end
                movefile(tmpname,'burst_full_results.mat');
            catch
                % fallback: save directly if move fails
                save('burst_full_results.mat','results','-v7.3');
                if exist(tmpname,'file'), delete(tmpname); end
            end

            % Update GUI state and redraw
            guiLocal.bursts1 = fullBursts1;
            guiLocal.bursts2 = fullBursts2;
            guiLocal.spikeMask1 = fullMask1;
            guiLocal.spikeMask2 = fullMask2;
            guiLocal.coincSummary = fullCoinc;
            setappdata(figLocal,'gui',guiLocal);
            draw_all(figLocal);

            % Inform user
            msg = sprintf(['Full analysis saved to burst_full_results.mat\n\n' ...
                'Cell1 bursts: %d   Cell2 bursts: %d\n' ...
                'Cell1 bursts/sec: %.4g Hz   Cell2 bursts/sec: %.4g Hz\n' ...
                'Cell1 tonic spikes: %d   tonic rate: %.4g Hz\n' ...
                'Cell2 tonic spikes: %d   tonic rate: %.4g Hz\n\n' ...
                'Threshold used: %.4g s   StrictBoth: %d   MinSpikes: %d   CoincWin: %.4g s\n' ...
                'Threshold source: %s'], ...
                summary1.nBursts, summary2.nBursts, summary1.burstsPerSecond, summary2.burstsPerSecond, ...
                summary1.nTonicSpikes, summary1.tonicSpikeRate_Hz, summary2.nTonicSpikes, summary2.tonicSpikeRate_Hz, ...
                results.threshold_used, results.threshold_strict, results.threshold_minSpikes, results.threshold_coincWindow, results.threshold_source);
            msgbox(msg,'Full analysis saved');

        catch ME
            % ---------------- ERROR HANDLING ----------------
            errMsg = sprintf('Error: %s\n\nIn: %s at line %d', ME.message, ME.stack(1).file, ME.stack(1).line);
            try
                errordlg(errMsg,'Analysis error','modal');
            catch
                warning('Error in cb_onAnalyzeFull: %s', ME.message);
            end

            % Save error details to a log (folder-local)
            try
                errlog.results = [];
                errlog.ME = ME;
                errlog.timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
                save('error_log.mat','errlog','-v7.3');

                % also write plain text for quick inspection
                fid = fopen('error_log.txt','w');
                if fid ~= -1
                    fprintf(fid,'Error occurred at %s\n\n', errlog.timestamp);
                    fprintf(fid,'Message:\n%s\n\n', ME.message);
                    fprintf(fid,'Stack (top 5 frames):\n');
                    for si = 1:min(5,numel(ME.stack))
                        fprintf(fid,'%d) File: %s  Function: %s  Line: %d\n', si, ME.stack(si).file, ME.stack(si).name, ME.stack(si).line);
                    end
                    fclose(fid);
                end
            catch
                % ignore logging errors
            end
        end
    end % cb_onAnalyzeFull

end % main function end


% ----------------- Helper routines (file-local) -----------------
%%
function dur = recording_duration(times, raw, Fs)
    dur = NaN;
    if ~isempty(times) && numel(times) > 1
        dur = max(times) - min(times);
    end
    if ~(isfinite(dur) && dur>0)
        if ~isempty(raw) && ~isempty(Fs) && Fs>0
            dur = numel(raw)/Fs;
        end
    end
    if ~(isfinite(dur) && dur>0)
        dur = NaN;
    end
end

function apply_threshold_and_update(figHandle)
    if nargin<1, figHandle = gcbf; end
    guiL = getappdata(figHandle,'gui');
    [guiL.bursts1, guiL.spikeMask1] = detect_bursts_threshold(guiL.times1, guiL.params.BurstISI, guiL.params.MinSpikes);
    [guiL.bursts2, guiL.spikeMask2] = detect_bursts_threshold(guiL.times2, guiL.params.BurstISI, guiL.params.MinSpikes);
    guiL.coincSummary = match_coincidence(guiL.bursts1, guiL.bursts2, guiL.params.CoincWin);
    setappdata(figHandle,'gui',guiL);
end

function draw_all(figHandle)
    guiL = getappdata(figHandle,'gui');
    draw_raw_axes(figHandle, guiL);
    draw_ISI(figHandle, guiL);
    draw_raster(figHandle, guiL);
    update_stats_list(figHandle, guiL);
end

function draw_raw_axes(figHandle, guiL)
    h = guiL.h;
    % Raw1
    axes(h.axRaw1); cla;
    if ~isempty(guiL.raw1)
        Fs = guiL.Fs1; if isempty(Fs), Fs = 1; end
        t_raw = (0:numel(guiL.raw1)-1)/Fs;
        idx = t_raw >= guiL.params.PlotStart & t_raw <= (guiL.params.PlotStart + guiL.params.PlotDuration);
        if any(idx)
            plot(t_raw(idx), guiL.raw1(idx),'k-','LineWidth',0.6); hold on;
            sp = guiL.times1(guiL.times1>=guiL.params.PlotStart & guiL.times1<=guiL.params.PlotStart+guiL.params.PlotDuration);
            yl = ylim();
            for s = sp(:)'; plot([s s],[yl(1)+0.02*diff(yl) yl(1)+0.12*diff(yl)],'r'); end
            for b = guiL.bursts1
                a = max(b.start, guiL.params.PlotStart); bE = min(b.end, guiL.params.PlotStart+guiL.params.PlotDuration);
                if a < bE
                    patch([a bE bE a],[yl(1) yl(1) yl(2) yl(2)],[0 0.6 0],'FaceAlpha',0.18,'EdgeColor','none');
                end
            end
            title(sprintf('Raw1: %s (Fs=%.1f Hz)', guiL.file1, Fs));
            xlim([guiL.params.PlotStart guiL.params.PlotStart+guiL.params.PlotDuration]);
            hold off;
        else
            text(0.1,0.5,'Plot window outside raw1 range'); axis off;
        end
    else
        text(0.1,0.5,'No raw1 trace'); axis off;
    end
    % Raw2
    axes(h.axRaw2); cla;
    if ~isempty(guiL.raw2)
        Fs = guiL.Fs2; if isempty(Fs), Fs = 1; end
        t_raw = (0:numel(guiL.raw2)-1)/Fs;
        idx = t_raw >= guiL.params.PlotStart & t_raw <= (guiL.params.PlotStart + guiL.params.PlotDuration);
        if any(idx)
            plot(t_raw(idx), guiL.raw2(idx),'k-','LineWidth',0.6); hold on;
            sp = guiL.times2(guiL.times2>=guiL.params.PlotStart & guiL.times2<=guiL.params.PlotStart+guiL.params.PlotDuration);
            yl = ylim();
            for s = sp(:)'; plot([s s],[yl(1)+0.02*diff(yl) yl(1)+0.12*diff(yl)],'b'); end
            for b = guiL.bursts2
                a = max(b.start, guiL.params.PlotStart); bE = min(b.end, guiL.params.PlotStart+guiL.params.PlotDuration);
                if a < bE
                    patch([a bE bE a],[yl(1) yl(1) yl(2) yl(2)],[0.7 0 0],'FaceAlpha',0.16,'EdgeColor','none');
                end
            end
            title(sprintf('Raw2: %s (Fs=%.1f Hz)', guiL.file2, Fs));
            xlim([guiL.params.PlotStart guiL.params.PlotStart+guiL.params.PlotDuration]);
            hold off;
        else
            text(0.1,0.5,'Plot window outside raw2 range'); axis off;
        end
    else
        text(0.1,0.5,'No raw2 trace'); axis off;
    end
end

function draw_ISI(figHandle, guiL)
    h = guiL.h;
    axes(h.axISI); cla; hold on;
    if guiL.activeISI==1
        ISI = guiL.ISI1; ISI_t = guiL.ISI1_t; mask = guiL.spikeMask1;
        title(sprintf('Trace 1 ISIs (thresh=%.3f s)', guiL.params.BurstISI));
    else
        ISI = guiL.ISI2; ISI_t = guiL.ISI2_t; mask = guiL.spikeMask2;
        title(sprintf('Trace 2 ISIs (thresh=%.3f s)', guiL.params.BurstISI));
    end
    if isempty(ISI)
        text(0.1,0.5,'Not enough spikes to show ISIs'); axis off; return;
    end
    scatter(ISI_t, ISI, 18, [0.5 0.5 0.5], 'filled'); hold on;
    % Mark ISIs that belong to detected bursts (pairs: previous or next spike part of mask)
    pairMask = false(size(ISI));
    N = numel(mask);
    for k = 1:numel(ISI)
        s1 = (k <= N) && mask(k);
        s2 = ((k+1) <= N) && mask(k+1);
        pairMask(k) = s1 || s2;
    end
    scatter(ISI_t(pairMask), ISI(pairMask), 30, 'g', 'filled');
    plot([min(ISI_t) max(ISI_t)], [guiL.params.BurstISI guiL.params.BurstISI], 'r--','LineWidth',1.2);
    ylabel('ISI (s)'); xlabel('Time (s)');
    set(gca,'YScale','log'); grid on;
    xlim([min(ISI_t) max(ISI_t)]);
    hold off;
end

function draw_raster(figHandle, guiL)
    h = guiL.h;
    axes(h.axRaster); cla; hold on;
    yl1 = 2; yl2 = 1;
    plot(guiL.times1, yl1*ones(size(guiL.times1)),'k.','MarkerSize',6);
    plot(guiL.times2, yl2*ones(size(guiL.times2)),'k.','MarkerSize',6);
    for b = guiL.bursts1
        if b.end >= guiL.params.PlotStart && b.start <= guiL.params.PlotStart + guiL.params.PlotDuration
            patch([max(b.start,guiL.params.PlotStart) min(b.end,guiL.params.PlotStart+guiL.params.PlotDuration) ...
                   min(b.end,guiL.params.PlotStart+guiL.params.PlotDuration) max(b.start,guiL.params.PlotStart)],...
                  [yl1-0.2 yl1-0.2 yl1+0.2 yl1+0.2],[0 0.6 0],'FaceAlpha',0.25,'EdgeColor','none');
        end
    end
    for b = guiL.bursts2
        if b.end >= guiL.params.PlotStart && b.start <= guiL.params.PlotStart + guiL.params.PlotDuration
            patch([max(b.start,guiL.params.PlotStart) min(b.end,guiL.params.PlotStart+guiL.params.PlotDuration) ...
                   min(b.end,guiL.params.PlotStart+guiL.params.PlotDuration) max(b.start,guiL.params.PlotStart)],...
                  [yl2-0.2 yl2-0.2 yl2+0.2 yl2+0.2],[0.7 0 0],'FaceAlpha',0.18,'EdgeColor','none');
        end
    end
    % coincident onset lines
    if isfield(guiL,'coincSummary') && ~isempty(guiL.coincSummary) && isfield(guiL.coincSummary,'pairs')
        for k = 1:size(guiL.coincSummary.pairs,1)
            i = guiL.coincSummary.pairs(k,1); j = guiL.coincSummary.pairs(k,2);
            s1 = guiL.bursts1(i).start; s2 = guiL.bursts2(j).start;
            if (s1 >= guiL.params.PlotStart && s1 <= guiL.params.PlotStart+guiL.params.PlotDuration) || ...
               (s2 >= guiL.params.PlotStart && s2 <= guiL.params.PlotStart+guiL.params.PlotDuration)
                plot([s1 s2],[2 1],'r-','LineWidth',1.2);
            end
        end
    end
    ylim([0.5 2.5]); yticks([1 2]); yticklabels({'Trace 2','Trace 1'});
    xlim([guiL.params.PlotStart guiL.params.PlotStart + guiL.params.PlotDuration]);
    xlabel('Time (s)'); title('Raster and detected bursts');
    hold off;
end

function update_stats_list(figHandle, guiL)
    h = guiL.h;
    lines = {};
    lines{end+1} = sprintf('File1: %s', guiL.file1);
    lines{end+1} = sprintf(' spikes: %d   bursts: %d', numel(guiL.times1), numel(guiL.bursts1));
    if ~isempty([guiL.bursts1.duration])
        lines{end+1} = sprintf(' mean dur: %.3f s   mean spikes/burst: %.2f', nanmean([guiL.bursts1.duration]), nanmean([guiL.bursts1.nSpikes]));
    else
        lines{end+1} = ' mean dur: N/A';
    end
    lines{end+1} = ' ';
    lines{end+1} = sprintf('File2: %s', guiL.file2);
    lines{end+1} = sprintf(' spikes: %d   bursts: %d', numel(guiL.times2), numel(guiL.bursts2));
    if ~isempty([guiL.bursts2.duration])
        lines{end+1} = sprintf(' mean dur: %.3f s   mean spikes/burst: %.2f', nanmean([guiL.bursts2.duration]), nanmean([guiL.bursts2.nSpikes]));
    else
        lines{end+1} = ' mean dur: N/A';
    end
    lines{end+1} = ' ';
    if isfield(guiL,'coincSummary')
        lines{end+1} = sprintf('Coinc pairs: %d', guiL.coincSummary.nPairs);
        lines{end+1} = sprintf(' pct file1 coinc: %.1f%%', guiL.coincSummary.pctBursts1Coincident);
        lines{end+1} = sprintf(' pct file2 coinc: %.1f%%', guiL.coincSummary.pctBursts2Coincident);
    else
        lines{end+1} = 'Coinc pairs: 0';
    end
    lines{end+1} = ' ';
    lines{end+1} = sprintf('Params: BurstISI=%.3f  MinSpikes=%d  CoincWin=%.3f', guiL.params.BurstISI, guiL.params.MinSpikes, guiL.params.CoincWin);
    lines{end+1} = sprintf('Plot window: %.2f - %.2f s', guiL.params.PlotStart, guiL.params.PlotStart+guiL.params.PlotDuration);
    set(h.lstStats,'String',lines);
end

function [bursts, spikeMask] = detect_bursts_threshold(times, thresh, minSpikes)
    bursts = struct('start',[],'end',[],'duration',[],'nSpikes',[],'spikeTimes',{});
    spikeMask = false(size(times));
    if isempty(times) || numel(times) < 2, return; end
    isis = diff(times);
    short = isis < thresh;
    if ~any(short), return; end
    dshort = diff([0; short(:); 0]);
    runStarts = find(dshort==1);
    runEnds = find(dshort==-1)-1;
    cnt = 0;
    for r = 1:numel(runStarts)
        s = runStarts(r); e = runEnds(r);
        Lshort = e - s + 1;
        nSpikes = Lshort + 1;
        if nSpikes >= minSpikes
            spikeIdxs = s:(e+1);
            cnt = cnt + 1;
            bursts(cnt).start = times(spikeIdxs(1));
            bursts(cnt).end = times(spikeIdxs(end));
            bursts(cnt).duration = bursts(cnt).end - bursts(cnt).start;
            bursts(cnt).nSpikes = numel(spikeIdxs);
            bursts(cnt).spikeTimes = times(spikeIdxs);
            spikeMask(spikeIdxs) = true;
        end
    end
    if cnt==0
        bursts = struct('start',[],'end',[],'duration',[],'nSpikes',[],'spikeTimes',{});
        spikeMask = false(size(times));
    end
end

function coincSummary = match_coincidence(b1, b2, coincWin)
    coincSummary.pairs = [];
    if isempty(b1) || isempty(b2)
        coincSummary.nPairs = 0;
        coincSummary.pctBursts1Coincident = 0;
        coincSummary.pctBursts2Coincident = 0;
        return;
    end
    pairs = [];
    for i=1:numel(b1)
        for j=1:numel(b2)
            s1 = b1(i).start; e1 = b1(i).end;
            s2 = b2(j).start; e2 = b2(j).end;
            onsetDiff = abs(s1 - s2);
            overlap = (s1 <= e2) && (s2 <= e1);
            if (onsetDiff <= coincWin) || overlap
                pairs(end+1,:) = [i j]; %#ok<AGROW>
            end
        end
    end
    if isempty(pairs)
        coincSummary.pairs = [];
        coincSummary.nPairs = 0;
        coincSummary.pctBursts1Coincident = 0;
        coincSummary.pctBursts2Coincident = 0;
    else
        coincSummary.pairs = pairs;
        coincSummary.nPairs = size(pairs,1);
        coincSummary.pctBursts1Coincident = 100 * numel(unique(pairs(:,1))) / max(1,numel(b1));
        coincSummary.pctBursts2Coincident = 100 * numel(unique(pairs(:,2))) / max(1,numel(b2));
    end
end

function stats = compute_burst_stats(bursts)
    if isempty(bursts)
        stats.spikesPerBurst = [];
        stats.durations = [];
        stats.periods = [];
        stats.duty = [];
        return;
    end
    spikesPerBurst = [bursts.nSpikes];
    durations = [bursts.duration];
    starts = [bursts.start];
    periods = diff(starts);
    if isempty(periods)
        duty = NaN(size(durations));
    else
        duty_temp = durations(1:end-1) ./ periods;
        duty = NaN(size(durations));
        duty(1:numel(duty_temp)) = duty_temp;
    end
    stats.spikesPerBurst = spikesPerBurst;
    stats.durations = durations;
    stats.periods = periods;
    stats.duty = duty;
end

function [ISI, tmid, idx] = compute_isi_values(times)
    ISI = []; tmid = []; idx = [];
    if isempty(times) || numel(times) < 2, return; end
    times = times(:);
    ISI = diff(times);
    tmid = (times(1:end-1) + times(2:end))/2;
    idx = (1:numel(ISI))';
end

function [S, times, spikeidx] = load_struct_and_spikes(fname, preferredSname)
    tmp = load(fname);
    S = []; times = []; spikeidx = [];
    fn = fieldnames(tmp);
    chosen = '';
    if any(strcmp(fn,preferredSname)); chosen = preferredSname;
    else
        for k=1:numel(fn)
            if ischar(fn{k}) && startsWith(fn{k},'S')
                chosen = fn{k}; break;
            end
        end
    end
    if ~isempty(chosen)
        S = tmp.(chosen);
        if isfield(S,'spiketimesec'), times = double(S.spiketimesec(:)); end
        if isfield(S,'spikeidx'), spikeidx = double(S.spikeidx(:)); end
    else
        if isfield(tmp,'spiketimesec'), times = double(tmp.spiketimesec(:)); end
        if isfield(tmp,'spikeidx'), spikeidx = double(tmp.spikeidx(:)); end
    end
    if isempty(times) && isfield(tmp,'spiketimes'), times = double(tmp.spiketimes(:)); end
end

function [raw, Fs, nameUsed] = try_load_abf_or_mat_raw(matfile, explicitABF, chan, useABFflag, Sstruct)
    raw = []; Fs = []; nameUsed = '(none)';
    if useABFflag
        candidates = {};
        if ~isempty(explicitABF), candidates{end+1} = explicitABF; end
        [pdir, base, ~] = fileparts(matfile);
        candidates{end+1} = fullfile(pdir, [base '.abf']);
        listing = dir(fullfile(pdir,'*.abf'));
        for k=1:numel(listing)
            candidates{end+1} = fullfile(pdir, listing(k).name);
        end
        candidates = unique(candidates,'stable');
        for c = 1:numel(candidates)
            f = candidates{c};
            if ~exist(f,'file'), continue; end
            try
                [data, si, h] = abfload(f); %#ok<NASGU,ASGLU>
                if ndims(data) == 3
                    if size(data,2) >= chan, data_sel = squeeze(data(:,chan,1));
                    else data_sel = squeeze(data(:,1,1)); end
                elseif ndims(data) == 2
                    if size(data,2) >= chan, data_sel = data(:,chan); else data_sel = data(:,1); end
                else
                    data_sel = data(:);
                end
                raw = double(data_sel(:));
                nameUsed = f;
                if exist('si','var') && ~isempty(si) && isnumeric(si) && si>0
                    Fs = 1e6/si;
                elseif exist('h','var') && isstruct(h) && isfield(h,'SI')
                    Fs = 1e6/h.SI;
                end
                return;
            catch %#ok<CTCH>
                continue;
            end
        end
    end
    % fallback: try to get a large numeric vector from Sstruct fields (if present)
    if ~isempty(Sstruct) && isstruct(Sstruct)
        fn = fieldnames(Sstruct);
        for k=1:numel(fn)
            val = Sstruct.(fn{k});
            if isnumeric(val) && isvector(val) && numel(val)>1000
                raw = double(val(:)); nameUsed = sprintf('from .mat (%s)', fn{k}); break;
            end
        end
    end
end