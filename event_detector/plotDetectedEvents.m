function plotDetectedEvents(output, path_plot)
% Plots a summary of the results returned by detectEvents. This includes
% the power spectrum, sleep scoring, length of data, detected slow waves
% and spindles, their phase coupling etc., depending on what is present in
% the data.
%
% Note: This was mainly tested using 4 channels. Placing several plots on
% top of each other is tricky, so every help is needed to make this more
% robust for other numnbers of channels.
%
% INPUT VARIABLES:
% output						Output of detectEvents
% path_plot						String; full path for saving the plot
%								(e.g., 'C:\onepath\file.png')
%
% OUTPUT VARIABLES:
% -
%
% AUTHOR:
% Jens Klinzing, klinzing@princeton.edu

%% INITIALIZATION AND SETTINGS
Fs							= output.info.Fs;			% sample rate of the artifact definition
if isfield(output, 'slo')
    slo						= 1;
    det_slo					= output.slo.neg_peaks;
else
    slo						= 0;
end
if isfield(output, 'spi')
    spi						= 1;
    det_spi					= output.spi.events;
    % Find center of each spindle
    for iCh = 1:numel(det_spi)
        det_spi{iCh} = det_spi{iCh}(1,:) + (det_spi{iCh}(2,:) - det_spi{iCh}(1,:))/2;
    end
else
    spi						= 0;
end
if isfield(output, 'rip')
    rip						= 1;
    det_rip					= output.rip.events;
    % Find center of each ripple
    for iCh = 1:numel(det_rip)
        if ~isempty(det_rip{iCh})
            det_rip{iCh} = det_rip{iCh}(1,:) + (det_rip{iCh}(2,:) - det_rip{iCh}(1,:))/2;
        end
    end
else
    rip						= 0;
end
if isfield(output, 'spectrum')
    spectrum				= 1;
else
    spectrum				= 0;
end
chans						= output.info.channel;
num_chans					= numel(chans);
if isfield(output, 'slo') && isfield(output, 'spi') && (numel(det_spi) ~= numel(det_slo) || numel(det_spi) ~= num_chans)
    error('Number of channels is not the same in all event types or doesnt match channel number promised by event detection structure.')
end
if num_chans > 6, warning('Not sure if that many channels will plot well.'), end

% Load hypnogram
hypnogram                   = output.info.scoring;
epochLength                 = output.info.scoring_epoch_length;

%% PREPARATIONS
% Translate hypnogram to seconds
% TODO: Why are we doing this? Legacy reason?
hypnogram_plot = [];
hypnogram_plot.Wake = [];
for iStNREM = 1:length(output.info.cfg.code_NREM) %create based on number of NREM stages
    hypnogram_plot.(strcat('S',num2str(iStNREM))) = [];
end
hypnogram_plot.REM = [];

for iSt = 1:length(hypnogram)
    sec = (iSt-1)*epochLength+1;    % translate the current position into seconds
    if hypnogram(iSt) == output.info.cfg.code_WAKE(1)
        for iEp = sec:sec+epochLength-1
            hypnogram_plot.Wake(end+1)      = iEp;
        end
    end
    if hypnogram(iSt) == output.info.cfg.code_NREM(1)
        for iEp = sec:sec+epochLength-1
            hypnogram_plot.S1(end+1)    = iEp;
        end
    end
    if length(output.info.cfg.code_NREM)>=2
        if hypnogram(iSt) == output.info.cfg.code_NREM(2)
            for iEp = sec:sec+epochLength-1
                hypnogram_plot.S2(end+1)    = iEp;
            end
        end
        if hypnogram(iSt) == output.info.cfg.code_NREM(3)
            for iEp = sec:sec+epochLength-1
                hypnogram_plot.S3(end+1)    = iEp;
            end
        end
        if hypnogram(iSt) == output.info.cfg.code_NREM(4)
            for iEp = sec:sec+epochLength-1
                hypnogram_plot.S4(end+1)    = iEp;
            end
        end
    end
    if hypnogram(iSt) == output.info.cfg.code_REM(1)
        for iEp = sec:sec+epochLength-1
            hypnogram_plot.REM(end+1)   = iEp;
        end
    end
end

% Delete all empty fields/sleep stages from array
sleepStage               = fieldnames(hypnogram_plot);
for iSt=1:numel(sleepStage)
    if isempty(hypnogram_plot.(sleepStage{iSt}))
        hypnogram_plot = rmfield(hypnogram_plot,sleepStage{iSt});
    end
end
sleepStage               = fieldnames(hypnogram_plot);

% Some settings
lineHeight_hyp           = 0;                % draw it a bit under the x axis
cmap.Wake                = [.85 .16 .0];        % some nice colors for the lines
cmap.S1                  = [1 .6 .2];
cmap.S2                  = [.7 .9 1];
cmap.S3                  = [.4 .6 1];
cmap.S4                  = [.1 .1 1];
cmap.REM                 = [0 .5 0];

%% PREPARE PLOTTING
num_col		= 6; % flexibly usable columns
num_pnl		= 2; % event and spectrum panel

%% PLOT EVENT PANEL
pnl			= 1; % first panel
splts		= 1:num_chans*pnl*num_col;
splts_r		= sort([num_col*[0:num_chans-1]+num_col-1 num_col*[0:num_chans-1]+num_col]);
splts_l		= splts(~ismember(splts, splts_r)); % and right side

h = figure;
% set(h,'Position',[20 100 1400 num_chans*130+500]);
set(h, 'Position', get(0, 'Screensize')); % fullscreen

% Left: Scoring and event distribution
subplot(num_chans*num_pnl,num_col,splts_l)
if isfield(output.info, 'name') && ~isempty(output.info.name)
    title(['Dataset: ' output.info.name ', recording length: ' num2str(output.info.length / output.info.Fs/60) ' min'])
else
    title(['Recording length: ' num2str(output.info.length / output.info.Fs/60) ' min'])
end

hold on
sleepStage_legend        = cell(0);  % for the legend, collect names of sleep stages that have actually occured
handlevector = [];
for iSt = sleepStage' % Needs ' for whatever reason
    t                       = hypnogram_plot.(iSt{1});
    t_diffs                 = diff(t) == 1; % find consecutive occurrences of this stage (= episodes)
    t_diffs                 = [0 t_diffs 0]; % so we also find the first and last one
    ep_starts               = strfind(t_diffs,[0 1]); % where does scoring flip to this stage
    ep_ends                 = strfind(t_diffs,[1 0]);
    for iEp = 1:length(ep_starts)
        temphandle			= line([t(ep_starts(iEp))/60 t(ep_ends(iEp))/60],[lineHeight_hyp lineHeight_hyp], 'LineWidth', 30, 'Color', cmap.(iSt{1}));
    end
    handlevector			= [handlevector temphandle];
    sleepStage_legend		= [sleepStage_legend iSt{1}]; % add handle of one of the lines
end

temphandle = line([0 output.info.length/Fs/60],[lineHeight_hyp+0.05 lineHeight_hyp+0.05], 'LineWidth', 1, 'Color', 'black');
handlevector = [temphandle handlevector];

% Plot slow oscillations and spindles
if slo
    for iCh = 1:num_chans
        temphandle      = plot(det_slo{iCh}/Fs/60, lineHeight_hyp - 0.1*iCh, 'ko'); hold on
        if iCh == 1
            handlevector(end+1) = temphandle(1); % each event is a separate plot, only want one legend entry
        end
        clear temphandle
    end
end
if spi
    for iCh = 1:num_chans
        temphandle      = plot(det_spi{iCh}/Fs/60, lineHeight_hyp - (0.1*iCh)-0.02, 'ro'); hold on
        if iCh == 1
            handlevector(end+1) = temphandle(1); % each event is a separate plot, only want one legend entry
        end
        clear temphandle
    end
end
if rip
    for iCh = 1:num_chans
        if ~isempty(det_rip{iCh})
            temphandle      = plot(det_rip{iCh}/Fs/60, lineHeight_hyp - (0.1*iCh)-0.04, 'go'); hold on
            if iCh == 1
                handlevector(end+1) = temphandle(1); % each event is a separate plot, only want one legend entry
            end
            clear temphandle
        end
    end
end
% Labels, title, legend
ylim([lineHeight_hyp - (0.1*num_chans)-0.1 .2])
xlabel('Time (in min)')
xlim([0 inf])
set(gca,'YTick', [])
set(gca,'YTickLabel', [])
ylabel('Channels (in the order provided)');
if slo && spi && rip
    leg_str = {'SOs', 'Spindles', 'Ripples'};
elseif slo && spi
    leg_str = {'SOs', 'Spindles'};
elseif slo
    leg_str = {'SOs'};
elseif spi
    leg_str = {'Spindles'};
elseif rip
    leg_str = {'Ripples'};
else
    leg_str = {};
end
legend(handlevector(:), 'Rec length', sleepStage_legend{:}, leg_str{:}, 'Location', 'SouthWest')

if slo
    % Average SOs
    color_line				= [0, .447, .741];
    color_area				= [0, .447, .741];
    alpha					= .4;
    line_width				= 1.5;
    err_metric				= 'std';
    win_width				= floor(size(output.slo.waveform{iCh}(1,:),2)/2); % in samples
    for iCh = 1:num_chans
        p = subplot(num_chans*num_pnl,num_col,splts_r([(iCh-1)*2+1 (iCh-1)*2+2]));
        y_mean = mean(output.slo.waveform{iCh},1);
        y_std = std(output.slo.waveform{iCh},0,1);
        switch(err_metric)
            case 'std', err = y_std;
            case 'sem', err = (y_std./sqrt(size(output.slo.waveform{iCh},1)));
            case 'var', err = (y_std.^2);
            case 'c95', err = (y_std./sqrt(size(output.slo.waveform{iCh},1))).*1.96;
        end
        x_axis = (-win_width:win_width)/Fs; % in seconds
        x_vector = [x_axis, fliplr(x_axis)];
        yyaxis right
        if size(output.slo.waveform{iCh}, 2) > 1
            patch = fill(x_vector, [y_mean+err,fliplr(y_mean-err)], color_area);
            set(patch, 'edgecolor', 'none');
            set(patch, 'FaceAlpha', alpha);
            hold on;
        end
        plot(x_axis, y_mean, 'color', color_line, 'LineWidth', line_width, 'LineStyle', '-', 'Color', 'black');
        xlim([x_axis(1) x_axis(end)])
        ylabel(output.info.channel{iCh}, 'Interpreter', 'none', 'Color', 'black');
        ax = gca;
        ax.YAxis(1).Visible = 0; % hide left unused axis
        ax.YAxis(2).Color = 'k'; % turn right axis black
        y_lim = ylim;
        spi_y = y_lim(1) + (y_lim(2)-y_lim(1))/3; % plot spindles at 1/3 of the ylim height
        
        if spi
            % Draw spindles into SO plot
            for iEv = 1:size(output.spi.events{iCh},2)
                spi_center = round(output.spi.events{iCh}(1,iEv) + (output.spi.events{iCh}(2,iEv)-output.spi.events{iCh}(1,iEv))/2);
                spi_cooc = find(spi_center > output.slo.neg_peaks{iCh}-win_width & spi_center < output.slo.neg_peaks{iCh}+win_width);
                % in case there is at least one co-occurrence, mark the one closest to the downpeak in SO plot
                if ~isempty(spi_cooc)
                    spi_time = (spi_center - output.slo.neg_peaks{iCh}(spi_cooc)) / Fs;
                    spi_time = spi_time(nearest(spi_time, 0)); % choose only the one closest to 0
                    plot(spi_time, spi_y, 'ko', 'MarkerSize', 5, 'MarkerEdgeColor', 'red') % -50 is in microvolt, this could be chosen better
                end
            end
            
            % Add polar histograms
            if isfield(output, 'SloSpiAmpCoupling')
                axes('pos',[p.Position(1) p.Position(2)+p.Position(4)*.11 p.Position(3)*.2 p.Position(3)*.2]) % left, bottom, width, height
                pol = polarhistogram(deg2rad(output.SloSpiAmpCoupling{iCh}));
                pol.Parent.ThetaAxis.Visible = 'off';
                pol.Parent.RAxis.Visible = 'off';
                if iCh == 1
                    title(sprintf(['SO Phase spindle band amp peak\n n=' num2str(numel(output.SloSpiAmpCoupling{iCh}))]),  'FontSize', 6)
                else
                    title(['n=' num2str(numel(output.SloSpiAmpCoupling{iCh}))],  'FontSize', 6)
                end
            end
            if isfield(output, 'SloSpiDetCoupling')
                axes('pos',[p.Position(1)+p.Position(3)*.8 p.Position(2)+p.Position(4)*.11 p.Position(3)*.2 p.Position(3)*.2]) % left (first argument) had +p.Position(3)*.7
                pol = polarhistogram(deg2rad(output.SloSpiDetCoupling{iCh}));
                pol.Parent.ThetaAxis.Visible = 'off';
                pol.Parent.RAxis.Visible = 'off';
                if iCh == 1
                    title(sprintf(['SO Phase detected spindle peak\n n=' num2str(numel(output.SloSpiDetCoupling{iCh}))]),  'FontSize', 6, 'Color', 'red')
                else
                    title(['n=' num2str(numel(output.SloSpiDetCoupling{iCh}))],  'FontSize', 6, 'Color', 'red')
                end
            end
            %plot spindle ripple coupling
            if isfield(output, 'SpiRipDetCoupling')
                axes('pos',[p.Position(1)+p.Position(3)*.8 p.Position(2)+p.Position(4)*.71 p.Position(3)*.2 p.Position(3)*.2]) % left (first argument) had +p.Position(3)*.7
                pol = polarhistogram(deg2rad(output.SpiRipDetCoupling{iCh}));
                pol.Parent.ThetaAxis.Visible = 'off';
                pol.Parent.RAxis.Visible = 'off';
                pol.FaceColor ='g';
                if iCh == 1
                    title(sprintf(['Spi Phase detected ripple peak\n n=' num2str(numel(output.SloSpiDetCoupling{iCh}))]),  'FontSize', 6, 'Color', 'green')
                else
                    title(['n=' num2str(numel(output.SpiRipDetCoupling{iCh}))],  'FontSize', 6, 'Color', 'green')
                end
            end
        end
        
    end
end
drawnow

%% PLOT SPECTRUM PANEL
if spectrum
    pnl			= 2; % first panel
    splts		= (pnl-1)*(num_chans*num_col)+1:num_chans*num_col*pnl; % this subplot business is getting confusing, are there better ways to flexibly manage subplots?
    splts_l = splts(logical(mod(ceil((1:length(splts)) / (num_col/2)), 2)));
    splts_r = splts(~logical(mod(ceil((1:length(splts)) / (num_col/2)), 2)));
    
    % Left: NREM
    if isfield(output.spectrum, 'rel_nrem')
        s1 = subplot(num_chans*num_pnl,num_col,splts_l);
        plot(output.spectrum.freq, output.spectrum.rel_nrem, '-', 'LineWidth', 2)
        ylim([-1 max(max(output.spectrum.rel_nrem(:, output.spectrum.freq > 2 & output.spectrum.freq < 30))) * 2])
        ylabel('Oscillatory relative to fractal component')
        title(sprintf('\n Spectrum NREM sleep'))
        legend(output.info.channel{:})
        xlim([1 35]) % thats where the really interesting stuff happens
    end
    
    % Right: REM
    if isfield(output.spectrum, 'rel_rem')
        s2 = subplot(num_chans*pnl,num_col,splts_r);
        plot(output.spectrum.freq, output.spectrum.rel_rem, '-', 'LineWidth', 2)
        ylim([-1 max(max(output.spectrum.rel_rem(:, output.spectrum.freq > 2 & output.spectrum.freq < 30))) * 2])
        ylabel('Oscillatory relative to fractal component')
        legend(output.info.channel{:})
        title(sprintf('\n Spectrum REM sleep'))
    end
end

% Shift both subplots a but down to make everything look better
s1.Position = [s1.Position(1) s1.Position(2)-.04 s1.Position(3) s1.Position(4)];
s2.Position = [s2.Position(1) s2.Position(2)-.04 s1.Position(3) s2.Position(4)];

%% Save the plot
if nargin > 1 && ~isempty(path_plot)
    [~,name,ext] = fileparts(path_plot);
    saveas(gcf,path_plot)
    disp(['Plot has been saved at provided location (''' name ext ''').']);
    close all
end

