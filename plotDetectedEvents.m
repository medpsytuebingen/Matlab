function plotDetectedEvents(output, path_plot)
% Plots a summary (scoring, length of data, detected slow waves and
% spindles) of the results returned by detectEvents.
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

%% TODO
% . currently only plots SOs and spindles

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
hypnogram_plot.S1 = [];
hypnogram_plot.S2 = [];
hypnogram_plot.S3 = [];
hypnogram_plot.S4 = [];
hypnogram_plot.REM = [];

for iSt = 1:length(hypnogram)
	sec = (iSt-1)*epochLength+1;    % translate the current position into seconds
	if hypnogram(iSt) == 0
		for iEp = sec:sec+epochLength-1
			hypnogram_plot.Wake(end+1)      = iEp;
		end
	end
	if hypnogram(iSt) == 1
		for iEp = sec:sec+epochLength-1
			hypnogram_plot.S1(end+1)    = iEp;
		end
	end
	if hypnogram(iSt) == 2
		for iEp = sec:sec+epochLength-1
			hypnogram_plot.S2(end+1)    = iEp;
		end
	end
	if hypnogram(iSt) == 3
		for iEp = sec:sec+epochLength-1
			hypnogram_plot.S3(end+1)    = iEp;
		end
	end
	if hypnogram(iSt) == 4
		for iEp = sec:sec+epochLength-1
			hypnogram_plot.S4(end+1)    = iEp;
		end
	end
	if hypnogram(iSt) == 5
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

%% PLOT EVENT PANELS
pnl			= 1; % first panel
splts		= 1:num_chans*num_col*pnl;
splts_l		= sort([num_col*[0:num_chans-1]+1 num_col*[0:num_chans-1]+2]); % all subplot indices for the left side
splts_r		= splts(~ismember(splts, splts_l)); % and right side

h = figure;
set(h,'Position',[50 100 1400 num_chans*140]);
subplot(num_chans,num_col,splts_r)
hold on
sleepStage_legend        = cell(0);  % for the legend, collect names of sleep stages that have actually occured
handlevector = [];
for iSt = sleepStage' % Needs ' for whatever reason
	t					= hypnogram_plot.(iSt{1});
	t_diffs				= diff(t) == 1; % find consecutive occurrences of this stage (= episodes)
	t_diffs				= [0 t_diffs 0]; % so we also find the first and last one
	ep_starts			= strfind(t_diffs,[0 1]); % where does scoring flip to this stage
	ep_ends				= strfind(t_diffs,[1 0]);
	for iEp = 1:length(ep_starts)
		temphandle			= line([t(ep_starts(iEp))/60 t(ep_ends(iEp))/60],[lineHeight_hyp lineHeight_hyp], 'LineWidth', 20, 'Color', cmap.(iSt{1}));
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

% Labels, title, legend
ylim([lineHeight_hyp - (0.1*num_chans)-0.05 .1])
xlabel('Time (in min)')
xlim([0 inf])
set(gca,'YTick', [])
ylabel('Channels (in the order provided)');
sgtitle(sprintf('Quality check'),  'FontSize', 14);
if slo && spi
	leg_str = {'Slow oscillations/waves', 'Spindles'};
elseif slo
	leg_str = {'Slow oscillations/waves'};
elseif spi
	leg_str = {'Spindles'};
else
	leg_str = {};
end
legend(handlevector(:), 'Recording Length', sleepStage_legend{:}, leg_str{:}, 'Location', 'EastOutside')

if slo
	% Average SOs
	color_line				= [0, .447, .741];
	color_area				= [0, .447, .741];
	alpha					= .4;
	line_width				= 1.5;
	err_metric				= 'std';
	win_width				= floor(size(output.slo.waveform{iCh}(1,:),2)/2); % in samples
	for iCh = 1:num_chans
		p = subplot(num_chans,num_col,splts_l([(iCh-1)*2+1 (iCh-1)*2+2]));
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
		if size(output.slo.waveform{iCh}, 2) > 1
			patch = fill(x_vector, [y_mean+err,fliplr(y_mean-err)], color_area);
			set(patch, 'edgecolor', 'none');
			set(patch, 'FaceAlpha', alpha);
			hold on;
		end
		plot(x_axis, y_mean, 'color', color_line, 'LineWidth', line_width, 'LineStyle', '-', 'Color', 'black');
		xlim([x_axis(1) x_axis(end)])
		ylabel(output.info.channel{iCh}, 'Interpreter', 'none');
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
				axes('pos',[p.Position(1)-p.Position(3)*.1 p.Position(2)+p.Position(4)*.2 p.Position(3)*.4 p.Position(3)*.4]) % left, bottom, width, height
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
				axes('pos',[p.Position(1)+p.Position(3)*.7 p.Position(2)+p.Position(4)*.2 p.Position(3)*.4 p.Position(3)*.4]) % left, bottom, width, height
				pol = polarhistogram(deg2rad(output.SloSpiDetCoupling{iCh}));
				pol.Parent.ThetaAxis.Visible = 'off';
				pol.Parent.RAxis.Visible = 'off';
				if iCh == 1
					title(sprintf(['SO Phase detected spindle peak\n n=' num2str(numel(output.SloSpiDetCoupling{iCh}))]),  'FontSize', 6, 'Color', 'red')
				else
					title(['n=' num2str(numel(output.SloSpiDetCoupling{iCh}))],  'FontSize', 6, 'Color', 'red')
				end
			end
		end
	end
end
drawnow

% Save the plot
if nargin > 1 && ~isempty(path_plot)
	[~,name,ext] = fileparts(path_plot);
	saveas(gcf,path_plot)
	disp(['Plot has been saved at provided location (''' name ext ''').']);
	close all
end

