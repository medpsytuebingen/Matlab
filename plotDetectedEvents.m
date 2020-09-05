function plotDetectedEvents(detectedEvents, path_plot)
% Plots a summary (scoring, length of data, detected slow waves and
% spindles) of the results returned by detectEvents.
%
% INPUT VARIABLES:
% detectedEvents				Output of detectEvents
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
% . allow providing an output file for saving result

%% INITIALIZATION AND SETTINGS
Fs							= detectedEvents.info.Fs;			% sample rate of the artifact definition
det_slo						= detectedEvents.slo.neg_peaks;
det_spi						= detectedEvents.spi.events;
chans						= detectedEvents.info.channel;
num_chans					= numel(chans);
if numel(det_spi) ~= numel(det_slo) || numel(det_spi) ~= num_chans
	error('Number of channels is not the same in all event types or doesnt match channel number promised by event detection structure.')
end
if num_chans > 6, warning('Not sure if that many channels will plot well.'), end

% Load hypnogram
hypnogram                   = detectedEvents.info.scoring;
epochLength                 = detectedEvents.info.scoring_epoch_length;

% Find center of each spindle
for iCh = 1:numel(det_spi)
	det_spi{iCh} = det_spi{iCh}(1,:) + (det_spi{iCh}(2,:) - det_spi{iCh}(1,:))/2;
end

%% PREPARATIONS

% Translate hypnogram to seconds
% Why are we doing this?
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
% - TODO: Can this be replaced by an isempty later on?
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

%% START
% TODO: Spindle and SO legend entries might be switched
% Calculate subplots to use for each of the two panels
splts = 1:num_chans * 3;
splts_r = splts(~any(splts == (3*[0:num_chans-1]+1)',1));
splts_l = splts(any(splts == (3*[0:num_chans-1]+1)',1));

h = figure;
set(h,'Position',[50 100 1400 num_chans*120]);
subplot(num_chans,3,splts_r)
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

temphandle = line([0 detectedEvents.info.length/Fs/60],[lineHeight_hyp+0.05 lineHeight_hyp+0.05], 'LineWidth', 1, 'Color', 'black');
handlevector = [temphandle handlevector];

% Plot slow oscillations and spindles
for iCh = 1:num_chans
	temphandle      = plot(det_slo{iCh}/Fs/60, lineHeight_hyp - 0.1*iCh, 'ko'); hold on
	if iCh == 1
		handlevector(end+1) = temphandle(1); % each event is a separate plot, only want one legend entry
	end
	clear temphandle
end
for iCh = 1:num_chans
	temphandle      = plot(det_spi{iCh}/Fs/60, lineHeight_hyp - (0.1*iCh)-0.02, 'ro'); hold on
	if iCh == 1
		handlevector(end+1) = temphandle(1); % each event is a separate plot, only want one legend entry
	end
	clear temphandle
end

% Plot detected NREM and REM episodes % change output to detectedEvents
% for iEp = 1:size(output.NREMepisode,2)
% 	temphandle			= line([output.NREMepisode(1,iEp) output.NREMepisode(2,iEp)],[lineHeight_hyp+.1 lineHeight_hyp+.1], 'LineWidth', 20, 'Color', 'black');
% end
% for iEp = 1:size(output.REMepisode,2)
% 	temphandle			= line([output.REMepisode(1,iEp) output.REMepisode(2,iEp)],[lineHeight_hyp+.2 lineHeight_hyp+.2], 'LineWidth', 20, 'Color', [.4 .4 .4]);
% end

% --- Labels, title, legend
ylim([lineHeight_hyp - (0.1*iCh)-0.05 .1])
xlabel('Time (in min)')
xlim([0 inf])
set(gca,'YTick', [])
ylabel('Channels (in the order provided)');
sgtitle(sprintf('Quality check'),  'FontSize', 14);
legend(handlevector(:), 'Recording Length', sleepStage_legend{:}, 'Slow oscillations', 'Spindles', 'Location', 'EastOutside')

% Average SOs
color_line				= [0, .447, .741];
color_area				= [0, .447, .741];
alpha					= .4;
line_width				= 1.5;
err_metric				= 'std';
win_width				= floor(size(detectedEvents.slo.waveform{iCh}(1,:),2)/2); % in samples
for iCh = 1:num_chans
	subplot(num_chans,3,splts_l(iCh));
	y_mean = mean(detectedEvents.slo.waveform{iCh},1);
	y_std = std(detectedEvents.slo.waveform{iCh},0,1);
	switch(err_metric)
		case 'std', err = y_std;
		case 'sem', err = (y_std./sqrt(size(detectedEvents.slo.waveform{iCh},1)));
		case 'var', err = (y_std.^2);
		case 'c95', err = (y_std./sqrt(size(detectedEvents.slo.waveform{iCh},1))).*1.96;
	end
	x_axis = (-win_width:win_width)/Fs; % in seconds
	x_vector = [x_axis, fliplr(x_axis)];
	if size(detectedEvents.slo.waveform{iCh}, 2) > 1
		patch = fill(x_vector, [y_mean+err,fliplr(y_mean-err)], color_area);
		set(patch, 'edgecolor', 'none');
		set(patch, 'FaceAlpha', alpha);
		hold on;
	end
	plot(x_axis, y_mean, 'color', color_line, 'LineWidth', line_width, 'LineStyle', '-', 'Color', 'black');
	xlim([x_axis(1) x_axis(end)])
	ylabel(detectedEvents.info.channel{iCh}, 'Interpreter', 'none');
	y_lim = ylim;
	spi_y = y_lim(1) + (y_lim(2)-y_lim(1))/3; % plot spindles at 1/3 of the ylim height
	
	% Draw spindles into SO plot
	for iEv = 1:size(detectedEvents.spi.events{iCh},2)
		spi_center = round(detectedEvents.spi.events{iCh}(1,iEv) + (detectedEvents.spi.events{iCh}(2,iEv)-detectedEvents.spi.events{iCh}(1,iEv))/2);
		spi_cooc = find(spi_center > detectedEvents.slo.neg_peaks{iCh}-win_width & spi_center < detectedEvents.slo.neg_peaks{iCh}+win_width);
		% in case there is at least one co-occurrence, mark the one closest to the downpeak in SO plot
		if ~isempty(spi_cooc)
			spi_time = (spi_center - detectedEvents.slo.neg_peaks{iCh}(spi_cooc)) / Fs;
			spi_time = spi_time(nearest(spi_time, 0)); % choose only the one closest to 0
			plot(spi_time, spi_y, 'ko', 'MarkerSize', 5, 'MarkerEdgeColor', 'red') % -50 is in microvolt, this could be chosen better
		end
	end
end
drawnow
% Add polar histograms (requires a second loop because Matlab...)
for iCh = 1:num_chans
	p = subplot(num_chans,3,splts_l(iCh));
	axes('pos',[p.Position(1)-p.Position(3)*.1 p.Position(2)+p.Position(4)*.2 p.Position(3)*.4 p.Position(3)*.4]) % left, bottom, width, height
	pol = polarhistogram(deg2rad(detectedEvents.SloSpiAmpCoupling{iCh}));
	pol.Parent.ThetaAxis.Visible = 'off';
	pol.Parent.RAxis.Visible = 'off';
	if iCh == 1
		title(sprintf(['SO Phase spindle band amp peak\n n=' num2str(numel(detectedEvents.SloSpiAmpCoupling{iCh}))]),  'FontSize', 6)
	else
		title(['n=' num2str(numel(detectedEvents.SloSpiAmpCoupling{iCh}))],  'FontSize', 6)
	end
	axes('pos',[p.Position(1)+p.Position(3)*.7 p.Position(2)+p.Position(4)*.2 p.Position(3)*.4 p.Position(3)*.4]) % left, bottom, width, height
	pol = polarhistogram(deg2rad(detectedEvents.SloSpiDetCoupling{iCh}));
	pol.Parent.ThetaAxis.Visible = 'off';
	pol.Parent.RAxis.Visible = 'off';
	if iCh == 1
		title(sprintf(['SO Phase detected spindle peak\n n=' num2str(numel(detectedEvents.SloSpiDetCoupling{iCh}))]),  'FontSize', 6, 'Color', 'red')
	else
		title(['n=' num2str(numel(detectedEvents.SloSpiDetCoupling{iCh}))],  'FontSize', 6, 'Color', 'red')
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

