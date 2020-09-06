function output = detectEvents(cfg, data)
% Detects sleep-associated events in EEG / intracranial data (SOs,
% spindles...). Takes a fieldtrip structure with a single trial, returns
% detections.
%
% INPUT VARIABLES:
% cfg
% .scoring						int array (num_epochs x 1)
% .scoring_epoch_length			int; length of scoring epochs in sec
% .code_NREM					int or int array; NREM sleep stages to use for detection
%								(usually [2 3 4] for humans, 2 for animals)
% .code_REM						same for REM sleep
% .code_WAKE					same for wake 
% .artfctdef					artifacts, either a) as an array (num_arts
%								x 2), each row providing [startsample endsample]
% 								of one artifact window, or b) for the lazy
% 								user, artifact definitions as returned by
% 								fieldtrip artifact functions, e.g.
% 								cfg.artfctdef.visual.artifact =	[234 242; 52342 65234];
%								cfg.artfctdef.zvalue.artifact =	[111 222]; ..and so on.
%								Data in these artifactual time windows will
%								not be used. Artifacts may be overlapping.
% .artfctpad					padding (in sec) of segments before and after
%								each artifact to discard (default: 0.5)
% .spi_freq						frequency range for spindle detection
%								(default: [12 16]); if cfg.spi_indiv == 1,
%								this will be the range in which the
%								individual spindle frequency peak is
%								termined; filtering range will then be +/-
%								cfg.spi_indiv_win
% .spi_indiv					logical; turns on the use of individual
%								spindle peak frequencies (calculated within
%								time window provided in *** and on IRASA-
%								computed oscillatory component); default: 0
% .cfg.spi_indiv_win			int (in Hz); signal will be filtered +/-
%								spi_indiv_win around individual spindle
%								peak frequency (default: 2)
% .spi_indiv_chan				cell array with string; channels for
%								estimating spindle peak frequency; you
%								probably dont want to mix frontal and
%								central channels here
% data							Fieldtrip raw data structure, should contain a single trial
%								Should adhere to https://github.com/fieldtrip/fieldtrip/blob/release/utilities/ft_datatype_raw.m
%
% OUTPUT VARIABLES:
% output            ...
%
%
% Todo:
% . also add a hypno function (adding for each sample the sleep stage) + a
% function that tells you for each sample whether its solid in a sleep
% stage: isStage(FT-Structure, sample, boundary_in_sec)
% . Polaritaet checkcen (EEG vs. invasive)
% . describe input and output in documentation
% . add dataset information to .info
% . slow spindles 8-12 hz
% . rework variable naming inside function and output
% . rework output: all data in one row per channel, also per ep and for
% entire recording
% . implement hongi's ripple search for alpha / ripples?
% . add from hongis code: merging of close events?
%
% AUTHORS:
% Niels Niethard, niels.niethard@medizin.uni-tuebingen.de
% Jens Klinzing, klinzing@princeton.edu

%% INPUT VALIDATION AND SETUP
% Check for required fields in cfg
requiredFields = {'scoring', 'scoring_epoch_length', 'code_NREM', 'code_REM', 'code_WAKE'};
for i = requiredFields
	if ~isfield(cfg,i)
		error(['Required field missing in cfg: ' i{1} '.']);
	end
end

% Check data
if length(data.trial) ~= 1, error('Function only accepts single-trial data.'), end
if any(size(data.sampleinfo) ~= [1 2]), error('Sampleinfo looks like data does not contain exactly one trial.'), end
Fs = data.fsample;

% Set default values
if ~isfield(cfg, 'debugging')
	cfg.debugging				= 0; % in s
end
if ~isfield(cfg, 'spi_dur_min')
	cfg.spi_dur_min				= [0.5 0.25]; % in s
end
if ~isfield(cfg, 'spi_dur_max')
	cfg.spi_dur_max				= [2.5 2.5];
end
if ~isfield(cfg, 'spi_thr')
	cfg.spi_thr(1,1)			= 1.5; % nn: 1, 1.5, 2; hongi 1.5 (with rms)
	cfg.spi_thr(2,1)			= 2;
	cfg.spi_thr(3,1)			= 2.5;
elseif isfield(cfg, 'spi_thr') && length(cfg.spi_thr) ~= 3
	error('Spindle thresholds not provided properly (should be 3x1 vector in cfg.spi_thr).')
end
if ~isfield(cfg, 'spi_freq')
	cfg.spi_freq				= [12 16]; % Hz; filtering range for spindle detection
end
if ~isfield(cfg, 'spi_filt_ord')
	cfg.spi_filt_ord			= 6;
end
if ~isfield(cfg, 'spi_indiv') % relatively costly computation
	cfg.spi_indiv		= 0; % if 1, will look for peak between freqs defined in cfg.spi_freq and use +/- cfg.spi_indiv_win instead
end
if size(cfg.scoring, 1) == 1
	cfg.scoring = cfg.scoring';
end
if isfield(cfg, 'artfctdef') && ~isfield(cfg, 'artfctpad')
	cfg.artfctpad = 0.5;
end
if cfg.spi_indiv == 1 && (~isfield(cfg, 'spi_indiv_chan') || isempty(cfg.spi_indiv_chan))
	error('If you want individual spindle peak magic, you gotta provide the channels to use.')
end
if cfg.spi_indiv == 1 && ~isfield(cfg, 'spi_indiv_win')
	cfg.spi_indiv_win = 2; % signal will be filtered +/- spi_indiv_win around individual spindle peak frequency
end
if ~isfield(cfg, 'the_freq')
	cfg.the_freq				= [4 8];
end
if ~isfield(cfg, 'the_filt_ord')
	cfg.the_filt_ord			= 6;
end
if ~isfield(cfg, 'slo_dur_min')
	cfg.slo_dur_min				= 0.5; % in s
end
if ~isfield(cfg, 'slo_dur_max')
	cfg.slo_dur_max				= 2.0;
end
if ~isfield(cfg, 'slo_thr')
	cfg.slo_thr					= 1.5; % in SD; nn: 1, 1.5, 2; hongi 1.5 (with rms)
end
if ~isfield(cfg, 'slo_freq')
	cfg.slo_freq				= [0.1 3.5]; % in Hz
end
if ~isfield(cfg, 'slo_filt_ord')
	cfg.slo_filt_ord			= 3;
end

% Start filling the output
output						= [];
output.info.cfg				= cfg;
output.info.Fs				= Fs;
output.info.length			= size(data.trial{1},2);
output.info.scoring			= cfg.scoring;
output.info.scoring_epoch_length = cfg.scoring_epoch_length;
output.slo					= [];
output.the					= [];
output.spi					= [];

%% PREPARATIONS
chans						= data.label;
multi						= cfg.scoring_epoch_length*Fs;

% Extract artifact time windows in case they are provided in fieldtrip
% format
if isfield(cfg, 'artfctdef')
	if isstruct(cfg.artfctdef)
		fnames = fieldnames(cfg.artfctdef);
		tmp = [];
		for iFn = 1:numel(fnames)
			tmp = [tmp; cfg.artfctdef.(fnames{iFn}).artifact];
		end
		cfg.artfctdef = tmp;
	end
end

% Compensate if scoring and data don't have the same length
tmp_diff					= size(data.trial{1},2) - length(cfg.scoring)*multi;
if tmp_diff < 0 % this should not happen or only be -1
	warning(['Data is shorter than scoring by ' num2str(tmp_diff * -1) ' sample(s). Will act as if I hadn''t seen this.'])
elseif tmp_diff > 0 % scoring is shorter than data (happens e.g., with SchlafAUS)
	data.trial{1}(:, end-(tmp_diff-1):end)	= [];
	data.time{1}(end-(tmp_diff-1):end)	= [];
	data.sampleinfo(2) = data.sampleinfo(2) - tmp_diff;
end
clear tmp_diff

% Create upsampled scoring vector
scoring_fine	= zeros(size(data.trial{1},2),1);
for iEp = 1:length(cfg.scoring)
	scoring_fine((iEp-1)*multi+1 : (iEp)*multi) = cfg.scoring(iEp);
end
output.info.scoring_fine	= scoring_fine; % Let's return the scoring without artifacts

% Mark artifacts in sleep scoring
if isfield(cfg, 'artfctdef')
	warning('Artifact handling is new and should be double-checked.')
	for iArt = 1:size(cfg.artfctdef, 1)
		a_beg = cfg.artfctdef(iArt, 1) -  cfg.artfctpad*Fs;
		if a_beg < 1 % padding shouldnt go to far
			a_beg = 1; 
		end
		if cfg.artfctdef(iArt, 2) > length(scoring_fine)
			error(['Artifact ' num2str(iArt) ' extends outside the data.'])
		end
		a_end = cfg.artfctdef(iArt, 2) +  cfg.artfctpad*Fs;
		if a_end > length(scoring_fine) % shouldnt be too long also after padding
			a_end = length(scoring_fine);
		end
		scoring_fine(a_beg:a_end) = 99;
	end
end


output.info.scoring_artsrem	= scoring_fine;
% also return scoring with artifacts removed

% Extract episodes (save in seconds)
% NREM
NREMBegEpisode = strfind(any(cfg.scoring==cfg.code_NREM,2)',[0 1]); % where does scoring flip to S2
NREMEndEpisode = strfind(any(cfg.scoring==cfg.code_NREM,2)',[1 0]); % where does scoring flip from S2 to something else?
NREMBegEpisode = NREMBegEpisode+1; % because it always finds the epoch before
if any(cfg.scoring(1,1)==cfg.code_NREM,2)
	NREMBegEpisode = [1 NREMBegEpisode];
end
if any(cfg.scoring(end,1)==cfg.code_NREM,2)
	NREMEndEpisode = [NREMEndEpisode length(cfg.scoring)];
end
NREMEpisodes = [(NREMBegEpisode-1)*cfg.scoring_epoch_length+1; NREMEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec

% REM
REMBegEpisode = strfind(any(cfg.scoring==cfg.code_REM,2)',[0 1]);
REMEndEpisode = strfind(any(cfg.scoring==cfg.code_REM,2)',[1 0]);
REMBegEpisode = REMBegEpisode+1;
if any(cfg.scoring(1,1)==cfg.code_REM,2)
	REMBegEpisode = [1 REMBegEpisode];
end
if any(cfg.scoring(end,1)==cfg.code_REM,2)
	REMEndEpisode = [REMEndEpisode length(cfg.scoring)];
end
REMEpisodes = [(REMBegEpisode-1)*cfg.scoring_epoch_length+1; REMEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec

% Wake
WAKBegEpisode = strfind((cfg.scoring==cfg.code_WAKE)',[0 1]);
WAKEndEpisode = strfind((cfg.scoring==cfg.code_WAKE)',[1 0]);
WAKBegEpisode = WAKBegEpisode+1;
if cfg.scoring(1,1) == cfg.code_WAKE
	WAKBegEpisode = [1 WAKBegEpisode];
end
if cfg.scoring(end,1) == cfg.code_WAKE
	WAKEndEpisode = [WAKEndEpisode length(cfg.scoring)];
end
WAKEpisodes = [(WAKBegEpisode-1)*cfg.scoring_epoch_length+1; WAKEndEpisode*cfg.scoring_epoch_length]; %create Matrix with NRem on and offset time in sec

% Fill the output
output.info.channel			= chans;
output.NREMepisode			= NREMEpisodes;
output.REMepisode			= REMEpisodes;
output.WAKEpisodes			= WAKEpisodes;

%% Spindles
% Find individual spindle peaks
if cfg.spi_indiv
	% Cut out NREM episodes
	cfg_tmp				= [];
	cfg_tmp.trl			= [NREMEpisodes'*Fs zeros(size(NREMEpisodes,2),1)];
	data_nrem			= ft_redefinetrial(cfg_tmp, data);
	
	% Partition data into segments (gives more reliable power estimates)
	cfg_tmp				= [];
	cfg_tmp.length		= 4; % should suffice for good spectral resolution
	cfg_tmp.overlap     = 0;
	data_nrem			= ft_redefinetrial(cfg_tmp, data);
	
	% Determine resampling frequency	TODO: Test this
	if Fs > 250
		if mod(Fs, 200) == 0
			res_freq = 200;
		elseif mod(Fs, 250) == 0
			res_freq = 250;
		elseif mod(Fs, 256) == 0
			res_freq = 256;
		else
			error('You got some weird sampling frequency, check out this part of the code and make your own decisions.')
		end
		
		% Resample data (improves performance)
		cfg_tmp				= [];
		cfg_tmp.resamplefs  = res_freq;
		data_nrem			= ft_resampledata(cfg_tmp, data_nrem);
	end
	
	% Perform spectral analysis (both IRASA fractal component and regular spectrum)
	cfg_tmp				= [];
	cfg_tmp.foi			= cfg.spi_freq(1):0.1:cfg.spi_freq(2); % .1 Hz sampling should be fine
	% 	cfg_tmp.foi			= 0.5:0.1:20; % .1 Hz sampling should be fine
	
	cfg_tmp.taper		= 'hanning';
	cfg_tmp.pad			= 'nextpow2';
	cfg_tmp.keeptrials	= 'no';
	cfg_tmp.channel		= cfg.spi_indiv_chan;
	cfg_tmp.method		= 'irasa';
	fra					= ft_freqanalysis(cfg_tmp, data_nrem);
	cfg_tmp.method		= 'mtmfft';
	mix					= ft_freqanalysis(cfg_tmp, data_nrem);
	
	% Average over channels
	cfg_tmp				= [];
	cfg_tmp.avgoverchan = 'yes';
	mix					= ft_selectdata(cfg_tmp,mix);
	fra					= ft_selectdata(cfg_tmp,fra);
	
	% Subtract fractal component from mixed power spectrum
	cfg_tmp				= [];
	cfg_tmp.parameter	= 'powspctrm';
	cfg_tmp.operation   = 'x2-x1';
	osc					= ft_math(cfg_tmp, fra, mix);
	
	% Calculate relative change
	cfg_tmp.operation	= 'divide';
	chan				= ft_math(cfg_tmp, osc, fra);
	
	% Find peaks
	[m,mi] = max(chan.powspctrm);
	spi_freq_indiv = [chan.freq(mi)-cfg.spi_indiv_win chan.freq(mi)+cfg.spi_indiv_win];
	
	% Debugging plots
	if cfg.debugging == 1
		% Use a wider range for power estimate for these plots to make more
		% sense (see comments)
		figure
		subplot(4,3,[1 2])
		plot(fra.freq, mix.powspctrm(1,:))
		title('mixed spectrum')
		xlim([fra.freq(1) fra.freq(end)])
		ylim([0 15])
		subplot(4,3,[4 5])
		plot(fra.freq, fra.powspctrm(1,:))
		title('fractal component')
		xlim([fra.freq(1) fra.freq(end)])
		ylim([0 25])
		subplot(4,3,[7 8])
		plot(fra.freq, osc.powspctrm(1,:))
		title('oscillatory component')
		xlim([fra.freq(1) fra.freq(end)])
		ylim([0 10])
		subplot(4,3,[10 11])
		plot(fra.freq, chan.powspctrm(1,:))
		title('oscillatory component / fractal component')
		xlim([fra.freq(1) fra.freq(end)])
		
		zoom = fra.freq>10 & fra.freq < 18;
		zoomed = fra.freq(zoom);
		
		subplot(4,3,[3])
		plot(zoomed, mix.powspctrm(1,zoom)), hold on
		[m,mi] = max(mix.powspctrm(zoom));
		xline(zoomed(mi))
		title('mixed spectrum')
		xlim([zoomed(1) zoomed(end)])
		subplot(4,3,[6])
		plot(zoomed, fra.powspctrm(1,zoom))
		xlim([zoomed(1) zoomed(end)])
		title('fractal component')
		subplot(4,3,[9])
		plot(zoomed, osc.powspctrm(1,zoom)), hold on
		xlim([zoomed(1) zoomed(end)])
		[m,mi] = max(osc.powspctrm(zoom));
		xline(zoomed(mi))
		title('oscillatory component')
		subplot(4,3,[12])
		plot(zoomed, chan.powspctrm(1,zoom)), hold on
		[m,mi] = max(chan.powspctrm(zoom));
		xline(zoomed(mi))
		xlim([zoomed(1) zoomed(end)])
		title('oscillatory component / fractal component')
	end
end

cfg_pp				= [];
cfg_pp.bpfilter		= 'yes';
if cfg.spi_indiv
	cfg_pp.bpfreq	= spi_freq_indiv;
	output.spi.freq = spi_freq_indiv;
else
	cfg_pp.bpfreq	= cfg.spi_freq;
	output.spi.freq = cfg.spi_freq;
end
cfg_pp.bpfiltord	= cfg.spi_filt_ord;
data_spi			= ft_preprocessing(cfg_pp, data);

spi_amp				= abs(hilbert(data_spi.trial{1}'))'; % needs to be transposed for hilbert, then transposed back...
spi_amp_mean		= mean(spi_amp(:,any(scoring_fine==cfg.code_NREM,2))');
spi_amp_std			= std(spi_amp(:,any(scoring_fine==cfg.code_NREM,2))');

% Detect spindles
spi = cell(size(NREMEpisodes,2),numel(chans)); % each cell will contain a two-row vector with beginning and ends of detected spindles
for iEpoch = 1:size(NREMEpisodes,2)
	spi_amp_tmp = spi_amp(:, NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs);
	for iCh = 1:numel(chans)
		% First threshold criterion
		% Where does the smoothed envelope cross the threshold?
		FastSpiAmplitudeTmp = smooth(spi_amp_tmp(iCh, :),0.1 * Fs); % get smoothed instantaneous amplitude (integer is the span of the smoothing) - !! does almost nothing
		above_threshold = FastSpiAmplitudeTmp > cfg.spi_thr(1,1)*spi_amp_std(iCh); % long column showing threshold crossings
		isLongEnough = bwareafilt(above_threshold, [cfg.spi_dur_min(1)*Fs, cfg.spi_dur_max(1)*Fs]); % find spindle within duration range
		isLongEnough = [0; isLongEnough]; %compensate that spindle might start in the beginning
		SpiBeginning =  strfind(isLongEnough',[0 1]); %find spindle Beginning line before compensates that it find last 0
		SpiEnd = strfind(isLongEnough',[1 0])-1; %find spindle Ending subtract 1 because of added 0 in the beginning
		
		% Some plots for debugging
		if cfg.debugging
			win = 1:5000;
			spi_raw = data_spi.trial{1}(iCh, NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs);
			plot(win/Fs, spi_raw(1,win)), hold on			% raw signal
			plot(win/Fs, spi_amp_tmp(iCh,win), 'r')			% envelope
			plot(win/Fs, FastSpiAmplitudeTmp(win), 'r')		% smoothed envelope
			line([win(1)/Fs win(end)/Fs],[cfg.spi_thr(1,1)*spi_amp_std(iCh) cfg.spi_thr(1,1)*spi_amp_std(iCh)]) % threshold
			plot(win/Fs, above_threshold(win))				% threshold crossed
			plot(win/Fs, isLongEnough(win))					% crosses min-length criterion
		end
		% Delete spindle if it is cut by beginning / end of epoch
		if ~isempty(SpiBeginning) || ~isempty(SpiEnd)
			if length(SpiEnd)<length(SpiBeginning)
				SpiBeginning(:,end)=[];
			end
			if ~isempty(SpiBeginning) || ~isempty(SpiEnd) && SpiBeginning(1,1)==1
				SpiBeginning(:,1) = [];
				SpiEnd(:,1) = [];
			end
			FastSpindles = [SpiBeginning;SpiEnd];
			spi{iEpoch,iCh} = FastSpindles+(NREMEpisodes(1,iEpoch)*Fs);%include beginning of NREMEpoch
		else
			spi{iEpoch,iCh} = [];
		end
		
		CurrentSpindles = spi{iEpoch,iCh};
		TempIdx = [];
		for iSpi = 1: size (CurrentSpindles,2)
			window_size = 5 * Fs; % in sec
			DataTmpSpi = data_spi.trial{1}(iCh, CurrentSpindles(1,iSpi)-window_size : CurrentSpindles(2,iSpi)+window_size); %get filteres spindle signal for eachspindle + - 5sec
			FastSpiAmplitudeTmp = smooth(abs(hilbert(DataTmpSpi)),40);%get smoothed instantaneous amplitude
			
			% Second threshold criterion
			above_threshold = FastSpiAmplitudeTmp(window_size:end-window_size) > cfg.spi_thr(2,1)*spi_amp_std(iCh);
			isLongEnough = bwareafilt(above_threshold, [cfg.spi_dur_min(2)*Fs, cfg.spi_dur_max(2)*Fs]); %find spindle within duration range
			
			% Third threshold criterion
			above_Max = FastSpiAmplitudeTmp(window_size:end-window_size) > cfg.spi_thr(3,1)*spi_amp_std(iCh);
			MaxIsThere = bwareafilt(above_Max, [1, cfg.spi_dur_max(1)*Fs]); %find spindle within duration range
			[pks,locs] = findpeaks(DataTmpSpi(1, window_size:end-window_size),'MinPeakProminence', cfg.spi_thr(1,1)*spi_amp_std(iCh));
			if sum(double(isLongEnough))>1 && sum(double(MaxIsThere))>1 && max(diff(locs))<100 %check if long enough spindle is present and check that no peak to peak distance is more than 125ms
				% do nothing
			else %if criteria not fullfilled store index of Spindles and kill it later
				TempIdx = [TempIdx iSpi];
			end
		end
		spi{iEpoch,iCh}(:,TempIdx)=[];%if not criteriy fullfilled delete detected spindle
	end
end

% Calculate spindel density
output.spi.density = zeros(numel(chans),1);
for iCh = 1:numel(chans)
	TotalNumberOfSpi = 0;
	EpisodeDurations = 0;
	for iEpoch = 1:size(spi,1)
		CurrentSpindles = spi{iEpoch,iCh};
		TotalNumberOfSpi = TotalNumberOfSpi +size(CurrentSpindles,2);
		EpisodeDurations = EpisodeDurations + NREMEpisodes(2,iEpoch)-NREMEpisodes(1,iEpoch);
	end
	output.spi.density(iCh) = TotalNumberOfSpi/(EpisodeDurations/60); %spindle density in spindles per minute
end

% Fill the output
output.spi.events			= cell(numel(chans), 1);
for iCh = 1:size(spi, 2)
	output.spi.events{iCh} = [spi{:,iCh}];
end
output.spi.events_perNREMep	= spi';
output.spi.amp_std			= spi_amp_std;
output.spi.amp_mean			= spi_amp_mean;

clear spi_amp_tmp TotalNumberOfSpi EpisodeDurations spi data_spi

%% SOs
cfg_pp				= [];
cfg_pp.bpfilter		= 'yes';
cfg_pp.bpfreq		= cfg.slo_freq;
cfg_pp.bpfiltord	= cfg.slo_filt_ord;
data_slo			= ft_preprocessing(cfg_pp, data);

% Calculate std for thresholds
slo_raw				= data_slo.trial{1};
slo_std				= std(slo_raw(:,any(scoring_fine==cfg.code_NREM,2))');
slo_mean			= mean(abs(slo_raw(:,any(scoring_fine==cfg.code_NREM,2))'));

% Find negative amplitudes bigger than Threshold
SOEpisodes = cell(numel(chans),1);
for iCh = 1:numel(chans)
	SoThreshold = slo_mean(iCh) + cfg.slo_thr * slo_std(iCh);
	for iEpoch = 1:size(NREMEpisodes,2)
		% Find potential SOs ('episodes')
		slo_tmp			= slo_raw(iCh, NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs)';
		SOBegEpisode	= strfind((slo_tmp<-SoThreshold)',[0 1])-1;
		SOEndEpisode	= strfind((slo_tmp<-SoThreshold)',[1 0]);
		
		% Double-check found events
		if size(SOEndEpisode,1)>0
			if SOEndEpisode(1,1) < SOBegEpisode(1,1)
				SOEndEpisode(:,1) = [];
			end
			if length(SOBegEpisode) > length(SOEndEpisode)
				SOBegEpisode(:,end) = [];
			end
			% Turn within-episode sample into recording sample and add it
			% to result
			SOEpisodes{iCh,1} = [SOEpisodes{iCh,1} [SOBegEpisode+NREMEpisodes(1,iEpoch)*Fs; SOEndEpisode+NREMEpisodes(1,iEpoch)*Fs]];
		end
	end
end

% Check for further characteristics based on zero crossings
ZeroCrossings = cell(numel(chans),1);
for iCh = 1:numel(chans)
	ZeroCrossings{iCh,1} = zeros(3,size(SOEpisodes{iCh,1},2));
	for iEvent = 1:size(SOEpisodes{iCh,1},2)
		X = 0;  % marker for left zero crossing found
		Y = 0;  % marker for right zero crossing found (1) + right plus-to-minus crossing after the upstate (2)
		for iSearchCrossing = 1:2*Fs
			if X == 0 && slo_raw(iCh, SOEpisodes{iCh,1}(1,iEvent)-iSearchCrossing)>0
				ZeroCrossings{iCh,1} (1,iEvent) = SOEpisodes{iCh,1}(1,iEvent)-iSearchCrossing;
				X = 1;
			end
			if Y == 0 && slo_raw(iCh, SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing)>0
				ZeroCrossings{iCh,1} (2,iEvent) = SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing;
				Y = 1;
			end
			if Y ==1 && slo_raw(iCh, SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing)<0
				ZeroCrossings{iCh,1}(3,iEvent) = SOEpisodes{iCh,1}(2,iEvent)+iSearchCrossing;
				Y = 2;
			end
		end
	end
	
	% Delete those events where not all crossings could have been found
	ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(1,:)==0))=[];
	ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(2,:)==0))=[];
	ZeroCrossings{iCh,1}(:,find(ZeroCrossings{iCh,1}(3,:)==0))=[];
	
	% Compensate for cases in which two down peaks lead to the same zero
	% crossings (which are therefore nore unique)
	tmp1 = unique(ZeroCrossings{iCh,1}(1,:));
	tmp2 = unique(ZeroCrossings{iCh,1}(2,:));
	tmp3 = unique(ZeroCrossings{iCh,1}(3,:));
	
	% Hot fix (delete this if above version is fixed): - probably not
	% needed...
	% 	tmp1 = (ZeroCrossings{iCh,1}(1,:));
	% 	tmp2 = (ZeroCrossings{iCh,1}(2,:));
	% 	tmp3 = (ZeroCrossings{iCh,1}(3,:));
	
	ZeroCrossings{iCh,1} = [tmp1; tmp2; tmp3];
	
	%remove SO with to long downstate
	ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(2,:)-ZeroCrossings{iCh,1}(1,:)) > cfg.slo_dur_max*Fs)=[];
	%remove SOs with to short duration
	ZeroCrossings{iCh,1}(:,(ZeroCrossings{iCh,1}(3,:)-ZeroCrossings{iCh,1}(1,:)) < cfg.slo_dur_min*Fs)=[];
	%remove SOs with to small peak to peak amplitude
	Peak2PeakAmp{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
	for iEvent = 1:size(ZeroCrossings{iCh,1},2)
		NegPeakValue = min(slo_raw(iCh, ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent)),[],2);
		PosPeakValue = max(slo_raw(iCh, ZeroCrossings{iCh,1}(2,iEvent):ZeroCrossings{iCh,1}(3,iEvent)),[],2);
		Peak2PeakAmp{iCh,1}(iEvent,1) = abs(NegPeakValue)+PosPeakValue;
	end
	ZeroCrossings{iCh,1}(:,Peak2PeakAmp{iCh,1}<0.07) = [];
	Peak2PeakAmp{iCh,1}(Peak2PeakAmp{iCh,1}<0.07,:) = [];
	
	% Find negative peaks
	NegativePeaks{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
	for iEvent = 1: size(ZeroCrossings{iCh,1},2)
		[M,I] = min(slo_raw(iCh, ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent)),[],2);
		NegativePeaks{iCh,1}(iEvent,1) = ZeroCrossings{iCh,1}(1,iEvent)+I;
	end
	
	% SO-spindle coupling
	if size(output.spi.events{iCh},2) > 0 % if there are spindles in this channel
		% Method 1: Extract SO phase at point of peak amplitude in spindle
		% band (Randolph method) - might fluctuate much and might require
		% more smoothing
		% ...also creates SO waveforms
		twindow						= 2.5; % data will be +/- twindow
		SloSpiAmpCoupling{iCh,1}	= [];
		SOGA{iCh,1}					= zeros(size(NegativePeaks{iCh,1},1),(twindow*2*Fs + 1));
		SOPhase						= [];
		for iSO = 1:size(NegativePeaks{iCh,1},1)
			SOGA{iCh,1}(iSO,:)		= slo_raw(iCh, NegativePeaks{iCh,1}(iSO,1)-(twindow*Fs):NegativePeaks{iCh,1}(iSO,1)+(twindow*Fs));
			SOPhase					= rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:))));
			[~, SpiAmpIndex]		= max(spi_amp(iCh, NegativePeaks{iCh,1}(iSO,1)-(twindow*Fs):NegativePeaks{iCh,1}(iSO,1)+(twindow*Fs))); % spindle maximum amp in samples from SO window start
			SloSpiAmpCoupling{iCh,1}(iSO,1) = SOPhase(SpiAmpIndex);
		end
		
		% Method 2: Extract SO phase based on spindle detection
		SloSpiDetCoupling{iCh,1}	= [];
		SOPhase						= [];
		cnt							= 1;
		for iSO = 1:size(NegativePeaks{iCh,1},1)
			spi_cur = [];
			% If there is a spindle fully inside SO plus minus time
			% window
			spi_ind = find(output.spi.events{iCh}(1,:) > NegativePeaks{iCh,1}(iSO,1)-twindow*Fs & output.spi.events{iCh}(2,:) < NegativePeaks{iCh,1}(iSO,1)+twindow*Fs);
			if size(spi_ind,2) > 0 % if at least one spindle was found
				spi_cur = output.spi.events{iCh}(:,spi_ind);
			end
			for iSpi = 1:size(spi_cur,2)
				SOPhase = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:)))); % SO phase along entire window
				[~, SpiAmpIndex] = max(spi_amp(iCh, spi_cur(1,iSpi):spi_cur(2,iSpi))); % find spindle maximum amp (samples from spindle start)
				tmp = spi_cur(1,iSpi) + SpiAmpIndex - 1; % spindle maximum amp in global samples
				tmp = tmp - (NegativePeaks{iCh,1}(iSO,1)-(twindow*Fs)); % spindle maximum amp in samples from SO window start
				SloSpiDetCoupling{iCh,1}(cnt,1) = SOPhase(tmp); % note down phase there
				cnt = cnt + 1;
			end
		end
	end
end

% add to putput
output.slo.events				= ZeroCrossings; % up-down, down-up, up-down crossings
output.slo.neg_peaks			= NegativePeaks;
output.slo.waveform				= SOGA;
output.SloSpiAmpCoupling		= SloSpiAmpCoupling; % based on spindle amplitude maximum around each slow wave
output.SloSpiDetCoupling		= SloSpiDetCoupling; % similar to above but only if a spindle event was detected

% clear data_slo SOEpisodes NegativePeaks SOGA slo_raw slo_std slo_mean

%% Theta - TODO
%calculated theta power during REM
cfg_pp				= [];
cfg_pp.bpfilter		= 'yes';
cfg_pp.bpfreq		= cfg.the_freq;
cfg_pp.bpfiltord	= cfg.the_filt_ord;
data_the			= ft_preprocessing(cfg_pp, data);
output.the.freq		= cfg.the_freq;

the_raw				= data_the.trial{1};
the_amp				= abs(hilbert(the_raw'))'; % needs to be transposed for hilbert, then transposed back...

output.the.amp_sum	= sum(the_amp(:,any(scoring_fine==cfg.code_REM,2))');
output.the.amp_mean = mean(the_amp(:,any(scoring_fine==cfg.code_REM,2))');
output.the.amp_std	= std(the_amp(:,any(scoring_fine==cfg.code_REM,2))');

output.the.amp_mean_perREMep	= {};
output.the.amp_sum_perREMep		= {};
for iEpoch = 1:size(REMEpisodes,2)
	for iCh = 1:numel(chans)
		output.the.amp_mean_perREMep{iCh,1}(iEpoch,1) = mean(the_amp(iCh, REMEpisodes(1,iEpoch):REMEpisodes(2,iEpoch)));
		output.the.amp_sum_perREMep{iCh,1}(iEpoch,1) = sum(the_amp(iCh, REMEpisodes(1,iEpoch):REMEpisodes(2,iEpoch)));
	end
end

% clear SOGA SOPhase SOSpiCoupling REMThetaMeanAmp REMThetaEnergy spi_amp
