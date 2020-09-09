function output = detectEvents(cfg, data)
% Detects NREM sleep-associated events in EEG or intracranial data
% (SOs/slow waves, spindles...). Also returns the waveforms of detected
% slow waves, the phase in the slow wave band during detected spindles, and
% the theta band amplitude during REM (or any other band if desired). Some
% measures are provided per NREM/REM episode (i.e., per block of connected
% NREM/REM epochs). Takes a fieldtrip structure with a single trial,
% returns detections.
%
% For SOs and spindles, the detection algorithm features some tricks
% specific to these types of events. For example, SOs are detected based on
% 0-crossings in a slow frequency band, their distance etc. Spindles are
% detected not only based on their envelope in the spindle frequency band
% but also by their shape. Furthmore, the individual spindle peak frequency
% of a subject can automatically be determined and detection will be run
% based on that. Additionally, a more generic detection algorithm can be
% used, e.g. to detect alpha bursts during wakefulness or hippocampal
% ripple events. Of course, this algorithm can also be used to detect
% spindles. If you do so, please tell us how these detections compare with
% detections by the other algorithm.
%
% An overview of a subset of the output of this function can be plotted
% using plotDetectedEvents(output).
%
% INPUT VARIABLES:
% data							Fieldtrip raw data structure, should contain a single trial
%								Should adhere to https://github.com/fieldtrip/fieldtrip/blob/release/utilities/ft_datatype_raw.m
% cfg
% .name							string (optional); dataset identifier, will be forwarded to output.info 
% .scoring						int array (num_epochs x 1)
% .scoring_epoch_length			int; length of scoring epochs in sec
% .code_NREM					int or int array; NREM sleep stages to use for detection (usually [2 3 4] for humans, 2 for animals)
% .code_REM						same for REM sleep
% .code_WAKE					same for wake
% .artfctdef					artifacts, either
%								a) as an array (num_arts x 2), each row providing [startsample endsample] of one artifact window, or
% 								b) for the lazy user, artifact definitions as returned by fieldtrip artifact functions, e.g.
% 								cfg.artfctdef.visual.artifact =	[234 242; 52342 65234];
%								cfg.artfctdef.zvalue.artifact =	[111 222]; ...and so on.
%								Data in these artifactual time windows will not be used. Artifacts may be overlapping.
% .artfctpad					int; padding (in sec) of segments to discard before and after to discard (default: 0.5)
% .spectrum						logical; turns estimation of power spectrum on (1) or off (0); default: 0
%								This returns, separately for artifact-free NREM and REM stages, the commonly used raw spectrum (or 'mixed spectrum', mix), as well as the IRASA-computed fractal component (fra), oscillatory component (osc), and their ratio (rel = osc/fra).
%								Note that if spi_indiv = 1, the spectrum is always returned since it has to be calculated for peak detection anyways
%
% Parameters SO/slow wave detection:
% .slo							logical; turns slow oscillation/slow wave detection on (1) or off (0); default: 0
% .slo_dur_min					lower duration threshold; default: 0.5
% .slo_dur_max					upper duration threshold; default: 2.0
% .slo_thr						the STD scaled by this factor will be the amplitude threshold; default: 1.5
% .slo_peak2peak_min            Minimum peak 2 peak amplitude; default: 0.07
% .slo_freq						frequency range in which to perform detection; default: [0.1 3.5]
% .slo_filt_ord					filter order; default: 3
%
% Parameters spindle detection:
% .spi							logical; turns spindle detection on (1) or off (0); default: 0
% .spi_dur_min					array (1 x 2); minimum duration (in sec) for which a spindle must cross the *first* and the *second* amplitude threshold (spi_thr(1,1)); default: [0.5 0.25];
%								Note: Using a second duration minimum that is slightly shorter than the first (together with a second amplitude threshold that is slightly higher than the first) enforces a spindle-typical waxing/waning shape
% .spi_dur_max					array (1 x 2); maximum (in sec) a spindle is allowed to cross the *first* and *second* amplitude threshold (usually these are the same); default: [2.5 2.5]
% .spi_thr(1,1)					the signal amplitude STD scaled by this factor will be the *first* amplitude threshold (events must cross this threshold for at least spi_dur_min(1,1) sec); default: 1.5
% .spi_thr(2,1)					...the *second* threshold (events must cross this threshold for at least spi_dur_min(2,1) sec); default: 2
% .spi_thr(3,1)					...the *third* threshold (events must cross this threshold at least once); default: 2.5
% .spi_freq						int array (1 x 2); frequency range in which to perform detection, default: [12 16])
%								Note: if cfg.spi_indiv == 1 this will be the range in which the individual spindle frequency peak is termined; filtering will then be performed at this frequency +/- cfg.spi_indiv_win
% .spi_filt_ord					order of spindle band filter; default: 6
% .spi_indiv					logical; use individual spindle peak frequencies (calculated based on IRASA-computed oscillatory component, within frequency window provided in spi_indiv_win and on channels in spi_indiv_chan); default: 0
% .spi_indiv_win				int (in Hz); width of frequency window over which to perform spindle detection (individual spindle peak frequency +/- spi_indiv_win); default: 2
% .spi_indiv_chan				cell array with string; channels for estimating spindle peak frequency; you probably dont want to mix far away channels here
%
% Parameters ripple detection:
% .rip_control_Chan             Channel name of the channel used as control for ripple detection.
%
% Parameters theta amplitude:
% .the							logical; turns computation of theta amplitude on (1) or off (0); default: 0
% .the_freq						frequency range in which to perform detection; default: [4 8]
% .the_filt_ord					filter order; default: 3
%
% Parameters generic event detection:
% .gen.***						*** = name of another event type (your choice)
%								For this event type, this function will mostly act as a wrapper for another, more generic event detection algorithm (generic_detector at the end of this file). More than one generic event type can be provided.
%								Working example for ripple detection (results would be added to the output under output.rip.***):
% 								cfg.gen.rip.param.bpfreq			= [70 140];
% 								cfg.gen.rip.param.stageoi			= [2 3 4];
% 								cfg.gen.rip.param.thresType			= 'channelwise';
% 								cfg.gen.rip.param.filtType			= 'fir';
% 								cfg.gen.rip.param.envType			= 'rms';
% 								cfg.gen.rip.param.envWin			= 0.01;
% 								cfg.gen.rip.criterion.duration		= [0.01 0.5];
%								cfg.gen.rip.criterion.variance		= 'centerSD'; % 'centerSD' (default), 'centerMAD' (median absolute derivation), 'scaleCenter', 'scaleEnvSD' or 'percentile' or 'scaleFiltSD'
% 								cfg.gen.rip.criterion.thres			= 2.2; % threshold criterion given as scaling applied to variance parameter
% 								cfg.gen.rip.paramOpt.scndThres		= 2.5; % if > 0 (default = 0), introduces a second (higher) threshold criterion an event has to pass at least once
%								cfg.gen.rip.paramOpt.smoothEnv      = 0.04;
% 								cfg.gen.rip.paramOpt.minCycles      = 1;
% 								cfg.gen.rip.paramOpt.minCyclesNum   = 4;
% 								cfg.gen.rip.paramOpt.mergeEvts		= .02;
%								For more options and documentation, see nested function at end of this file.
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
% . add further dataset information to .info
% . slow spindles 8-12 hz
% . rework variable naming inside function and output
% . rework output: all data in one row per channel, also per ep and for
% entire recording
% . add from hongis code: merging of close events?
% . let people choose whether to compute threshold for each or all channels
% . input range for data must be defined (micro or milli volts) Neuralynx
% creates files with mV!!!
% . SO-check should delete SOs completely (see todo comment)
% . artifact handling! currently, events are detected based on NREM episodes, which are unaffected by artifacts. only std/amp calculations exclude artifact since they are based on scoring_fine, in which artifacts are marked (99).
%   one solution possible: add after each event detection another check for any overlaps with artifacts
%
% AUTHORS:
% Jens Klinzing, klinzing@princeton.edu
% Niels Niethard, niels.niethard@medizin.uni-tuebingen.de (main detection algorithms)
% Hong-Viet V. Ngo, h.ngo@donders.ru.nl & Til Ole Bergmann, til-ole.bergmann@lir-mainz.de (generic detection algorithm)


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
if size(cfg.scoring, 1) == 1, cfg.scoring = cfg.scoring'; end
if size(data.label, 1) == 1, data.label = data.label'; end
Fs = data.fsample;

% Set default values - general
if ~isfield(cfg, 'debugging') % undocumented debugging option
	cfg.debugging				= 0; % in s
end
if isfield(cfg, 'artfctdef') && ~isfield(cfg, 'artfctpad')
	cfg.artfctpad = 0.5;
end

% Set default values - spectrum
if ~isfield(cfg, 'spectrum') % 
	cfg.spectrum				= 0; 
end
if isfield(cfg, 'spi_indiv') && cfg.spi_indiv
	disp('If individual spindle peaks are requested, the spectrum is always calculated.')
	cfg.spectrum				= 1; 
end

% Set default values - slow oscillations/slow waves
if ~isfield(cfg, 'slo')
	cfg.slo						= 0;
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
if ~isfield(cfg, 'slo_peak2peak_min')
	cfg.slo_peak2peak_min       = 0.07; % in same scaling as recoprding!
end
if ~isfield(cfg, 'slo_freq')
	cfg.slo_freq				= [0.1 3.5]; % in Hz
end
if ~isfield(cfg, 'slo_filt_ord')
	cfg.slo_filt_ord			= 3;
end

% Set default values - spindles
if ~isfield(cfg, 'spi')
	cfg.spi						= 0;
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
	cfg.spi_indiv				= 0; % if 1, will look for peak between freqs defined in cfg.spi_freq and use +/- cfg.spi_indiv_win instead
end
if cfg.spi_indiv == 1 && (~isfield(cfg, 'spi_indiv_chan') || isempty(cfg.spi_indiv_chan))
	error('If you want individual spindle peak magic, you gotta provide the channels to use.')
end
if cfg.spi_indiv == 1 && ~isfield(cfg, 'spi_indiv_win')
	cfg.spi_indiv_win			= 2; % signal will be filtered +/- spi_indiv_win around individual spindle peak frequency
end

% Set default values - ripples
if ~isfield(cfg, 'rip')
	cfg.rip						= 0;
end
if ~isfield(cfg, 'rip_indiv')
	cfg.rip_indiv               = 0;
end
if ~isfield(cfg, 'rip_freq')
	cfg.rip_freq				= [150 250];
end
if ~isfield(cfg, 'rip_filt_ord')
	cfg.rip_filt_ord			= 3;
end
if ~isfield(cfg, 'rip_thr')
	cfg.rip_thr					= [2; 5]; % to do: for consistency, all input arrays should be row vectors (also e.g., spi_thr)
end
if ~isfield(cfg, 'rip_dur_min')
	cfg.rip_dur_min				= 0.03;
end
if ~isfield(cfg, 'rip_dur_max')
	cfg.rip_dur_max				= 0.3;
end

% Set default values - theta band
if ~isfield(cfg, 'the')
	cfg.the						= 0;
end
if ~isfield(cfg, 'the_freq')
	cfg.the_freq				= [4 8];
end
if ~isfield(cfg, 'the_filt_ord')
	cfg.the_filt_ord			= 3;
end

% Start filling the output
output						= [];
output.info.cfg				= cfg;
output.info.Fs				= Fs;
output.info.length			= size(data.trial{1},2);
output.info.scoring			= cfg.scoring;
output.info.scoring_epoch_length = cfg.scoring_epoch_length;
output.info.name			= cfg.name;

%% PREPARATIONS
chans						= data.label;
multi						= cfg.scoring_epoch_length*Fs;

% Extract artifact time windows in case they are provided in fieldtrip
% format
if isfield(cfg, 'artfctdef')
	if isstruct(cfg.artfctdef)
		art_types = fieldnames(cfg.artfctdef);
		tmp = [];
		for iAt = 1:numel(art_types)
			tmp = [tmp; cfg.artfctdef.(art_types{iAt}).artifact];
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
output.info.scoring_artsrem	= scoring_fine; % also return scoring with artifacts removed

% Extract episodes (save in seconds)
% NREM
NREMBegEpisode = strfind(any(cfg.scoring==cfg.code_NREM,2)',[0 1]); % where does scoring flip to NREM
NREMEndEpisode = strfind(any(cfg.scoring==cfg.code_NREM,2)',[1 0]); % where does scoring flip from NREM to something else
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

%% Spectrum
% Will be computed on artifact-free data, and separately for NREM and REM
if cfg.spectrum
	disp('Calculating spectrum...')
	spec_freq = [1 45]; % let's not ask the user (to make sure the spindle range is included in this range)
	
	% NREM episodes in sample resolution
	nrem_begs = strfind(any(scoring_fine==cfg.code_NREM,2)',[0 1]); % where does scoring flip to NREM
	nrem_ends = strfind(any(scoring_fine==cfg.code_NREM,2)',[1 0]); % where does scoring flip from NREM to something else
	nrem_begs = nrem_begs+1; % because it always finds the sample before
	if any(scoring_fine(1,1)==cfg.code_NREM,2) % in case recording starts with this stage
		nrem_begs = [1 nrem_begs];
	end
	if any(scoring_fine(end,1)==cfg.code_NREM,2) % in case recording starts with this stage
		nrem_ends = [nrem_ends length(scoring_fine)];
	end
	
	% REM episodes in sample resolution
	rem_begs = strfind(any(scoring_fine==cfg.code_REM,2)',[0 1]); % where does scoring flip to NREM
	rem_ends = strfind(any(scoring_fine==cfg.code_REM,2)',[1 0]); % where does scoring flip from NREM to something else
	rem_begs = rem_begs+1; % because it always finds the sample before
	if any(scoring_fine(1,1)==cfg.code_REM,2) % in case recording starts with this stage
		rem_begs = [1 rem_begs];
	end
	if any(scoring_fine(end,1)==cfg.code_REM,2) % in case recording starts with this stage
		rem_ends = [rem_ends length(scoring_fine)];
	end
	
	% Cut out NREM and REM segemnts
	cfg_tmp						= [];
	cfg_tmp.trl					= [nrem_begs' nrem_ends' zeros(length(nrem_ends'),1)];
	tmp_nrem					= ft_redefinetrial(cfg_tmp, data);
	cfg_tmp.trl					= [rem_begs' rem_ends' zeros(length(rem_ends'),1)];
	tmp_rem						= ft_redefinetrial(cfg_tmp, data);
	
	% Downsample data to speed up spectral estimates (done before cutting
	% in smaller segments is orders of magnitudes faster)
	if Fs > 128
		if mod(Fs, 100) == 0
			res_freq = 100;
		elseif mod(Fs, 125) == 0
			res_freq = 125;
		elseif mod(Fs, 128) == 0
			res_freq = 128;
		else
			error('You got some weird sampling frequency, check out this part of the code and make your own decisions.')
		end
		
		% Resample data (improves performance)
		cfg_tmp				= [];
		cfg_tmp.resamplefs  = res_freq;
		tmp_nrem			= ft_resampledata(cfg_tmp, tmp_nrem);
		tmp_rem				= ft_resampledata(cfg_tmp, tmp_rem);
	end
	
	% Cut into small segments (improves and smoothens spectral estimates)
	cfg_tmp						= [];
	cfg_tmp.length				= 4;  % cut data into segments of this length (in sec)
	cfg_tmp.overlap				= 0;  % with this overlap
	tmp_nrem					= ft_redefinetrial(cfg_tmp, tmp_nrem);
	tmp_rem						= ft_redefinetrial(cfg_tmp, tmp_rem);
	
	% Calculate spectra
	cfg_tmp						= [];
	cfg_tmp.foi					= spec_freq(1):0.05:spec_freq(2);
	cfg_tmp.method				= 'irasa';
	cfg_tmp.pad					= 'nextpow2';
	fra_nrem					= ft_freqanalysis(cfg_tmp, tmp_nrem);
	fra_rem						= ft_freqanalysis(cfg_tmp, tmp_rem);
	
	cfg_tmp.method 				= 'mtmfft';
	cfg_tmp.taper 				= 'hanning';
	mix_nrem					= ft_freqanalysis(cfg_tmp, tmp_nrem);
	mix_rem						= ft_freqanalysis(cfg_tmp, tmp_rem);
	
	% Calculate the oscillatory component by subtracting the fractal from the
	% mixed component
	cfg_tmp						= [];
	cfg_tmp.parameter			= 'powspctrm';
	cfg_tmp.operation			= 'subtract';
	osc_nrem					= ft_math(cfg_tmp, mix_nrem, fra_nrem);
	osc_rem						= ft_math(cfg_tmp, mix_rem, fra_rem);
	
	% Use percent change for even more obvious peaks
	cfg_tmp.operation			= 'divide';
	rel_nrem					= ft_math(cfg_tmp, osc_nrem, fra_nrem);
	rel_rem						= ft_math(cfg_tmp, osc_rem, fra_rem);
	
	output.spectrum.fra_nrem	= fra_nrem.powspctrm;
	output.spectrum.fra_rem		= fra_rem.powspctrm;
	output.spectrum.mix_nrem	= mix_nrem.powspctrm;
	output.spectrum.mix_rem		= mix_rem.powspctrm;
	output.spectrum.osc_nrem	= osc_nrem.powspctrm;
	output.spectrum.osc_rem		= osc_rem.powspctrm;
	output.spectrum.rel_nrem	= rel_nrem.powspctrm;
	output.spectrum.rel_rem		= rel_rem.powspctrm;
	output.spectrum.freq		= fra_nrem.freq; % add frequency vector
	
	if cfg.debugging
		figure
		subplot(3,1,1)
		plot(output.spectrum.freq, output.spectrum.fra_nrem(1,:)), hold on
		plot(output.spectrum.freq, output.spectrum.mix_nrem(1,:))
		subplot(3,1,2)
		plot(output.spectrum.freq, output.spectrum.osc_nrem(1,:))
		subplot(3,1,3)
		plot(output.spectrum.freq, output.spectrum.rel_nrem(1,:))
	end
end


%% Spindles
if cfg.spi
	disp('Starting spindle detection...')
	
	% Find individual spindle peaks
	if cfg.spi_indiv
		% Extract search window, average existing relative spectrum over channels, find maximum
		tmp_freqi		= output.spectrum.freq >= cfg.spi_freq(1) & output.spectrum.freq <= cfg.spi_freq(2);
		tmp_rel			= mean(output.spectrum.rel_nrem(contains(chans, cfg.spi_indiv_chan), tmp_freqi), 1);
		[~,mi]			= max(tmp_rel);
		spi_freq_indiv	= [chan.freq(mi)-cfg.spi_indiv_win chan.freq(mi)+cfg.spi_indiv_win];
		clear tmp_freqi tmp_rel mi
	end
	
	% Filter data in spindle band, extract mean and std
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
				win = 1:50000;
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
			TempIdx = []; % these spindle candidates will be eliminated
			for iSpi = 1: size (CurrentSpindles,2)
				window_size = 5 * Fs; % in sec
				if CurrentSpindles(2,iSpi)+window_size < length(data_spi.trial{1}(iCh,:)) %delete Spi to close to recording end
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
				else
					TempIdx = [TempIdx iSpi];
				end
			end
			spi{iEpoch,iCh}(:,TempIdx)=[];%if not criteriy fullfilled delete detected spindle
		end
	end
	
	% Calculate spindle density
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
end

%% SOs
if cfg.slo
	disp('Starting slow oscillation/slow wave detection...')
	
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
		SOEpisodes{iCh,1} = round(SOEpisodes{iCh,1});%compensate if Fs is not integer
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
		ZeroCrossings{iCh,1}(:,Peak2PeakAmp{iCh,1}<cfg.slo_peak2peak_min) = [];
		Peak2PeakAmp{iCh,1}(Peak2PeakAmp{iCh,1}<cfg.slo_peak2peak_min,:) = [];
		
		% Find negative peaks
		NegativePeaks{iCh,1} = zeros(size(ZeroCrossings{iCh,1},2),1);
		for iEvent = 1: size(ZeroCrossings{iCh,1},2)
			[M,I] = min(slo_raw(iCh, ZeroCrossings{iCh,1}(1,iEvent):ZeroCrossings{iCh,1}(2,iEvent)),[],2);
			NegativePeaks{iCh,1}(iEvent,1) = ZeroCrossings{iCh,1}(1,iEvent)+I;
		end
		
		% Calculate waveforms
		twindow						= 2.5; % data will be +/- twindow
		SOGA{iCh,1}					= zeros(size(NegativePeaks{iCh,1},1),round(twindow*2*Fs + 1));
		for iSO = 1:size(NegativePeaks{iCh,1},1)
			SOGA{iCh,1}(iSO,:)		= slo_raw(iCh, NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs):NegativePeaks{iCh,1}(iSO,1)+round(twindow*Fs));
		end
		
		% SO-spindle coupling
		if cfg.spi && size(output.spi.events{iCh},2) > 0 % if there are spindles in this channel
			% Method 1: Extract SO phase at point of peak amplitude in spindle
			% band (Randolph method) - might fluctuate much and might require
			% more smoothing
			SloSpiAmpCoupling{iCh,1}	= [];
			SOPhase						= [];
			for iSO = 1:size(NegativePeaks{iCh,1},1)
				% Todo: The check whether SOs are too far at the end should be performed further up, such SOs should completely be deleted and show up nowhere
				if NegativePeaks{iCh,1}(iSO,1)+round(twindow*Fs) < length(slo_raw) %only consider SOs far enough from recording end
					SOPhase					= rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:))));
					[~, SpiAmpIndex]		= max(spi_amp(iCh, NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs):NegativePeaks{iCh,1}(iSO,1)+round(twindow*Fs))); % spindle maximum amp in samples from SO window start
					SloSpiAmpCoupling{iCh,1}(iSO,1) = SOPhase(SpiAmpIndex);
				end
			end
			
			% Method 2: Extract SO phase based on spindle detection
			SloSpiDetCoupling{iCh,1}	= [];
			SOPhase						= [];
			cnt							= 1;
			for iSO = 1:size(NegativePeaks{iCh,1},1)
				spi_cur = [];
				% If there is a spindle fully inside SO plus minus time
				% window
				spi_ind = find(output.spi.events{iCh}(1,:) > NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs) & output.spi.events{iCh}(2,:) < NegativePeaks{iCh,1}(iSO,1)+round(twindow*Fs));
				if size(spi_ind,2) > 0 % if at least one spindle was found
					spi_cur = output.spi.events{iCh}(:,spi_ind);
				end
				for iSpi = 1:size(spi_cur,2)
					SOPhase = rad2deg(angle(hilbert(SOGA{iCh,1}(iSO,:)))); % SO phase along entire window
					[~, SpiAmpIndex] = max(spi_amp(iCh, spi_cur(1,iSpi):spi_cur(2,iSpi))); % find spindle maximum amp (samples from spindle start)
					tmp = spi_cur(1,iSpi) + SpiAmpIndex - 1; % spindle maximum amp in global samples
					tmp = tmp - (NegativePeaks{iCh,1}(iSO,1)-round(twindow*Fs)); % spindle maximum amp in samples from SO window start
					SloSpiDetCoupling{iCh,1}(cnt,1) = SOPhase(round(tmp)); % note down phase there
					cnt = cnt + 1;
				end
			end
		end
	end
	
	% add to output
	output.slo.events				= ZeroCrossings; % up-down, down-up, up-down crossings
	output.slo.neg_peaks			= NegativePeaks;
	output.slo.waveform				= SOGA;
	if cfg.spi
		output.SloSpiAmpCoupling		= SloSpiAmpCoupling; % based on spindle amplitude maximum around each slow wave
		output.SloSpiDetCoupling		= SloSpiDetCoupling; % similar to above but only if a spindle event was detected
	end
	% clear data_slo SOEpisodes NegativePeaks SOGA slo_raw slo_std slo_mean
end

%% Ripple detection
if cfg.rip
	disp('Starting ripple detection...')
	
	cfg_pp				= [];
	cfg_pp.bpfilter		= 'yes';
	cfg_pp.bpfreq	= cfg.rip_freq;
	output.rip.freq = cfg.rip_freq;
	cfg_pp.bpfiltord	= cfg.rip_filt_ord;
	data_rip			= ft_preprocessing(cfg_pp, data);
	
	rip_amp				= abs(hilbert(data_rip.trial{1}'))'; % needs to be transposed for hilbert, then transposed back...
	rip_amp_mean		= mean(rip_amp(:,any(scoring_fine==cfg.code_NREM,2))');
	rip_amp_std			= std(rip_amp(:,any(scoring_fine==cfg.code_NREM,2))');
	
	rip = cell(size(NREMEpisodes,2),numel(chans)); % each cell will contain a two-row vector with beginning and ends of detected ripples
	for iEpoch = 1:size(NREMEpisodes,2)
		rip_amp_tmp = rip_amp(:, NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs);
		for iCh = 1:numel(chans)
			% First threshold criterion for min duration
			% Where does the smoothed envelope cross the threshold?
			RipAmplitudeTmp = smooth(rip_amp_tmp(iCh, :),0.004 * Fs); % get smoothed instantaneous amplitude (integer is the span of the smoothing) - !! does almost nothing
			above_threshold = RipAmplitudeTmp > cfg.rip_thr(1,1)*rip_amp_std(iCh); % long column showing threshold crossings
			isLongEnough = bwareafilt(above_threshold, [cfg.rip_dur_min(1)*Fs, cfg.rip_dur_max(1)*Fs]); % find ripple within duration range
			isLongEnough = [0; isLongEnough]; %compensate that ripple might start in the beginning
			ripBeginning =  strfind(isLongEnough',[0 1]); %find ripple Beginning line before compensates that it find last 0
			ripEnd = strfind(isLongEnough',[1 0])-1; %find ripple Ending subtract 1 because of added 0 in the beginning
			
			% Some plots for debugging
			if cfg.debugging
				win = 1:50000;
				rip_raw = data_rip.trial{1}(iCh, NREMEpisodes(1,iEpoch)*Fs : NREMEpisodes(2,iEpoch)*Fs);
				plot(win/Fs, rip_raw(1,win)), hold on			% raw signal
				plot(win/Fs, rip_amp_tmp(iCh,win), 'r')			% envelope
				plot(win/Fs, RipAmplitudeTmp(win), 'r')		% smoothed envelope
				line([win(1)/Fs win(end)/Fs],[cfg.rip_thr(1,1)*rip_amp_std(iCh) cfg.rip_thr(1,1)*rip_amp_std(iCh)]) % threshold
				plot(win/Fs, above_threshold(win))				% threshold crossed
				plot(win/Fs, isLongEnough(win))					% crosses min-length criterion
			end
			% Delete ripple if it is cut by beginning / end of epoch
			if ~isempty(ripBeginning) || ~isempty(ripEnd)
				if length(ripEnd)<length(ripBeginning)
					ripBeginning(:,end)=[];
				end
				if ~isempty(ripBeginning) || ~isempty(ripEnd) && ripBeginning(1,1)==1
					ripBeginning(:,1) = [];
					ripEnd(:,1) = [];
				end
				ripples = [ripBeginning;ripEnd];
				rip{iEpoch,iCh} = ripples+(NREMEpisodes(1,iEpoch)*Fs);%include beginning of NREMEpoch
			else
				rip{iEpoch,iCh} = [];
			end
			
			CurrentRipples = rip{iEpoch,iCh};
			TempIdx = [];
			for irip = 1: size (CurrentRipples,2)
				window_size = 0.5 * Fs; % in sec
				if  CurrentRipples(2,irip)+window_size < length(data_rip.trial{1}(iCh,:)) %check for distance to recording end
					DataTmprip = data_rip.trial{1}(iCh, CurrentRipples(1,irip)-window_size : CurrentRipples(2,irip)+window_size); %get filteres ripple signal for eachripple + - 5sec
					RipAmplitudeTmp = smooth(abs(hilbert(DataTmprip)),40);%get smoothed instantaneous amplitude
					
					% Second Peak threshold criterion
					above_Max = RipAmplitudeTmp(window_size:end-window_size) > cfg.rip_thr(2,1)*rip_amp_std(iCh);
					MaxIsThere = bwareafilt(above_Max, [1, cfg.rip_dur_max(1)*Fs]); %find ripple within duration range
					[pks,locs] = findpeaks(DataTmprip(1, window_size:end-window_size),'MinPeakProminence', cfg.rip_thr(1,1)*rip_amp_std(iCh));
					if sum(double(isLongEnough))>1 && sum(double(MaxIsThere))>1 && max(diff(locs))<100 %check if long enough ripple is present and check that no peak to peak distance is more than 125ms
						% do nothing
					else %if criteria not fullfilled store index of ripples and kill it later
						TempIdx = [TempIdx irip];
					end
					if isfield(cfg,'rip_control_Chan') && strcmp(data_rip.label{iCh},cfg.rip_control_Chan)%check for detected common noise in control channel
						CurrentControlRipples = rip{iEpoch,strcmp(data_rip.label,cfg.rip_control_Chan)};
						if any(ismember(CurrentControlRipples(1,:),CurrentRipples(1,irip):CurrentRipples(2,irip)))||... %check if control ripple Beginning is inside detected ripple
								any(ismember(CurrentControlRipples(2,:),CurrentRipples(1,irip):CurrentRipples(2,irip)))  %check if control ripple Ending is inside detected ripple
							TempIdx = [TempIdx irip];
						end
					end
				else %if ripple to close to recording end
					TempIdx = [TempIdx irip];
				end
			end
			
			rip{iEpoch,iCh}(:,TempIdx)=[];%if not criteriy fullfilled delete detected ripple
		end
	end
	
	% Calculate ripple density
	output.rip.density = zeros(numel(chans),1);
	for iCh = 1:numel(chans)
		TotalNumberOfRip = 0;
		EpisodeDurations = 0;
		for iEpoch = 1:size(rip,1)
			CurrentRipples = rip{iEpoch,iCh};
			TotalNumberOfRip = TotalNumberOfRip +size(CurrentRipples,2);
			EpisodeDurations = EpisodeDurations + NREMEpisodes(2,iEpoch)-NREMEpisodes(1,iEpoch);
		end
		output.rip.density(iCh) = TotalNumberOfRip/(EpisodeDurations/60); %ripple density in ripples per minute
	end
	
	% Fill the output
	output.rip.events			= cell(numel(chans), 1);
	for iCh = 1:size(rip, 2)
		output.rip.events{iCh} = [rip{:,iCh}];
	end
	output.rip.events_perNREMep	= rip';
	output.rip.amp_std			= rip_amp_std;
	output.rip.amp_mean			= rip_amp_mean;
	
	clear rip_amp_tmp TotalNumberOfrip EpisodeDurations rip data_rip
end

%% Theta
% Calculates theta amplitude during REM
if cfg.the
	disp('Starting computation of theta amplitude...')
	
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
end

%% Additional generic events
if isfield(cfg, 'gen') && ~isempty(cfg.gen)
	data.staging = {scoring_fine};
	event_types = fieldnames(cfg.gen);
	for iEt = 1:numel(event_types)
		if cfg.gen.(event_types{iEt}).param.bpfreq(2) > data.fsample / 2 % check this here because it may lead to less intelligable errors later
			error(['Don''t ask for events that extends beyond the Nyquist frequency (' event_types{iEt} ').'])
		end
		output.(event_types{iEt}) = generic_detector(cfg.gen.(event_types{iEt}), data);
	end
end
end

%% Local functions for generic event detection
function cfg_gen = generic_detector(cfg_gen, inData)
% Generic detection function, based on code by Hong-Viet V. Ngo and Til
% Bergmann.
%
% Configuration parameters
% cfg.chThres           = (N x 1) string listing the N channels to determine the thresholds for; default: all
% cfg.chDetct           = (M x 1) string listing the M (< N) channels to detect SO in (based on the previously determined thresholds); default: all
% cfg.param.stageoi     = (N_SlSt x 1) vector containing sleep stages of interest for detection
% cfg.param.bpfreq      = [freqMin, freqMax]: Lower and upper limit for bandpass filter
% cfg.param.filtType    = 'fir' or 'but' (default) - determines whether data processingis based on a FIR or butterworth filter
% cfg.param.thresType   = 'channelwise' (default), 'average' or 'fixed' (not implemented)
%                         Use thresholds corresponding to each channel, a mean threshold
%                         across all channels or given fixed values. If channelwise is chosen
%                         detectCh must be a subset of chThres.
% cfg.param.envType     = 'rms' (default), 'analytic' (Hilbert envelope); perform thresholding on rms signal or envelope (via spline interpolation of the peaks in the rectified signal) of bandpassed signal
% cfg.param.envWin      = Window (in s) to use for envelope calculation (default 0.2 s)
% cfg.param.artfctPad   = [prePad, postPad] (in s) additional padding around each
%                         artifact, negative values = pre-artifact padding (default = [0 0])
%
% cfg.criterion.duration    = [minLen, maxLen]: Minimal and maximal allowed duration of slow oscillation events
% cfg.criterion.center      = 'mean' (default) or 'median'
% cfg.criterion.variance    = 'centerSD' (default), 'centerMAD' (median absolute derivation), 'scaleCenter', 'scaleEnvSD' or 'percentile' or 'scaleFiltSD'
% cfg.criterion.thres     = scalar representing the scaling factor applied to the 'variance' parameter
% cfg.criterion.padding     = [prePad, postPad] (in s) additional padding around
%                             each event candidate, negative value = pre-event
%                             padding (default = [0 0])

% cfg.paramOpt.smoothEnv    = If > 0 (default = 0), window (in s) to use for smoothing of envelope signal
% cfg.paramOpt.upperCutoff  = Scaling factor for an additional threshold, which discard any exceeding events to avoid artifacts (default = inf)
%                             Corresponding threshold is saved in cfg.thres.upperCutoff
% cfg.paramOpt.mergeEvts    = if > 0, maximal gap (in s) until which close events are merged
% cfg.paramOpt.scndThres    = if > 0 (default = 0), introduces a second (higher) threshold criterion an event has to pass at least once
% cfg.paramOpt.minCycles    = Option to activate minimum number of required cycles per event
%                             0: Deactivated (default)
%                             1: Based on raw signal
%                             2: Based on bandpass filtered signal
% cfg.paramOpt.minCyclesNum = Scalar representing the minimum number of required cycles
%
% cfg.doFalseposRjct  = set to 1 if detected events are checked for their
%                         frequency profile, i.e. a spectral peak with a
%                         specific prominence within a specified frequency range
% cfg.falseposRjct.freqlim    = [freqMin freqMax], frequency range event
%                                 should have a spectral maximum
% cfg.falseposRjct.timePad    = [timeMin timeMax], time padded around
%                                 events to calculate TFRs
% cfg.falseposRjct.tfrFreq    = Frequencies for TFR calculation. Must
%                                 include freqlim set above
% cfg.falseposRjct.tfrTime    = Time range of TFR
% cfg.falseposRjct.tfrWin     = Length of time window for TFR calculation
% cfg.falseposRjct.avgWin     = Time window used for averaging TFR, usual
%                                 narrowly set around the event, i.e. t = 0 s
% cfg.falseposRjct.prominence = Threshold the spectral peak amplitude has to exceed
%
% example input structure for two channel analyis:
% inData.label:       {'Fz'; 'Cz'}
% inData.time:        {[2?37320000 double]}
% inData.trial:       {[2?37320000 double]}
% inData.staging:     {[1?37320000 int8]} -- dont worry, this is done by wrapper function
% inData.fsample:     1000
% inData.artifacts:   {[46020?2 double]; [50020?2 double]} -- dont worry, this is done by wrapper function
% inData.sampleinfo:  [1 37320000]

%% Housekeeping
%.. check input
if nargin ~= 2;                     error('Wrong number of arguments');     end
if ~isfield(inData, 'staging');     error('No sleep staging available');    end

%.. check input channels
if ~isfield(cfg_gen,'chThres') || ismember('all', cfg_gen.chThres);     cfg_gen.chThres = inData.label;     end
if ~isfield(cfg_gen,'chDetct') || ismember('all', cfg_gen.chDetct);     cfg_gen.chDetct = inData.label;     end

if strcmp(cfg_gen.param.thresType, 'channelwise') %% Check if detectCh is a subset of chThres
	if sum(ismember(cfg_gen.chThres, cfg_gen.chDetct)) ~= size(cfg_gen.chDetct)
		error('Channelwise detection requires detectCh to be a subset of chThres');
	end
end

%.. check criterion struct
if ~isfield(cfg_gen.criterion,'duration');  error('no duration specfied');          end
if ~isfield(cfg_gen.criterion,'center');    cfg_gen.criterion.center    = 'mean';       end
if ~isfield(cfg_gen.criterion,'variance');  cfg_gen.criterion.variance  = 'centerSD';   end
if ~isfield(cfg_gen.criterion,'scaling');   cfg_gen.criterion.thres   = 1;            end
if ~isfield(cfg_gen.criterion,'padding');   cfg_gen.criterion.padding   = [0 0];        end

%.. check param struct
if ~isfield(cfg_gen.param,'bpfreq');        error('band pass filter information is missing!');  end
if ~isfield(cfg_gen.param,'stageoi');       error('please specify sleep stages of interest!');  end
if ~isfield(cfg_gen.param,'filtType');      cfg_gen.param.filtType  = 'fir';                        end
if ~isfield(cfg_gen.param,'envType');       cfg_gen.param.envType   = 'rms';                        end
if ~isfield(cfg_gen.param,'envWin');        cfg_gen.param.envWin    = 0.2;                          end
if ~isfield(cfg_gen.param,'artfctPad');     cfg_gen.param.artfctPad = [0 0];                        end

%.. check paramOpt struct
if ~isfield(cfg_gen.paramOpt,'smoothEnv');     cfg_gen.paramOpt.smoothEnv = 0;      end
if ~isfield(cfg_gen.paramOpt,'scndThres');     cfg_gen.paramOpt.scndThres = 0;      end
if ~isfield(cfg_gen.paramOpt,'upperCutoff');   cfg_gen.paramOpt.upperCutoff = inf;  end
if ~isfield(cfg_gen.paramOpt,'mergeEvts');      cfg_gen.paramOpt.mergeEvts = 0;     end

% Check minCycle option
if ~isfield(cfg_gen.paramOpt,'minCycles')
	cfg_gen.paramOpt.minCycles      = 0;
	cfg_gen.paramOpt.minCyclesNum   = 0;
elseif cfg_gen.paramOpt.minCycles > 0 && ~isfield(cfg_gen.paramOpt, 'minCyclesNum')
	error('Specification of minCyclesNum is missing');
end

%.. false positive rejection
if ~isfield(cfg_gen,'doFalseposRjct'); cfg_gen.doFalseposRjct = 0; end

% Prepare output structure and some important variables
cfg_gen.thres       = [];
cfg_gen.evtIndiv    = struct;
cfg_gen.summary     = struct;

fsample     = round(inData.fsample);
numTrl      = size(inData.trial,2);
lenTrl      = diff(inData.sampleinfo,1,2)+1;
numThres    = size(cfg_gen.chThres,1);
numDetct    = size(cfg_gen.chDetct,1);
artfctPad   = round(cfg_gen.param.artfctPad * fsample);
critPad     = round(cfg_gen.criterion.padding * fsample);

%% Prepare artfctFilter -- redundant, taken care of by wrapper function
godfltr = cellfun(@(x) ones(numThres,size(x,2),'logical'),inData.trial,'UniformOutput',0);
if isfield(inData, 'artfcts')
	for iTrl = 1 : numTrl
		for iCh = 1 : numThres
			currCh = ismember(inData.label,cfg_gen.chThres{iCh});
			tmpArt = inData.artifacts{currCh,iTrl} + artfctPad;
			tmpArt(tmpArt(:,1) < 1,1)                         = 1;                              %% Ensure padding does not...
			tmpArt(tmpArt(:,2) > inData.sampleinfo(iTrl,2),2) = inData.sampleinfo(iTrl,2);      %% exceed data range
			godfltr{iTrl}(iCh,:) = all([ismember(inData.staging{iTrl},cfg_gen.param.stageoi);...
				~hvn_createBnrySignal(tmpArt,lenTrl(iTrl))]);
			clear tmpArt
		end
	end
end

%% Calculate envelope and determine amplitude threshold
thres       = arrayfun(@(x) nan(numThres,x),lenTrl,'UniformOutput',0);
thresEnv    = arrayfun(@(x) nan(numThres,x),lenTrl,'UniformOutput',0);

% Filtering within desired frequency range
for iTrl = 1 : numTrl                                                   % loop over trials
	for iCh = 1 : numThres                                            % loop over channels
		fprintf('Filter channel %d/%d (%s)\n', iCh, numThres, cfg_gen.chThres{iCh});
		currCh = ismember(inData.label,cfg_gen.chThres{iCh});
		%.. filtering
		switch cfg_gen.param.filtType
			case 'fir'
				thres{iTrl}(iCh,:) = ft_preproc_bandpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg_gen.param.bpfreq, 3*fix(fsample/cfg_gen.param.bpfreq(1))+1, 'fir', 'twopass');
			case 'but'
				tmpHP = ft_preproc_highpassfilter(inData.trial{iTrl}(currCh,:), fsample, cfg_gen.param.bpfreq(1), 5, 'but', 'twopass','reduce');
				thres{iTrl}(iCh,:) = ft_preproc_lowpassfilter(tmpHP, fsample, cfg_gen.param.bpfreq(2), 5, 'but', 'twopass','reduce');
				clear tmpHP
		end
		
		%.. envelope
		thresEnv{iTrl}(iCh,:) = envelope(thres{iTrl}(iCh,:),round(cfg_gen.param.envWin * fsample),cfg_gen.param.envType);
		
		%.. optional smoothing
		if cfg_gen.paramOpt.smoothEnv > 0 %% Smooth envelope if requested
			thresEnv{iTrl}(iCh,:) = smoothdata(thresEnv{iTrl}(iCh,:),2,'movmean',round(cfg_gen.paramOpt.smoothEnv * fsample));
		end
	end
end

% Collect data from all trials
vardata = cell(numThres,1);
for iTrl = 1 : numTrl                           % loop over trials
	if strcmp(cfg_gen.criterion.variance,'scaleFiltSD')
		for iCh = 1 : numThres
			vardata{iCh} = [vardata{iCh}, thres{iTrl}(iCh,godfltr{iTrl}(iCh,:))];
		end
	else
		for iCh = 1: numThres
			vardata{iCh} = [vardata{iCh}, thresEnv{iTrl}(iCh,godfltr{iTrl}(iCh,:))];
		end
	end
end

% Determine center and threshold in one of many different ways
switch cfg_gen.criterion.center
	case 'median'
		centerfun = @(x) nanmedian(x,2);
	case 'mean'
		centerfun = @(x) nanmean(x,2);
end
cfg_gen.thres.center = cell2mat(cellfun(@(x) centerfun(x), vardata,'UniformOutput',0));

if strcmp(cfg_gen.criterion.variance,'centerMAD')
	varfun = @(x) mad(x,1,2,'omitnan');
elseif any(contains(cfg_gen.criterion.variance,{'centerSD','scaleEnvSD','scaleFiltSD'}))
	varfun = @(x) nanstd(x,1,2);
else
	varfun = @(x) [];
end
cfg_gen.thres.variance = cell2mat(cellfun(@(x) varfun(x),vardata,'UniformOutput',0));

if any(contains(cfg_gen.criterion.variance,{'centerMAD','centerSD'}))
	thresfun = @(x,y,z) x + (y .* z);
elseif any(contains(cfg_gen.criterion.variance,{'scaleEnvSD','scaleFiltSD'}))
	thresfun = @(x,y,z) y .* z;
elseif strcmp(cfg_gen.criterion.variance,'scaleCenter')
	thresfun = @(x,y,z) x .* z;
end

if strcmp(cfg_gen.criterion.variance,'percentile')
	cfg_gen.thres.main          = cellfun(@(x) prctile(x, cfg_gen.criterion.thres,2), vardata);
	cfg_gen.thres.upperCutoff   = cellfun(@(x) prctile(x, cfg_gen.paramOpt.upperCutoff,2), vardata);
	if cfg_gen.paramOpt.scndThres > 0
		cfg_gen.thres.second = cellfun(@(x) prctile(x, cfg_gen.paramOpt.scndThres,2), vardata);
	end
else
	cfg_gen.thres.main           = thresfun(cfg_gen.thres.center, cfg_gen.thres.variance, cfg_gen.criterion.thres);
	cfg_gen.thres.upperCutoff    = thresfun(cfg_gen.thres.center, cfg_gen.thres.variance, cfg_gen.paramOpt.upperCutoff);
	if cfg_gen.paramOpt.scndThres > 0
		cfg_gen.thres.second = thresfun(cfg_gen.thres.center, cfg_gen.thres.variance, cfg_gen.paramOpt.scndThres);
	end
end

% Optional: Deal with upper cutoff
% In case an additional upper cutoff was set to prevent artifacts (set
% samples upperCutoff to Nan and re-calculate the threshold)
if cfg_gen.paramOpt.upperCutoff < Inf
	for iCh = 1 : numThres
		vardata{iCh}(vardata{iCh} > cfg_gen.thres.upperCutoff(iCh)) = nan;
	end
	
	cfg_gen.thres.center      = cell2mat(cellfun(@(x) centerfun(x),vardata,'UniformOutput',0));
	cfg_gen.thres.variance    = cell2mat(cellfun(@(x) varfun(x),vardata,'UniformOutput',0));
	
	if strcmp(cfg_gen.criterion.variance,'percentile')
		cfg_gen.thres.main = cellfun(@(x) prctile(x, cfg_gen.criterion.thres,2), vardata);
		
		if cfg_gen.paramOpt.scndThres > 0
			cfg_gen.thres.second = cellfun(@(x) prctile(x, cfg_gen.paramOpt.scndThres,2), vardata);
		end
	else
		cfg_gen.thres.main = thresfun(cfg_gen.thres.center, cfg_gen.thres.variance, cfg_gen.criterion.thres);
		
		if cfg_gen.paramOpt.scndThres > 0
			cfg_gen.thres.second = thresfun(cfg_gen.thres.center, cfg_gen.thres.variance, cfg_gen.paramOpt.scndThres);
		end
	end
end

% If desired calculate and use average threshold across channels
if strcmp(cfg_gen.param.thresType, 'average')
	cfg_gen.thres.main = repmat(mean(cfg_gen.thres.main,1),numThres,1);
	if cfg_gen.paramOpt.scndThres > 0
		cfg_gen.thres.second = repmat(mean(cfg_gen.thres.second,1),numThres,1);
	end
end

%% Event detection
% Collects all samples that are above the threshold and under the upper cutoff
supThres = arrayfun(@(x) nan(numDetct,x),lenTrl,'UniformOutput',0);
for iTrl = 1 : numTrl % loop over trials
	for iCh = 1 : numDetct
		currCh = ismember(cfg_gen.chThres, cfg_gen.chDetct{iCh});  % match current detection channel to chThres vector
		supThres{iTrl}(iCh,:) = all([thresEnv{iTrl}(currCh,:) >= cfg_gen.thres.main(currCh);...
			thresEnv{iTrl}(currCh,:) <= cfg_gen.thres.upperCutoff(currCh);...
			godfltr{iTrl}(currCh,:)]);
	end
end

% Check each event for additional requirements and calculate metrics
for iTrl = 1 : numTrl               % loop over trials
	for iCh = 1 : numDetct       % loop over detect channels
		currCh = ismember(cfg_gen.chThres, cfg_gen.chDetct{iCh});
		cfg_gen.evtIndiv(iCh,iTrl).label   = cfg_gen.chDetct{iCh};   % Channel name
		cfg_gen.evtIndiv(iCh,iTrl).tss     = sum(godfltr{iTrl}(currCh,:),2) / (fsample * 60);  % time spend asleep (in min), based on artifact-free sleep % ?? Why is that calculated per channel and why here
		
		tmpEvts = hvn_extrctBnryBouts(supThres{iTrl}(iCh,:)); % find events (bouts of connected threshold crossing)
		
		% Optional: Second thresholding
		if cfg_gen.paramOpt.scndThres > 0
			for kEvt = 1 : size(tmpEvts,1)
				if max(thresEnv{iTrl}(currCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2))) < cfg_gen.thres.second(currCh)
					supThres{iTrl}(iCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2)) = 0;
				end
			end
			tmpEvts = hvn_extrctBnryBouts(supThres{iTrl}(iCh,:)); % redo event search
		end
		
		% Discard intervals not fulfilling minimal length
		rmvIdx  = hvn_createBnrySignal(tmpEvts(diff(tmpEvts,1,2) < round(cfg_gen.criterion.duration(1) * fsample),:),lenTrl(iTrl));
		supThres{iTrl}(iCh,rmvIdx) = 0;
		tmpEvts = hvn_extrctBnryBouts(supThres{iTrl}(iCh,:)); % redo event search
		% ?? Why dont we discard too long events right here, too?
		
		% Optional: Merge intervals closer than specified margin
		if cfg_gen.paramOpt.mergeEvts > 0
			for kEvt = 1 : size(tmpEvts,1)-1
				if tmpEvts(kEvt+1,1) - tmpEvts(kEvt,2) <= round(cfg_gen.paramOpt.mergeEvts * fsample) && ...
						all(godfltr{iTrl}(currCh,tmpEvts(kEvt,1):tmpEvts(kEvt+1,2)))
					supThres{iTrl}(iCh,tmpEvts(kEvt,1):tmpEvts(kEvt+1,2)) = 1;
				end
			end
			tmpEvts = hvn_extrctBnryBouts(supThres{iTrl}(iCh,:)); % redo event search
		end
		
		% Optional: Discard events not fulfilling minimal number of cycles
		if cfg_gen.paramOpt.minCycles > 0
			switch cfg_gen.paramOpt.minCycles
				case 1
					currSig = smoothdata(inData.trial{iTrl}(ismember(inData.label,cfg_gen.chDetct{iCh}),:),'movmedian',3);
				case 2
					currSig = smoothdata(thres{iTrl}(currCh,:),'movmedian',3);
			end
			for kEvt = 1 : size(tmpEvts,1)
				[~,maxidx] = findpeaks(currSig(1,tmpEvts(kEvt,1):tmpEvts(kEvt,2)));
				[~,minidx] = findpeaks((-1) * currSig(1,tmpEvts(kEvt,1):tmpEvts(kEvt,2)));
				if (numel(maxidx) < cfg_gen.paramOpt.minCyclesNum) || (numel(minidx) < cfg_gen.paramOpt.minCyclesNum)
					supThres{iTrl}(iCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2)) = 0;
				end
			end
			clear currSig
			tmpEvts = hvn_extrctBnryBouts(supThres{iTrl}(iCh,:)); % redo event search
		end
		
		% Discard events longer than specified duration + event bookkeeping
		numEvt  = 0; % initialize event counter
		tmpLen  = diff(tmpEvts,1,2)+1; % length of each event
		for kEvt = 1: size(tmpEvts,1)
			% If the event isnt too long, note it down together with all kinds of properties
			if tmpLen(kEvt) <= round(cfg_gen.criterion.duration(2)*fsample) && ... % inclusion criterion fullfilled w/o artifacts
					all(godfltr{iTrl}(currCh,tmpEvts(kEvt,1)+critPad(1) : tmpEvts(kEvt,2)+critPad(2)))
				
				numEvt = numEvt + 1;
				
				% General properties
				cfg_gen.evtIndiv(iCh,iTrl).staTime(numEvt)		= tmpEvts(kEvt,1);                   % event start (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).midTime(numEvt)		= round(mean(tmpEvts(kEvt,:),2));    % event cent (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).endTime(numEvt)		= tmpEvts(kEvt,2);                   % event end (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).duration(numEvt)		= tmpLen(kEvt) / fsample;  % event duration (in seconds)
				cfg_gen.evtIndiv(iCh,iTrl).stage(numEvt)		= inData.staging{iTrl}(tmpEvts(kEvt,1));
				
				% Properties of filtered signal
				tmpWin											= thres{iTrl}(currCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2));
				[minAmp,minIdx]									= min(tmpWin);
				[maxAmp,maxIdx]									= max(tmpWin);
				cfg_gen.evtIndiv(iCh,iTrl).maxTime(numEvt)		= tmpEvts(kEvt,1) + maxIdx - 1; % time of event peak (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).minTime(numEvt)		= tmpEvts(kEvt,1) + minIdx - 1; % time of event trough (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).minAmp(numEvt)		= minAmp;
				cfg_gen.evtIndiv(iCh,iTrl).maxAmp(numEvt)		= maxAmp;
				
				% Properties of envelope
				tmpWin											= thresEnv{iTrl}(currCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2));
				[maxAmp,maxIdx]									= max(tmpWin);
				cfg_gen.evtIndiv(iCh,iTrl).envMaxAmp(numEvt)    = maxAmp;                     % RMS max
				cfg_gen.evtIndiv(iCh,iTrl).envMaxTime(numEvt)   = tmpEvts(kEvt,1) + maxIdx - 1;  % time of RMS max (in datapoints)
				cfg_gen.evtIndiv(iCh,iTrl).envMean(numEvt)      = mean(tmpWin,2);
				cfg_gen.evtIndiv(iCh,iTrl).envSum(numEvt)       = sum(tmpWin,2);
				
				% Add peaks and troughs
				tmpWin          = thres{iTrl}(currCh,tmpEvts(kEvt,1):tmpEvts(kEvt,2));
				[~, peaks]      = findpeaks(tmpWin,tmpEvts(kEvt,1):tmpEvts(kEvt,2));
				[~, troughs]    = findpeaks((-1) * tmpWin,tmpEvts(kEvt,1):tmpEvts(kEvt,2));
				
				cfg_gen.evtIndiv(iCh,iTrl).peaks{numEvt}    = peaks;
				cfg_gen.evtIndiv(iCh,iTrl).troughs{numEvt}  = troughs;
				cfg_gen.evtIndiv(iCh,iTrl).freq(numEvt)     = fsample / mean([diff(peaks),diff(troughs)]);
			end
			cfg_gen.evtIndiv(iCh,iTrl).numEvt = numEvt;   % Save number of detected events
		end
	end
end
clear thres thresEnv supThres

%% Optional: False positive rejection based on frequency profile
if cfg_gen.doFalseposRjct
	fprintf('----- Event rejection by frequency profile\n');
	cfg_gen.falseposRjct.rejects = cell(numDetct,numTrl);
	
	for iTrl = 1 : numTrl
		for iCh = 1 : numDetct
			fprintf('Channel %d/%d (%s): ', iCh, numDetct,cfg_gen.chDetct{iCh});
			
			tmpTic = tic;
			currCh = ismember(inData.label,cfg_gen.chDetct{iCh});
			
			%--- Segment data
			tfg         = [];
			tfg.trl     = [cfg_gen.evtIndiv(iCh,iTrl).envMaxTime' + round(cfg_gen.falseposRjct.timePad(1) * fsample),...
				cfg_gen.evtIndiv(iCh,iTrl).envMaxTime' + round(cfg_gen.falseposRjct.timePad(2) * fsample),...
				ones(cfg_gen.evtIndiv(iCh,iTrl).numEvt,1) * round(cfg_gen.falseposRjct.timePad(1) * fsample)];
			tmpTrls   = ft_redefinetrial(tfg,inData);
			
			%--- Calculate time frequency representation
			tfg             = [];
			tfg.channel     = tmpTrls.label(currCh);
			tfg.taper       = 'hanning';
			tfg.method      = 'mtmconvol';
			tfg.pad         = 'nextpow2';
			tfg.output      = 'pow';
			tfg.keeptrials  = 'yes';
			tfg.foi         = cfg_gen.falseposRjct.tfrFreq;
			tfg.toi         = cfg_gen.falseposRjct.tfrTime;
			tfg.t_ftimwin   = cfg_gen.falseposRjct.tfrWin;
			
			tmpTFR  = ft_freqanalysis(tfg,tmpTrls);
			tmpPow  = squeeze(tmpTFR.powspctrm);                                        %% Note: rpt x freq x time
			tmpTime = arrayfun(@(x) nearest(tmpTFR.time,x),cfg_gen.falseposRjct.avgWin);
			tmpFreq = arrayfun(@(x) nearest(tmpTFR.freq,x),cfg_gen.falseposRjct.freqlim);
			
			%--- Perform event rejection
			cfg_gen.falseposRjct.rejects{iCh,iTrl} = ones(size(tmpPow,1),1,'logical');
			
			tmpPow = squeeze(sum(tmpPow(:,:,tmpTime(1):tmpTime(2)),3));
			tmpMax = max(tmpPow,[],2);                                      %% Determine maximum per trial
			tmpPow = tmpPow ./ repmat(tmpMax,1,size(tmpPow,2));             %% Normalise by maximum value
			
			for iEvt = 1 : size(tmpPow,1)
				[~, tmpPks,~,tmpProm] = findpeaks(tmpPow(iEvt,:));
				
				hazMax = find(tmpPks >= tmpFreq(1) & tmpPks <= tmpFreq(2));
				if numel(hazMax) > 0 && any(tmpProm(hazMax) > cfg_gen.falseposRjct.prominence)
					cfg_gen.falseposRjct.rejects{iCh,iTrl}(iEvt) = 0;
				end
			end
			
			fprintf(' reject %d of %d (%.2f) - took %.2f s\n',...
				sum(cfg_gen.falseposRjct.rejects{iCh,iTrl}),...
				size(tmpPow,1),...
				100 * sum(cfg_gen.falseposRjct.rejects{iCh,iTrl}) / size(tmpPow,1),...
				toc(tmpTic));
			
			clear tmpTrls tmpTFR tmpPow
			
		end
	end
end


%% add summary statistics to cfg.summary
for iTrl = 1 : numTrl % loop over trials
	for iCh = 1 : numDetct
		cfg_gen.summary(iCh,iTrl).label = cfg_gen.chDetct{iCh};
		cfg_gen.summary(iCh,iTrl).tss   = cfg_gen.evtIndiv(iCh,iTrl).tss;
		
		if cfg_gen.doFalseposRjct
			cfg_gen.summary(iCh,iTrl).numEvt    = sum(~cfg_gen.falseposRjct.rejects{iCh,iTrl});
			tmpIdx                          = ~cfg_gen.falseposRjct.rejects{iCh,iTrl};
		else
			cfg_gen.summary(iCh,iTrl).numEvt    = cfg_gen.evtIndiv(iCh,iTrl).numEvt;
			tmpIdx                          = ones(cfg_gen.summary(iCh,iTrl).numEvt,1,'logical');
		end
		
		if cfg_gen.summary(iCh,iTrl).numEvt > 0
			cfg_gen.summary(iCh,iTrl).density   = cfg_gen.summary(iCh,iTrl).numEvt / cfg_gen.summary(iCh,iTrl).tss;
			cfg_gen.summary(iCh,iTrl).duration  = mean(cfg_gen.evtIndiv(iCh,iTrl).duration(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).freq      = mean(cfg_gen.evtIndiv(iCh,iTrl).freq(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).minAmp    = mean(cfg_gen.evtIndiv(iCh,iTrl).minAmp(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).maxAmp    = mean(cfg_gen.evtIndiv(iCh,iTrl).maxAmp(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).envMax    = mean(cfg_gen.evtIndiv(iCh,iTrl).envMaxAmp(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).envMean   = mean(cfg_gen.evtIndiv(iCh,iTrl).envMean(tmpIdx),2);
			cfg_gen.summary(iCh,iTrl).envSum    = mean(cfg_gen.evtIndiv(iCh,iTrl).envSum(tmpIdx),2);
		end
	end
end
end

function outBnry = hvn_createBnrySignal(inBnry, inLen)
% Created by H.-V.V. Ngo
%
% Creates a logical vector with '1'-bouts given by inBinary
%
% Usage: outBnry = hvn_createBnrySignal(inBnry, inLen)
%
% Input parameters
% inBnry    = [N x 2] array contraining N 'active'-bouts, with their on-
%             and offset depicted by the 1st and 2nd column, respectively
% inLen     = Length of resulting vector in data points
outBnry = zeros(1,inLen,'logical');
if ~isempty(inBnry)
	binaryIdx = cell2mat(reshape(arrayfun(@(x,y)  x:y, inBnry(:,1), inBnry(:,2),'UniformOutput',0),1,size(inBnry,1)));
	outBnry(binaryIdx) = 1;
end
end

function outBnry = hvn_extrctBnryBouts(inBnry)
% Created by H.-V.V. Ngo
%
% Extracts on- and offset of 'active'-bouts within a logical vector
%
% Usage: outBnry = hvn_extrctBnryBouts(inBnry)
%
% Input parameters
% inBnry    = 1-D logical vector with active bouts (periods with value 1)
%
% Output
% outBnry   = [N x 2] vector of N bouts with on- and offset represented by
%             1st and 2nd column
if size(inBnry,1) ~= 1
	inBnry = inBnry';
end
bnryOn      = find(diff([0 inBnry]) ==  1);
bnryOff     = find(diff([inBnry 0]) == -1);
numBouts    = min(numel(bnryOn),numel(bnryOff));
outBnry = [bnryOn(1:numBouts)', bnryOff(1:numBouts)'];
end