function varargout = trial_selector(varargin)
% GUI allowing to manually exclude trials from a fieldtrip structure. Shows
% 20 trials at once. Input data is expected to be a standard fieldtrip
% dataset with one channel and the fields 'trials' and 'time'.
%
% Use as:
%  output = trial_selector(cfg)
%  output = trial_selector(cfg, data)
%  trial_selector(cfg, data)
%
% INPUT VARIABLES:
% cfg               structure; configuration with fields:
% .inputfile		string (optional); path to dataset
%					alternatively, data can be provided as a second input
%					argument
% .page				int; on which page to start on (default: 1)
% .latency			int array; time limits for displaying trials 
%					(e.g. [-1 4]
% .comment			cell array of strings with as many entries as there are
%					trials; will be shown in the corner of each plot; could
%					for instance be used to show trial conditions
%
% OUTPUT VARIABLES:
% output			array of logicals; is only returned if cfg.outputfile
%					is not specified or empty
%
% COMMENTS:
% All important stuff happens in:
% trial_selector_OpeningFcn		runs once at the beginning 
% initializePage				is called for displaying each new page
%
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trial_selector_OpeningFcn, ...
                   'gui_OutputFcn',  @trial_selector_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trial_selector is made visible.
function trial_selector_OpeningFcn(hObject, eventdata, handles, varargin)

% For debugging
% if numel(varargin) == 0
% 	varargin{1} = load_file('debug_cfg.mat');
% 	varargin{2} = varargin{1}.data;
% 	varargin{1} = rmfield(varargin{1}, 'data');
% end
% hObject.Visible = 'on';

% Make all fields in input available
if numel(varargin) > 0 && isstruct(varargin{1})
	f = fields(varargin{1});
	for iField = 1:numel(f)
		handles.(f{iField}) = varargin{1}.(f{iField});
	end
else
	error('No configuration given as input.')
end

% Check input and load data
if isfield(handles, 'inputfile') && isempty(handles.inputfile)
	handles = rmfield(handles, 'inputfile');
end
if numel(varargin) == 2
	if isfield(handles, 'inputfile')
		error('You cannot provide both cfg.inputfile and data as input.')
	else
		handles.data = varargin{2};
	end
else
	temp = load(handles.inputfile);
	names = fieldnames(temp);
	if length(names) ~= 1
		error('Unexpected content in preprocessed data file. File should contain one single structure.')
	end
	handles.data = temp.(names{1});
	clear temp names
end

% Check if data is the way we expect it
requiredFields = {'trial', 'time'};
for i = requiredFields
	if ~isfield(handles.data,i)
		error(['Required field missing in data: ' i{1} '.']);
	end
end
if size(handles.data.trial, 1) ~= 1
	error('Data should have only one channel.')
end

% Set defaults
if ~isfield(handles, 'page'), handles.page = 1; end
if ~isfield(handles, 'ylim'), handles.ylim = [-150 150]; end
if ~isfield(handles, 'xlim'), handles.xlim = [-2 6]; end
% if ~isfield(handles, 'outputfile'), handles.outputfile = []; end

% Bookkeeping
handles.ax_handles	= {handles.axes1, handles.axes2, handles.axes3, handles.axes4, handles.axes5, handles.axes6, handles.axes7, handles.axes8, handles.axes9, handles.axes10, handles.axes11, handles.axes12, handles.axes13, handles.axes14, handles.axes15, handles.axes16, handles.axes17, handles.axes18, handles.axes19, handles.axes20}; % TODO: Add all 
handles.num_axes	= numel(handles.ax_handles);
handles.num_trials	= size(handles.data.trial,2);
handles.num_pages	= ceil(handles.num_trials / handles.num_axes);
handles.output		= false(handles.num_trials,1); % eventual output

% Bring page selection in order
set(handles.pagetext2, 'String', [' / ' num2str(handles.num_pages)]);

% Make sure all shenanigans are save in the almighty cloud
guidata(hObject, handles);

% Trigger the first page (mostly plotting)
initializePage(hObject, handles)

% Output voodoo in case we'll return result to user
if ~isfield(handles, 'outputfile'), uiwait(handles.figure1); end
% uiwait(handles.figure1);

% --- Shows a new page (potting, bookkeeping etc.)
function initializePage(hObject, handles)
% Determine what to plot on the current page
handles.ax_trial		= zeros(20,1); 
first_trial				= (handles.page-1) * handles.num_axes + 1;
last_trial				= handles.page * handles.num_axes;
if last_trial > handles.num_trials, last_trial = handles.num_trials; end
cnt = 1;
for iTrial = first_trial:last_trial
	handles.ax_trial(cnt)	= iTrial;
	cnt						= cnt + 1;
end

% Update buttons
set(handles.pagesel, 'String', num2str(handles.page));
if handles.page == 1
	handles.previous.Enable		= 'off';
	handles.next.Enable			= 'on';
elseif handles.page == handles.num_pages
	handles.next.Enable			= 'off';
	handles.previous.Enable		= 'on';
else
	handles.next.Enable			= 'on';
	handles.previous.Enable		= 'on';
end

% Do the actual plotting
for iAx = 1:handles.num_axes
	ax		= handles.ax_handles{iAx}; % current axis
	tr		= handles.ax_trial(iAx); % which trial to plot
	cla(ax)
	if tr ~= 0
		% Do line plot
		[~,l1]	= min(abs(handles.data.time{tr} - handles.xlim(1))); % lower time limit
		[~,l2]	= min(abs(handles.data.time{tr} - handles.xlim(2))); % upper time limit
		h				= plot(ax, handles.data.time{tr}(l1:l2), handles.data.trial{tr}(l1:l2));
		h.LineWidth		= 1.0;
		h.Color			= 'k';
		hold(ax, 'on')
		h2				= plot(ax, [0 0], handles.ylim); % vertical line at 0
		h2.LineWidth	= .5;
		h2.Color		= 'r';
		if isfield(handles, 'comment') && ~isempty(handles.comment)
			current_annotag = ['anno' num2str(iAx)];
			delete(findall(gcf,'Tag',current_annotag)); % otherwise they stack up and look ugly
			annotation('textbox',ax.Position,'String', handles.comment{iAx},'FitBoxToText','on', 'EdgeColor','none', 'Tag', current_annotag); % + [.001 -.001 0 0]
		end
		
		% Change some axes properties
		ax.YLim			= handles.ylim;
		ax.YTickLabel	= '';
		if handles.output(handles.ax_trial(iAx))
			ax.Color			= [1 .7 .7];
		end		
		if iAx > numel(handles.ax_trial) - 5 && iAx <= numel(handles.ax_trial)
			ax.XTick		= handles.data.time{tr}(l1)-rem(handles.data.time{tr}(l1),2):2:handles.data.time{tr}(l2)-rem(handles.data.time{tr}(l2),2); % :P
		else
			ax.XTick		= [];
		end
		hold(ax, 'off')
		
		% Behavior when clicked
		set(h,'HitTest','off') % otherwise those plots will register a button down
		set(h2,'HitTest','off') %... and not our axes
		set(ax,'ButtonDownFcn',{@toggleAxes, iAx}) % for the axes, please call our own button-down function
	else
		cla(ax) % if there ain't no data
	end
end
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = trial_selector_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure. We need to check
% whether outputfile is provided because a) nargout is not available within
% a GUI for some reason (nargout is always 0) and b) because this output
% function is called right away if uiwait is not enabled (absolute
% nonsense, thanks matlab). 
% If you know of a better solution to the whole output disaster, please let
% me know!
if ~isfield(handles, 'outputfile')
	varargout{1} = handles.output;
	delete(hObject);
end


% --- Selects/unselects plots.
function toggleAxes(hObject, ~, iAx)
handles = guidata(hObject);
if ~handles.output(handles.ax_trial(iAx))
	handles.ax_handles{iAx}.Color			= [1 .7 .7];
	handles.output(handles.ax_trial(iAx))	= true;
else
	handles.ax_handles{iAx}.Color			= [1 1 1];
	handles.output(handles.ax_trial(iAx))	= false;
end
guidata(hObject, handles);


% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.page > 1
	handles.page = handles.page-1;
	guidata(hObject, handles);
	initializePage(hObject, handles);
end


% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.page < handles.num_pages
	handles.page = handles.page+1;
	guidata(hObject, handles);
	initializePage(hObject, handles);
end


function pagesel_Callback(hObject, ~, handles)
% hObject    handle to pagesel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pagesel as text
%        str2double(get(hObject,'String')) returns contents of pagesel as a double
in = str2double(get(hObject,'String'));
if in < 1 || in > handles.num_pages
	set(hObject,'String', num2str(handles.page));
	set(handles.pageinvalid, 'Visible', 'on');
	initializePage(hObject, handles);
else
	handles.page = in;
	guidata(hObject, handles);
	set(handles.pageinvalid, 'Visible', 'off');
	initializePage(hObject, handles);
end


% --- Executes during object creation, after setting all properties.
function pagesel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pagesel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 close(gcf); % this leads to a call to figure1_CloseRequestFcn, which closereq doesn't


% --- Executes when user attempts to close figure1.
% Only executed if X is pressed.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This is taken from the internet. If we resume here, the output function
% will still have access to the handles, otherwise its empty. Not really
% sure why otherwise we delete the object, but that code is never reached
% anyways..
if isequal(get(hObject,'waitstatus'),'waiting')
	uiresume(hObject);
    guidata(hObject,handles);
else
 	if ~isfield(handles, 'outputfile')|| isempty(handles.outputfile)
		warning('Something went wrong, no outputfile provided.')
	else
		out = handles.output;
		disp(['Saving output to ' handles.outputfile])
		save(handles.outputfile, 'out');
	end
    delete(hObject); % closes the figure
	
end
