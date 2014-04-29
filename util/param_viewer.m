function varargout = param_viewer(varargin)
% PARAM_VIEWER M-file for param_viewer.fig
%      PARAM_VIEWER, by itself, creates a new PARAM_VIEWER or raises the existing
%      singleton*.
%
%      H = PARAM_VIEWER returns the handle to a new PARAM_VIEWER or the handle to
%      the existing singleton*.
%
%      PARAM_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARAM_VIEWER.M with the given input arguments.
%
%      PARAM_VIEWER('Property','Value',...) creates a new PARAM_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before param_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to param_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help param_viewer

% Last Modified by GUIDE v2.5 26-Aug-2009 10:04:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @param_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @param_viewer_OutputFcn, ...
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


% --- Executes just before param_viewer is made visible.
function param_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to param_viewer (see VARARGIN)

% Choose default command line output for param_viewer
handles.output = hObject;

if length(varargin) > 0
    params = varargin{1};
    param_names = [];
    params_tmp = [];
    cnt = 1;
    if length(varargin) == 2
        paramList = varargin{2};
    else
        paramList = 1:length(params);
    end;
    for k=1:length(params)
        if ~isempty(params{k}) && ismember(params{k}.id,paramList)
            if ~isstruct(params{k}.value) && (numel(params{k}.value)==1)
                unit = '';
                if ~strcmp(params{k}.unit,'N/A')
                    unit = params{k}.unit;
                end;
                param_names{cnt} = sprintf('%s   : %i %s',params{k}.name,params{k}.value,unit);
            else
                param_names{cnt} = params{k}.name;
            end;
            params_tmp{cnt} = params{k};
            cnt = cnt + 1;
        end
    end
    
    [param_names indx] = sort(param_names);
   
    set(handles.param_list,'String',param_names);
    
    handles.params = params_tmp;
    handles.sorted_indx = indx;
    
    set(handles.param_list,'Value',1);
    if length(varargin) > 1
        set(hObject,'Name', ['Parameter viewer: ',varargin{2}]);
    end;
else
    handles.params = [];
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes param_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = param_viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in param_list.
function param_list_Callback(hObject, eventdata, handles)
% hObject    handle to param_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns param_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from param_list

if isempty(handles.params)
    return;
end
    
param = handles.params{handles.sorted_indx(get(hObject,'Value'))};
set(handles.id_label,'String',sprintf('ID : %d',param.id));
set(handles.unit_label,'String',sprintf('Unit : %s',param.unit));

labels = 1:8;
for k=1:length(labels)
    set(handles.(sprintf('value_label%d',labels(k))),'String','-');
end
switch param.typeid
    case 0
        % Int        
        set(handles.value_label1,'String',sprintf('%d',param.value));
    case 2
        % Rect       
        set(handles.value_label1,'String',sprintf('l : %d',param.value.left));
        set(handles.value_label2,'String',sprintf('r : %d',param.value.right));
        set(handles.value_label3,'String',sprintf('t : %d',param.value.top));
        set(handles.value_label4,'String',sprintf('b : %d',param.value.bottom));
    case 6
        % Curve        
        set(handles.value_label1,'String',sprintf('t : %d',param.value.t));
        set(handles.value_label2,'String',sprintf('m : %d',param.value.m));
        set(handles.value_label3,'String',sprintf('b : %d',param.value.b));
        set(handles.value_label4,'String',sprintf('vm : %d',param.value.vm));
    case 8
        % 
        for k=1:length(labels)
            set(handles.(sprintf('value_label%d',labels(k))),'String',sprintf('%d',param.value(k)));
        end
    otherwise
        param.typeid
end

% --- Executes during object creation, after setting all properties.
function param_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to param_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
