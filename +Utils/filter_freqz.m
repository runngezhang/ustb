function varargout = filter_freqz(varargin)
% FILTER_FREQZ M-file for filter_freqz.fig
%      FILTER_FREQZ, by itself, creates a new FILTER_FREQZ or raises the existing
%      singleton*.
%
%      H = FILTER_FREQZ returns the handle to a new FILTER_FREQZ or the handle to
%      the existing singleton*.
%
%      FILTER_FREQZ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FILTER_FREQZ.M with the given input arguments.
%
%      FILTER_FREQZ('Property','Value',...) creates a new FILTER_FREQZ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before filter_freqz_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to filter_freqz_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help filter_freqz

% Last Modified by GUIDE v2.5 22-Sep-2010 13:35:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @filter_freqz_OpeningFcn, ...
                   'gui_OutputFcn',  @filter_freqz_OutputFcn, ...
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


% --- Executes just before filter_freqz is made visible.
function filter_freqz_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to filter_freqz (see VARARGIN)

% Choose default command line output for filter_freqz
handles.output = hObject;

handles.choosen_indxs = [];

handles.filters = read_filter_file(varargin{1});

filter_str = {};
filter_indx = [];
cnt = 1;
for kk=1:length(handles.filters)
    if ~isempty(handles.filters{kk})
        fprintf('ID %d \n',handles.filters{kk}.id)
        filter_str{cnt} = sprintf('%s (%d)',handles.filters{kk}.type,handles.filters{kk}.num_taps);
        filter_indx = [filter_indx kk];
        cnt = cnt + 1;
    end
end

handles.filter_indx = filter_indx;
handles.filter_str = filter_str;
handles.fs = 40e6;
handles.Nfft = str2double(get(handles.nfft_edit,'String'));

df = handles.fs/handles.Nfft;
handles.f_ax = ifftshift((0:handles.Nfft-1)'-floor(handles.Nfft/2))*df;

handles.xlim = str2num(get(handles.xlim_edit,'String'));
handles.ylim = str2num(get(handles.ylim_edit,'String'));

set(handles.filters_list,'String',filter_str);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes filter_freqz wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = filter_freqz_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in filters_list.
function filters_list_Callback(hObject, eventdata, handles)
% hObject    handle to filters_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns filters_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filters_list
indx = get(hObject,'Value');
filt = handles.filters{handles.filter_indx(indx)};
set(handles.id_label,'String',sprintf('ID : %d',filt.id));
set(handles.type_label,'String',sprintf('%s',filt.type));
set(handles.ntaps_label,'String',sprintf('Ntaps : %d',filt.num_taps));
set(handles.ntaps_label,'String',sprintf('Desc : %s',filt.desc));

if get(handles.hold_on_check,'Value')
    handles.choosen_indxs = unique([handles.choosen_indxs handles.filter_indx(indx)]);
else
    handles.choosen_indxs = handles.filter_indx(indx);
end

handles = plot_data(handles);

guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function filters_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filters_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hold_on_check.
function hold_on_check_Callback(hObject, eventdata, handles)
% hObject    handle to hold_on_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold_on_check



function ylim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ylim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylim_edit as text
%        str2double(get(hObject,'String')) returns contents of ylim_edit as a double

value = handles.ylim;
if ~isempty(str2num(get(hObject,'String')))
    value = str2num(get(hObject,'String'));
end
handles.ylim = value;

handles = plot_data(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ylim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to xlim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlim_edit as text
%        str2double(get(hObject,'String')) returns contents of xlim_edit as a double
value = handles.xlim;
if ~isempty(str2num(get(hObject,'String')))
    value = str2num(get(hObject,'String'));
end
handles.xlim = value;

handles = plot_data(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function xlim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = plot_data(handles)
    
Nfft = handles.Nfft;
cmap = lines;
cla(handles.plot_axes);
for kk=1:length(handles.choosen_indxs)
    filt = handles.filters{handles.choosen_indxs(kk)};
    H = fft(filt.taps,Nfft);
    if get(handles.hold_on_check,'Value')
        hold(handles.plot_axes,'on');
    end
    ph = plot(handles.plot_axes,1e-6*handles.f_ax(1:Nfft/2),20*log10(abs(H(1:Nfft/2)/max(H(1:Nfft/2)))));
    if get(handles.hold_on_check,'Value')
        hold(handles.plot_axes,'off');
    end
    set(ph,'LineWidth',2);
    set(ph,'Color',cmap(kk,:));
end
ylim(handles.ylim),xlim(handles.xlim)
xlabel('[MHz]'),ylabel('[dB]')
grid on;



function nfft_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nfft_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nfft_edit as text
%        str2double(get(hObject,'String')) returns contents of nfft_edit as a double

value = handles.Nfft;
if ~isempty(str2double(get(hObject,'String')))
    value = str2double(get(hObject,'String'));
end
handles.Nfft = value;

df = handles.fs/handles.Nfft;
handles.f_ax = ifftshift((0:handles.Nfft-1)'-floor(handles.Nfft/2))*df;

handles = plot_data(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function nfft_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nfft_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
