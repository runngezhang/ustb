function varargout = make_maps(varargin)
% MAKE_MAPS M-file for make_maps.fig
%      MAKE_MAPS, by itself, creates a new MAKE_MAPS or raises the existing
%      singleton*.
%
%      H = MAKE_MAPS returns the handle to a new MAKE_MAPS or the handle to
%      the existing singleton*.
%
%      MAKE_MAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAKE_MAPS.M with the given input arguments.
%
%      MAKE_MAPS('Property','Value',...) creates a new MAKE_MAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before make_maps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to make_maps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help make_maps

% Last Modified by GUIDE v2.5 31-Mar-2009 14:43:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @make_maps_OpeningFcn, ...
                   'gui_OutputFcn',  @make_maps_OutputFcn, ...
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


% --- Executes just before make_maps is made visible.
function make_maps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to make_maps (see VARARGIN)

% Choose default command line output for make_maps
handles.output = hObject;

handles.image.gain_val = get(handles.gain_slider,'Value');
handles.envmap.gamma_val = get(handles.gamma_env_slider,'Value');

handles.imgmap.gamma_val = get(handles.gamma_map_slider,'Value');
handles.imgmap.contrast_val = get(handles.contrast_slider,'Value');
handles.imgmap.brightness_val = get(handles.brightness_slider,'Value');

if isempty(varargin)
    handles.image.iq = [];
else
    handles.image.iq = varargin{1};
end

handles = Update(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes make_maps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = make_maps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.envmap.map;
varargout{2} = handles.imgmap.map;


% --- Executes on slider movement.
function gamma_env_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_env_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.envmap.gamma_val = get(hObject,'Value');
handles = Update(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gamma_env_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_env_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function gain_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.image.gain_val = get(hObject,'Value');
handles = Update(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gain_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function gamma_map_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_map_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.imgmap.gamma_val = get(hObject,'Value');
handles = Update(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function gamma_map_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_map_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function contrast_slider_Callback(hObject, eventdata, handles)
% hObject    handle to contrast_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.imgmap.contrast_val = get(hObject,'Value');
handles = Update(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function contrast_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to contrast_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function brightness_slider_Callback(hObject, eventdata, handles)
% hObject    handle to brightness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.imgmap.brightness_val = get(hObject,'Value');
handles = Update(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function brightness_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brightness_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function handles = Update(handles)

%
max_val16 = 2^16 - 1;
max_val8 = 255;

envmap = floor(max_val8*((0:max_val16)/max_val16).^handles.envmap.gamma_val);

contrast = floor(handles.imgmap.contrast_val);
if mod(contrast,2) ~= 0
    contrast = contrast + 1;
end
brightness = handles.imgmap.brightness_val;
max_val = max_val8 - contrast;
gamma2 = handles.imgmap.gamma_val;

imgmap = [zeros(1,contrast/2) floor(max_val8*((0:max_val)/max_val).^gamma2) max_val8*ones(1,contrast/2)]+ brightness;
imgmap(imgmap > length(imgmap)) = length(imgmap);

% Image
if ~isempty(handles.image.iq)
    gain_dB = handles.image.gain_val;
    gain = 10^(gain_dB/20);
    
    env = ceil(gain*abs(handles.image.iq));
    env(env > max_val16) = max_val16;
    
    img1 = envmap(env);            
    img2 = imgmap(img1+1);
    
    axes(handles.img_axes);
    image(img2);
    colormap(gray(256));
else
    img1 = [];
    img2 = [];
end

% EnvMap
[H,Nax] = hist(env(:),max_val16/10);
H = H/max(H(:));

axes(handles.env_compression_axes)

[Ax,H1,H2] = plotyy(0:max_val16,envmap,Nax,H,'plot','plot');

ylim(Ax(1),[0 max_val8])
xlim(Ax(1),[0 max_val16])
xlim(Ax(2),[0 max_val16])
grid(Ax(1),'on');
grid(Ax(2),'on');

% ImgMap
[H,Nax] = hist(img1(:),max_val8/4);
H = H/max(H(:));

axes(handles.map_axes)
[Ax] = plotyy(0:max_val8,imgmap,Nax,H);
ylim(Ax(1),[0 max_val8])

xlim(Ax(1),[0 max_val8])
xlim(Ax(2),[0 max_val8])

grid(Ax(1),'on');
grid(Ax(2),'on');

handles.imgmap.map = imgmap;
handles.envmap.map = envmap;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','imgmap',handles.imgmap.map);
assignin('base','envmap',handles.envmap.map);
assignin('base','gain_dB',handles.image.gain_val);

