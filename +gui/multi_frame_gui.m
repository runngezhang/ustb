function varargout = multi_frame_gui(varargin)
% MULTI_FRAME_GUI MATLAB code for multi_frame_gui.fig
%
% A simple GUI to display ultrasound images with multiple frames in the 
% USTB. The GUI allows you to click through all frames using a "Previous
% Frame" and "Next Frame" button. It also allows to play the loop of all
% frames as a movie. 
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
%   Authors: Ole Marius Hoel Rindal (olemarius@olemarius.net)
%   $Date : 2017/03/20$



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @multi_frame_gui_OpeningFcn, ...
    'gui_OutputFcn',  @multi_frame_gui_OutputFcn, ...
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

global b_data;
global current_frame;
global image_handle;
global all_images_dB;
global play_loop;
global in_title;
play_loop = false;

% --- Executes just before multi_frame_gui is made visible.
function multi_frame_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to multi_frame_gui (see VARARGIN)

% Choose default command line output for multi_frame_gui
handles.output = hObject;

assert(strcmp(class(varargin{1}),'huff.beamformed_data'),'First input should be a beamformed_data object.')

global b_data;
global current_frame;
global image_handle;
global all_images_dB;
global in_title;
b_data = varargin{1};
in_title = varargin{2};
dynamic_range = varargin{3};
current_frame = 1;
[image_handle, all_images_dB] = b_data.draw_image(handles.image_handle,[in_title,', Frame = ',num2str(current_frame),'/',num2str(size(all_images_dB,3))],dynamic_range);

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes multi_frame_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = multi_frame_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonPrevious.
function pushbuttonPrevious_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPrevious (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global current_frame;
global image_handle;
global all_images_dB;
global in_title;
if current_frame-1 > 0
    current_frame = current_frame-1;
    set(image_handle,'CData',all_images_dB(:,:,current_frame));
    title([in_title,', Frame = ',num2str(current_frame),'/',num2str(size(all_images_dB,3))]);
end

% --- Executes on button press in pushbuttonNext.
function pushbuttonNext_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global current_frame;
global image_handle;
global all_images_dB;
global in_title;
if current_frame+1 <= size(all_images_dB,3)
    current_frame = current_frame+1;
    set(image_handle,'CData',all_images_dB(:,:,current_frame));
    title(handles.image_handle,[in_title,', Frame = ',num2str(current_frame),'/',num2str(size(all_images_dB,3))]);
end


% --- Executes on button press in pushbuttonPlayLoop.
function pushbuttonPlayLoop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonPlayLoop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global play_loop;
global current_frame;
global image_handle;
global all_images_dB;
global in_title;
play_loop = ~play_loop;
%current_idx = current_frame;
while(play_loop)
    current_frame = current_frame+1;
    if current_frame > size(all_images_dB,3)
        current_frame = 1;
    end
    set(image_handle,'CData',all_images_dB(:,:,current_frame));
    title(handles.image_handle,[in_title,', Frame = ',num2str(current_frame),'/',num2str(size(all_images_dB,3))]);
    drawnow();
    pause(0.05);
end