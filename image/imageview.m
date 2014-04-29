function varargout = imageview(varargin)
% IMAGEVIEW M-file for imageview.fig
%      IMAGEVIEW, by itself, creates a new IMAGEVIEW or raises the existing
%      singleton*.
%
%      H = IMAGEVIEW returns the handle to a new IMAGEVIEW or the handle to
%      the existing singleton*.
%
%      IMAGEVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEVIEW.M with the given input arguments.
%
%      IMAGEVIEW('Property','Value',...) creates a new IMAGEVIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageview_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imageview_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageview

% Last Modified by GUIDE v2.5 02-Oct-2009 08:43:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageview_OpeningFcn, ...
                   'gui_OutputFcn',  @imageview_OutputFcn, ...
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


% --- Executes just before imageview is made visible.
function imageview_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imageview (see VARARGIN)
%hObject = clf(hObject,'reset');

set(hObject,'Name','Image Viewer');

width_one = 568;
height = 675 + 250;

width_total = width_one;
num_img = 1;
if length(varargin) == 3 || length(varargin) == 5
    width_total = width_one;
    num_img = 1;
elseif length(varargin) == 6 || length(varargin) == 8
    width_total = 2*width_one;
    num_img = 2;
end
set(hObject,'units','pixels');
pos = get(hObject,'outerposition');
set(hObject,'outerposition',[pos(1) pos(2) width_total height]);


for k=1:num_img
    handles.data{k}.iq = varargin{(k-1)*3 + 1};
    handles.data{k}.env2 = abs(handles.data{k}.iq).^2;
    handles.data{k}.gain = varargin{(k-1)*3 + 2};
    handles.data{k}.dyn = varargin{(k-1)*3 + 3};    
    
    if isempty(handles.data{k}.gain)
        handles.data{k}.gain = -20;
    end
    
    if isempty(handles.data{k}.dyn)
        handles.data{k}.dyn = 60;
    end
    
    num_frames(k) = size(handles.data{k}.iq,3);
end

if length(varargin{end-1}) > 1
    % Axis structure
    handles.x_ax = varargin{end-1};    
else
    handles.x_ax = [0 size(handles.data{1}.iq,2)-1];    
end

if length(varargin{end}) > 1
    % Axis structure
    handles.y_ax = varargin{end};    
else
    handles.y_ax = [0 size(handles.data{1}.iq,1)-1];    
end

for k=1:num_img        
    % Configure axes
    if ~isfield(handles,'axes_handles')
        axes_handle = axes;
        axes_handle_hist = axes;
    else
        if length(handles.axes_handles) < k
            axes_handle = axes;
            axes_handle_hist = axes;
        else
            axes_handle = handles.axes_handles(k);
            axes_handle_hist = handles.axes_handles_hist(k);
        end
    end
        
    set(axes_handle,'parent',hObject);
    set(axes_handle,'units','pixels');
    
    set(axes_handle_hist,'parent',hObject);
    set(axes_handle_hist,'units','pixels');
    
    if k > 1
        axes_pos = handles.axes_pos(:,1);
    else
        axes_pos = get(axes_handle,'position');    
        axes_pos(4) = 300;    
    end
    
    axes_pos(1) = 73 + (k-1)*width_total/2;
    axes_pos(2) = 75 + 350;
    axes_pos(3) = 450;
    set(axes_handle,'position',axes_pos);

    axes_pos_hist = axes_pos;
    axes_pos_hist(1) = 73 + (k-1)*width_total/2;
    axes_pos_hist(2) = 75;
    axes_pos_hist(3) = 450;
    set(axes_handle_hist,'position',axes_pos_hist);
    
    if ~isfield(handles,'gain_slider') || length(handles.gain_slider) < k        
        handles.gain_slider(k) = uicontrol(hObject, 'Style', 'slider', 'Tag', ['gain_slider' num2str(k)]);    
        handles.gain_label(k) = uicontrol(hObject, 'Style', 'text', 'Tag', ['gain_label' num2str(k)]);
        
        handles.dyn_slider(k) = uicontrol(hObject, 'Style', 'slider', 'Tag', ['dyn_slider' num2str(k)]);
        handles.dyn_label(k) = uicontrol(hObject, 'Style', 'text', 'Tag', ['dyn_label' num2str(k)]);
        
        if k == 1
            handles.frame_slider(k) = uicontrol(hObject, 'Style', 'slider', 'Tag', ['frame_slider' num2str(k)]);
            handles.frame_label(k) = uicontrol(hObject, 'Style', 'text', 'Tag', ['frame_label' num2str(k)]);
        end
    end
    
    set(handles.gain_slider(k),'units','pixels');
    gain_slider_pos = get(handles.gain_slider(k),'position');
    gain_slider_pos(1) = axes_pos(1);
    gain_slider_pos(2) = axes_pos(2) + axes_pos(4) + 70;
    gain_slider_pos(3) = axes_pos(3);
    set(handles.gain_slider(k),'Position',gain_slider_pos);
    set(handles.gain_slider(k),'Max',20);
    set(handles.gain_slider(k),'Min',-60);
    set(handles.gain_slider(k),'Value',handles.data{k}.gain);
    set(handles.gain_slider(k),'SliderStep',[1 5]./(20+60));
    
    fhandle = @(hObject,eventdata)imageview(['gain_slider' num2str(k) '_Callback'],hObject,eventdata,guidata(hObject));
    set(handles.gain_slider(k),'Callback',fhandle);
    
    set(handles.dyn_slider(k),'units','pixels');
    dyn_slider_pos = get(handles.dyn_slider(k),'position');
    dyn_slider_pos(1) = axes_pos(1);
    dyn_slider_pos(2) = axes_pos(2) + axes_pos(4) + 20;
    dyn_slider_pos(3) = axes_pos(3);
    set(handles.dyn_slider(k),'Position',dyn_slider_pos);
    set(handles.dyn_slider(k),'Max',90);
    set(handles.dyn_slider(k),'Min',15);
    set(handles.dyn_slider(k),'Value',handles.data{k}.dyn);
    set(handles.dyn_slider(k),'SliderStep',[1 5]./(90-15));
    
    fhandle = @(hObject,eventdata)imageview(['dyn_slider' num2str(k) '_Callback'],hObject,eventdata,guidata(hObject));
    set(handles.dyn_slider(k),'Callback',fhandle);
    
    if k == 1
        set(handles.frame_slider(k),'units','pixels');
        frame_slider_pos = get(handles.frame_slider(k),'position');
        frame_slider_pos(1) = axes_pos(1);
        frame_slider_pos(2) = axes_pos(2) + axes_pos(4) + 120;
        frame_slider_pos(3) = axes_pos(3);
        set(handles.frame_slider(k),'Position',frame_slider_pos);

        if num_frames(k) > 1        
            set(handles.frame_slider(k),'Max',num_frames(k));
            set(handles.frame_slider(k),'Min',1);
            set(handles.frame_slider(k),'Value',1);
            set(handles.frame_slider(k),'SliderStep',[1 5]./(num_frames(k) - 1));    
            fhandle = @(hObject,eventdata)imageview(['frame_slider' num2str(k) '_Callback'],hObject,eventdata,guidata(hObject));
            set(handles.frame_slider(k),'Callback',fhandle);
        else
            set(handles.frame_slider(k),'Enable','off');
        end
    end
    
    set(handles.gain_label(k),'units','pixels');
    gain_label_pos = get(handles.gain_label(k),'position');
    gain_label_pos(1) = axes_pos(1);
    gain_label_pos(2) = axes_pos(2) + axes_pos(4) + 90;
    gain_label_pos(3) = axes_pos(3);
    set(handles.gain_label(k),'Position',gain_label_pos);
    set(handles.gain_label(k),'String',sprintf('Gain : %2.1f dB',handles.data{k}.gain));
    set(handles.gain_label(k),'FontSize',10);
    
    set(handles.dyn_label(k),'units','pixels');
    dyn_label_pos = get(handles.dyn_label(k),'position');
    dyn_label_pos(1) = axes_pos(1);
    dyn_label_pos(2) = axes_pos(2) + axes_pos(4) + 40;
    dyn_label_pos(3) = axes_pos(3);
    set(handles.dyn_label(k),'Position',dyn_label_pos);
    set(handles.dyn_label(k),'String',sprintf('Dynamic Range : %2.1f dB',handles.data{k}.dyn));
    set(handles.dyn_label(k),'FontSize',10);
    
    if k == 1
        set(handles.frame_label(k),'units','pixels');
        frame_label_pos = get(handles.frame_label(k),'position');
        frame_label_pos(1) = axes_pos(1);
        frame_label_pos(2) = axes_pos(2) + axes_pos(4) + 140;
        frame_label_pos(3) = axes_pos(3);
        set(handles.frame_label(k),'Position',frame_label_pos);
        set(handles.frame_label(k),'String',sprintf('Frame # : %d',1));
        set(handles.frame_label(k),'FontSize',10);
        
        handles.done_button = uicontrol(hObject, 'Style', 'pushbutton', 'Tag', 'DoneButton');    
        set(handles.done_button,'units','pixels');
        button_pos = get(handles.done_button,'position');
                
        button_pos(1) = axes_pos_hist(1);
        button_pos(2) = axes_pos_hist(2) - 60;
        button_pos(3) = axes_pos_hist(3);
        button_pos(4) = 1.1*button_pos(4);
        
        set(handles.done_button,'Position',button_pos,'String','Done','FontSize',10);
        fhandle = @(hObject,eventdata)imageview('DoneButton_Callback',hObject,eventdata,guidata(hObject));
        set(handles.done_button,'Callback',fhandle);
    end
    
    handles.axes_handles(k) = axes_handle;
    handles.axes_handles_hist(k) = axes_handle_hist;
    handles.axes_pos(:,k) = axes_pos;        
    handles.axes_pos_hist(:,k) = axes_pos_hist;        
end

handles.frame_indx = 1;

handles.num_img = num_img;
for k=1:num_img
    handles = UpdateImages(handles,k,handles.frame_indx);
end

% Choose default command line output for imageview
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imageview wait for user response (see UIRESUME)
%uiwait(handles.figure1)


% --- Outputs from this function are returned to the command line.
function varargout = imageview_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%cla(handles.axes_handles);
% Get default command line output from handles structure
%waitfor(handles.figure1,'BeingDeleted');
uiwait(handles.figure1)
if isempty(handles)
    for k=1:nargout
        varargout{k} = 0;        
    end
    return;
end

for k=1:handles.num_img
    varargout{(k-1)*2 + 1} = handles.data{k}.gain;
    varargout{(k-1)*2 + 2} = handles.data{k}.dyn;
end

function gain_slider1_Callback(hObject, eventdata, handles)
handles = SetGain(handles,hObject,1);
guidata(hObject, handles);

function gain_slider2_Callback(hObject, eventdata, handles)
handles = SetGain(handles,hObject,2);
guidata(hObject, handles);

function dyn_slider1_Callback(hObject, eventdata, handles)

handles = SetDynamicRange(handles,hObject,1);
guidata(hObject, handles);

function dyn_slider2_Callback(hObject, eventdata, handles)

handles = SetDynamicRange(handles,hObject,2);
guidata(hObject, handles);

function frame_slider1_Callback(hObject, eventdata, handles)

handles.frame_indx = get(hObject,'Value');
set(handles.frame_label(1),'String',sprintf('Frame # : %d',handles.frame_indx));

for k=1:handles.num_img
    UpdateImages(handles,k,handles.frame_indx);
end

guidata(hObject, handles);

function DoneButton_Callback(hObject, eventdata, handles)

uiresume(handles.figure1);
%guidata(hObject, handles);
%close(handles.figure1);

function handles = UpdateImages(handles,indx,frame_indx)

img = imagelog2(handles.data{indx}.env2(:,:,frame_indx),handles.data{indx}.gain,handles.data{indx}.dyn,256);    
axes(handles.axes_handles(indx))
image(handles.x_ax,handles.y_ax,img)
colormap(gray(256));

[histogram bins]= hist(img(:),(0:2:255) + 0.5);
%histogram = histogram./max(histogram(:));
axes(handles.axes_handles_hist(indx))
bar(bins(2:end),histogram(2:end));
axis tight


function handles = SetDynamicRange(handles,hObject,indx)

handles.data{indx}.dyn = get(hObject,'Value');
set(handles.dyn_label(indx),'String',sprintf('Dynamic Range : %2.1f dB',handles.data{indx}.dyn));
handles = UpdateImages(handles,indx,handles.frame_indx);

function handles = SetGain(handles,hObject,indx)

handles.data{indx}.gain = get(hObject,'Value');
set(handles.gain_label(indx),'String',sprintf('Gain : %2.1f dB',handles.data{indx}.gain));
handles = UpdateImages(handles,indx,handles.frame_indx);