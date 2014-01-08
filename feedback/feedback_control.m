function varargout = feedback_control(varargin)
% FEEDBACK_CONTROL M-file for feedback_control.fig
%      FEEDBACK_CONTROL, by itself, creates a new FEEDBACK_CONTROL or raises the existing
%      singleton*.
%
%      H = FEEDBACK_CONTROL returns the handle to a new FEEDBACK_CONTROL or the handle to
%      the existing singleton*.
%
%      FEEDBACK_CONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FEEDBACK_CONTROL.M with the given input arguments.
%
%      FEEDBACK_CONTROL('Property','Value',...) creates a new FEEDBACK_CONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before feedback_control_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to feedback_control_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help feedback_control

% Last Modified by GUIDE v2.5 26-Aug-2011 14:57:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @feedback_control_OpeningFcn, ...
                   'gui_OutputFcn',  @feedback_control_OutputFcn, ...
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


% --- Executes just before feedback_control is made visible.
function feedback_control_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to feedback_control (see VARARGIN)

% Choose default command line output for feedback_control
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes feedback_control wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = feedback_control_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global fbdata;
fbdata.pulseind = fbdata.buttonpls(1);

function pushbutton2_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.pulseind = fbdata.buttonpls(3);

function pushbutton3_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.pulseind = fbdata.buttonpls(5);
fbdata.tpihist={};
fbdata.tpi2hist={};

function pushbutton4_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.pulseind = fbdata.buttonpls(2);

function pushbutton5_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.pulseind = fbdata.buttonpls(4);

function pushbutton6_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.pulseind = [];

function slider1_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.ctrlval(end) = get(hObject,'Value');
set(handles.text1,'String',num2str(fbdata.ctrlval(end)));

%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global fbdata;
fbdata.ctrlval(end) = str2num(get(hObject,'String'));
set(handles.slider1,'Value',fbdata.ctrlval(end));
set(handles.text1,'String',num2str(fbdata.ctrlval(end)));

% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
pushbutton7_Callback(hObject,eventdata,handles);
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
pushbutton7_Callback(hObject,eventdata,handles);
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
pushbutton7_Callback(hObject,eventdata,handles);
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tauL.
function tauL_Callback(hObject, eventdata, handles)
% hObject    handle to tauL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tauL contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tauL
changeTau(handles);

% --- Executes during object creation, after setting all properties.
function tauL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tauL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global fbdata;
if isfield(fbdata,'tau')
  for i=1:length(fbdata.tau)
    str{i}=num2str(fbdata.tau(i));
  end
  set(hObject,'String',str);
else
  set(hObject,'String',{'1','2','3'});
end


% --- Executes on selection change in tauR.
function tauR_Callback(hObject, eventdata, handles)
% hObject    handle to tauR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns tauR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tauR
changeTau(handles);

% --- Executes during object creation, after setting all properties.
function tauR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tauR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global fbdata;
if isfield(fbdata,'tau')
  for i=1:length(fbdata.tau)
    str{i}=num2str(fbdata.tau(i));
  end
  set(hObject,'String',str);
else
  set(hObject,'String',{'1','2','3'});
end

function changeTau(handles)
  global fbdata;
  tl=get(handles.tauL,'value');
  tr=get(handles.tauR,'value');
  fbdata.taus=[fbdata.tau(tl) fbdata.tau(tr)];
  fbdata.buttonpls(5) = fbtau_plsind([fbdata.tau(tl),fbdata.tau(tr)]);
  fbdata.pulseind = fbdata.buttonpls(5);
  
  
