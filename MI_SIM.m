function varargout = MI_SIM(varargin)
% MI_SIM MATLAB code for MI_SIM.fig
%      MI_SIM, by itself, creates a new MI_SIM or raises the existing
%      singleton*.
%
%      H = MI_SIM returns the handle to a new MI_SIM or the handle to
%      the existing singleton*.
%
%      MI_SIM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MI_SIM.M with the given input arguments.
%
%      MI_SIM('Property','Value',...) creates a new MI_SIM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MI_SIM_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MI_SIM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MI_SIM

% Last Modified by GUIDE v2.5 27-Jan-2016 13:32:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MI_SIM_OpeningFcn, ...
    'gui_OutputFcn',  @MI_SIM_OutputFcn, ...
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


% --- Executes just before MI_SIM is made visible.
function MI_SIM_OpeningFcn(hObject, eventdata, handles, varargin)


%reset the plots
axes(handles.solutionplot);
cla reset
axes(handles.trajectoryplot);
cla reset

xlabel(handles.solutionplot,'Time (days)')
ylabel(handles.solutionplot,'Concentration (kgCOD m^{-3})')

%Default motif
set(handles.model_label,'String','4 ODE: Syntrophy');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','on'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
handles.motif='syn';
set(handles.twodimplot,'Checked','on'); set(handles.threedimplot,'Checked','off')
set(handles.timeplot,'Checked','on'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','off')
set(handles.sol_disp,'String','Display:')
handles.plotdim='two';
handles.plotsolution='time';
handles.growthmodel='Monod';
set(handles.overlay,'enable','off')
set(handles.s3_check,'value',0,'enable','off')

%reset all of the values
set(handles.growthmenu,'Value',1);
set(handles.km1_in,'String',13);
set(handles.y1_in,'String',0.04);
set(handles.kdec1_in,'String',0.02);
set(handles.km2_in,'String',35);
set(handles.y2_in,'String',0.06);
set(handles.kdec2_in,'String',0.02);
set(handles.s1in_in,'String',5);
set(handles.ki2_in,'String',0.0000035);
set(handles.ks1_in,'String',0.3);
set(handles.ks2_in,'String',0.000025);
set(handles.ks1a_in,'String',3);
set(handles.ks2a_in,'String',2);
set(handles.d_in,'String',0.1);
set(handles.time_in,'String',1000);
set(handles.s1_init,'String',0.1);
set(handles.x1_init,'String',0.1);
set(handles.s2_init,'String',0.1);
set(handles.x2_init,'String',0.1);
set(handles.fixed_points_s1,'String','-');
set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-');
set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_stability,'String','-');
set(handles.solver,'Value',3);
set(handles.abstol,'String','1e-8');
set(handles.reltol,'String','1e-8');
set(handles.error_val,'String','1e-4');
set(handles.pert_p,'String','0.0001');
set(handles.func_prog,'String','');
set(handles.lsanaly,'Value',1,'ForegroundColor',[0.078, 0.169, 0.549])
set(handles.routhcrit,'Value',0,'ForegroundColor','r')
set(handles.jacobian_but,'Value',0,'ForegroundColor','r')
set(handles.uipanel6,'Title','Plot of trajectory from initial conditions')
set(handles.spmatrix,'Visible','off'); set(handles.bifursp,'Visible','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off');
set(handles.normeig,'Visible','off','Value',0)
dir_fig=dir('temp_fig');
if length(dir_fig)<3
    set(handles.report_menu,'Enable','off')
    set(handles.reportpanel,'Visible','off'); set(handles.email_add,'Enable','off'); set(handles.server_add,'Enable','off')
    set(handles.email_txt,'Enable','off'); set(handles.serv_txt,'Enable','off');
else
    set(handles.report_menu,'Enable','on')
    set(handles.reportpanel,'Visible','on');set(handles.email_add,'Enable','on'); set(handles.server_add,'Enable','on')
    set(handles.email_txt,'Enable','on'); set(handles.serv_txt,'Enable','on');
end
handles.mpltopt=0;
handles.yout_phase=[];
handles.clcyc=1;
axes(handles.motif_image)
imshow('motifs_images/syn.png');

%Second organism parameter ON
set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

%Third organism parameter OFF
set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off'); set(handles.text81,'Enable','off')
set(handles.km3_in,'Enable','off','String',21);set(handles.y3_in,'Enable','off','String',0.04)
set(handles.kdec3_in,'Enable','off'); set(handles.ks3_in,'Enable','off','String',1e-3)
set(handles.ks32_in,'Enable','off','String',1e-6);

%Gamma values
set(handles.text42,'Enable','on'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off'); set(handles.gamma0,'Enable','on','String',0.43);
set(handles.gamma1,'Enable','off','String',0.1429); set(handles.gamma2,'Enable','off','String',0.0769);

%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1I_2$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1I_2 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma(1-Y_1)f_1X_1I_2 - f_2X_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';

f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';

handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2','KI2');
handles.var_names=strvcat('S1','X1','S2','X2');
handles.simtype='single_p';
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
fnctext = strvcat(f1tx,f2tx,I2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
% Choose default command line output for MI_SIM
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MI_SIM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MI_SIM_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


function x1_init_Callback(hObject, eventdata, handles)
X1_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function x1_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x2_init_Callback(hObject, eventdata, handles)
X2_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function x2_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s1_init_Callback(hObject, eventdata, handles)
S1_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s1_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function s2_init_Callback(hObject, eventdata, handles)
S2_init=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s2_init_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function s3_init_Callback(hObject, eventdata, handles)
S3_init=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s3_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function x3_init_Callback(hObject, eventdata, handles)
X3_init=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function x3_init_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in growthmenu.
function growthmenu_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
growth=contents{get(hObject,'Value')};
handles.growthmodel=growth;
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function growthmenu_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km1_in_Callback(hObject, eventdata, handles)

km1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function km1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function y1_in_Callback(hObject, eventdata, handles)
Y1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function y1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec1_in_Callback(hObject, eventdata, handles)
kdec1=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function kdec1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km2_in_Callback(hObject, eventdata, handles)
km2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function km2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y2_in_Callback(hObject, eventdata, handles)
Y2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function y2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec2_in_Callback(hObject, eventdata, handles)
kdec2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function kdec2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ki2_in_Callback(hObject, eventdata, handles)
KI2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function ki2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks1_in_Callback(hObject, eventdata, handles)
Ks1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks1_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks2_in_Callback(hObject, eventdata, handles)
Ks2=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function ks2_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks3_in_Callback(hObject, eventdata, handles)
Ks3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks1a_in_Callback(hObject, eventdata, handles)
Ks1a=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks1a_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks2a_in_Callback(hObject, eventdata, handles)
Ks2a=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks2a_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function km3_in_Callback(hObject, eventdata, handles)
km3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function km3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y3_in_Callback(hObject, eventdata, handles)
Y3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function y3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function kdec3_in_Callback(hObject, eventdata, handles)
kdec3=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function kdec3_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ks3a_in_Callback(hObject, eventdata, handles)
Ks3a=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks3a_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function d_in_Callback(hObject, eventdata, handles)
D=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function d_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s1in_in_Callback(hObject, eventdata, handles)
S1in=str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function s1in_in_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s2in_in_Callback(hObject, eventdata, handles)
S2in=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s2in_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function s3in_in_Callback(hObject, eventdata, handles)
S3in=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function s3in_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function time_in_Callback(hObject, eventdata, handles)
time1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function time_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gamma0_Callback(hObject, eventdata, handles)
gamma0=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma0_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ks32_in_Callback(hObject, eventdata, handles)
ks3c=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function ks32_in_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma1_Callback(hObject, eventdata, handles)
gamma1=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function gamma2_Callback(hObject, eventdata, handles)
gamma2=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function gamma2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in solver.
function solver_Callback(hObject, eventdata, handles)
solvstr=get(hObject,'String');
solvnum=get(hObject,'Value');
solver=solvstr(solvnum,:);
jacb=get(handles.jacobian_but,'Value');
switch strtrim(char(solver))
    case 'ode23s'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode15s'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode23t'
        set(handles.jacobian_but,'enable','on','Value',jacb)
    case 'ode23'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
    case 'ode45'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
    case 'ode113'
        set(handles.jacobian_but,'enable','off','Foregroundcolor','r','Value',0)
end

% --- Executes during object creation, after setting all properties.
function solver_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function abstol_Callback(hObject, eventdata, handles)
abstol=str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function abstol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function reltol_Callback(hObject, eventdata, handles)
reltol=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function reltol_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function error_val_Callback(hObject, eventdata, handles)
error=str2double(get(hObject,'string'));

% --- Executes during object creation, after setting all properties.
function error_val_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plot_button2.
function plot_button2_Callback(hObject, eventdata, handles)
which_plot='traj';
plot_results;


% --- Executes on button press in s1_check.
function s1_check_Callback(hObject, eventdata, handles)
s1_check=get(hObject,'Value');


% --- Executes on button press in x1_check.
function x1_check_Callback(hObject, eventdata, handles)
x1_check=get(hObject,'Value');


% --- Executes on button press in s2_check.
function s2_check_Callback(hObject, eventdata, handles)
s2_check=get(hObject,'Value');


% --- Executes on button press in x2_check.
function x2_check_Callback(hObject, eventdata, handles)
x2_check=get(hObject,'Value');

% --- Executes on button press in s3_check.
function s3_check_Callback(hObject, eventdata, handles)
s3_check=get(hObject,'Value');

% --- Executes on button press in x3_check.
function x3_check_Callback(hObject, eventdata, handles)
x3_check=get(hObject,'Value');

% --- Executes on button press in overlay.
function overlay_Callback(hObject, eventdata, handles)
multi=get(hObject,'Value');

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)

% Get user input
growth=get(handles.growthmenu,'Value');
motif=handles.motif;

km1 = str2double(get(handles.km1_in,'String'));
if km1<0 || isnumeric(km1)==0
    set(handles.km1_in,'String',13);
    km1 = str2double(get(handles.km1_in,'String'));
    msgbox('The value entered for km1 was not valid so the default value was reset','Error','error')
else
end
Y1 = str2double(get(handles.y1_in,'String'));
if Y1<0 || isnumeric(Y1)==0
    set(handles.y1_in,'String',0.04);
    Y1 = str2double(get(handles.y1_in,'String'));
    msgbox('The value entered for Y1 was not valid so the default value was reset','Error','error')
else
end
kdec1 = str2double(get(handles.kdec1_in,'String'));
if kdec1<0 || kdec1>1 || isnumeric(kdec1)==0
    set(handles.kdec1_in,'String',0.02);
    kdec1 = str2double(get(handles.kdec1_in,'String'));
    msgbox('The value entered for (0<)kdec1(<1) was not valid so the default value was reset','Error','error')
else
end
km2 = str2double(get(handles.km2_in,'String'));
if km2<0 || isnumeric(km2)==0
    set(handles.km2_in,'String',35);
    km2 = str2double(get(handles.km2_in,'String'));
    msgbox('The value entered for km2 was not valid so the default value was reset','Error','error')
else
end
Y2 = str2double(get(handles.y2_in,'String'));
if Y2<0 || isnumeric(Y2)==0
    set(handles.y2_in,'String',0.06);
    Y2 = str2double(get(handles.y2_in,'String'));
    msgbox('The value entered for Y2 was not valid so the default value was reset','Error','error')
else
end
kdec2 = str2double(get(handles.kdec2_in,'String'));
if kdec2<0 || kdec2>1 || isnumeric(kdec2)==0
    set(handles.kdec2_in,'String',0.02);
    kdec2 = str2double(get(handles.kdec2_in,'String'));
    msgbox('The value entered for (0<)kdec2(<1) was not valid so the default value was reset','Error','error')
else
end

km3 = str2double(get(handles.km3_in,'String'));
if km3<0 || isnumeric(km3)==0
    set(handles.km3_in,'String',21);
    km3 = str2double(get(handles.km3_in,'String'));
    msgbox('The value entered for km3 was not valid so the default value was reset','Error','error')
else
end
Y3 = str2double(get(handles.y3_in,'String'));
if Y3<0 || isnumeric(Y3)==0
    set(handles.y3_in,'String',0.04);
    Y3 = str2double(get(handles.y3_in,'String'));
    msgbox('The value entered for Y3 was not valid so the default value was reset','Error','error')
else
end
kdec3 = str2double(get(handles.kdec3_in,'String'));
if kdec3<0 || kdec3>1 || isnumeric(kdec3)==0
    set(handles.kdec3_in,'String',0.02);
    kdec3 = str2double(get(handles.kdec3_in,'String'));
    msgbox('The value entered for (0<)kdec3(<1) was not valid so the default value was reset','Error','error')
else
end

S1in = str2double(get(handles.s1in_in,'String'));
if S1in<0 || isnumeric(S1in)==0
    set(handles.s1in_in,'String',5);
    S1in = str2double(get(handles.s1in_in,'String'));
    msgbox('The value entered for S1in was not valid so the default value was reset','Error','error')
else
end
S2in = str2double(get(handles.s2in_in,'String'));
if S2in<0 || isnumeric(S2in)==0
    set(handles.s2in_in,'String',5);
    S2in = str2double(get(handles.s2in_in,'String'));
    msgbox('The value entered for S2in was not valid so the default value was reset','Error','error')
else
end

S3in = str2double(get(handles.s3in_in,'String'));
if S3in<0 || isnumeric(S3in)==0
    set(handles.s3in_in,'String',1);
    S3in = str2double(get(handles.s3in_in,'String'));
    msgbox('The value entered for S3in was not valid so the default value was reset','Error','error')
else
end
KI2 = str2double(get(handles.ki2_in,'String'));
if KI2<0 || isnumeric(KI2)==0
    set(handles.ki2_in,'String',0.0000035);
    KI2 = str2double(get(handles.ki2_in,'String'));
    msgbox('The value entered for KI was not valid so the default value was reset','Error','error')
else
end
Ks1 = str2double(get(handles.ks1_in,'String'));
if Ks1<0 || isnumeric(Ks1)==0
    set(handles.ks1_in,'String',0.3);
    Ks1 = str2double(get(handles.ks1_in,'String'));
    msgbox('The value entered for KS1) was not valid so the default value was reset','Error','error')
else
end
Ks2 = str2double(get(handles.ks2_in,'String'));
if Ks2<0 || isnumeric(Ks2)==0
    set(handles.ks2_in,'String',0.000025);
    Ks2 = str2double(get(handles.ks2_in,'String'));
    msgbox('The value entered for KS2 was not valid so the default value was reset','Error','error')
else
end

Ks3 = str2double(get(handles.ks3_in,'String'));
if Ks3<0 || isnumeric(Ks3)==0
    set(handles.ks3_in,'String',0.001);
    Ks3 = str2double(get(handles.ks3_in,'String'));
    msgbox('The value entered for KS3 was not valid so the default value was reset','Error','error')
else
end

Ks3c = str2double(get(handles.ks32_in,'String'));
if Ks3c<0 || isnumeric(Ks3c)==0
    set(handles.ks32_in,'String',1e-6);
    Ks3c = str2double(get(handles.ks32_in,'String'));
    msgbox('The value entered for KS3c was not valid so the default value was reset','Error','error')
else
end

Ks1a = str2double(get(handles.ks1a_in,'String'));
if Ks1a<0 || isnumeric(Ks1a)==0
    set(handles.ks1a_in,'String',3);
    Ks1a = str2double(get(handles.ks1a_in,'String'));
    msgbox('The value entered for KS1c was not valid so the default value was reset','Error','error')
else
end
Ks2a = str2double(get(handles.ks2a_in,'String'));
if Ks2a<0 || isnumeric(Ks2a)==0
    set(handles.ks2a_in,'String',2);
    Ks2a = str2double(get(handles.ks2a_in,'String'));
    msgbox('The value entered for KS2c was not valid so the default value was reset','Error','error')
else
end
Ks3a = str2double(get(handles.ks3a_in,'String'));
if Ks3a<0 || isnumeric(Ks3a)==0
    set(handles.ks3a_in,'String',2);
    Ks3a = str2double(get(handles.ks3a_in,'String'));
    msgbox('The value entered for KS3c was not valid so the default value was reset','Error','error')
else
end
gamma0 = str2double(get(handles.gamma0,'String'));
if gamma0<0 || isnumeric(gamma0)==0
    set(handles.gamma0,'String',0.43);
    gamma0 = str2double(get(handles.gamma0,'String'));
    msgbox('The value entered for \gamma_0 was not valid so the default value was reset','Error','error')
else
end
gamma1 = str2double(get(handles.gamma1,'String'));
if gamma1<0 || isnumeric(gamma1)==0
    set(handles.gamma1,'String',0.1429);
    gamma1 = str2double(get(handles.gamma1,'String'));
    msgbox('The value entered for \gamma_1 was not valid so the default value was reset','Error','error')
else
end
gamma2 = str2double(get(handles.gamma2,'String'));
if gamma2<0 || isnumeric(gamma2)==0
    set(handles.gamma2,'String',0.0769);
    gamma2 = str2double(get(handles.gamma2,'String'));
    msgbox('The value entered for \gamma_2 was not valid so the default value was reset','Error','error')
else
end
D = str2double(get(handles.d_in,'String'));
if D<0 || isnumeric(D)==0
    set(handles.d_in,'String',0.1);
    D = str2double(get(handles.d_in,'String'));
    msgbox('The value entered for dilution(>0) was not valid so the default value was reset','Error','error')
else
end
time1 = str2double(get(handles.time_in,'String'));
if time1<0 || isnumeric(time1)==0
    set(handles.time_in,'String',1000);
    time1 = str2double(get(handles.time_in,'String'));
    msgbox('The value entered for time(>0) was not valid so the default value was reset','Error','error')
else
end

S1_init=str2double(get(handles.s1_init,'String'));
if S1_init<0 || isnumeric(S1_init)==0
    set(handles.s1_init,'String',0.1);
    S1_init=str2double(get(handles.s1_init,'String'));
    msgbox('The value entered for S1_init was not valid so the default value was reset','Error','error')
else
end
X1_init=str2double(get(handles.x1_init,'String'));
if X1_init<0 || isnumeric(X1_init)==0
    set(handles.x1_init,'String',0.1);
    X1_init=str2double(get(handles.x1_init,'String'));
    msgbox('The value entered for X1_init was not valid so the default value was reset','Error','error')
else
end
S2_init=str2double(get(handles.s2_init,'String'));
if S2_init<0 || isnumeric(S2_init)==0
    set(handles.s2_init,'String',0.1);
    S2_init=str2double(get(handles.s2_init,'String'));
    msgbox('The value entered for S2_init was not valid so the default value was reset','Error','error')
else
end
X2_init=str2double(get(handles.x2_init,'String'));
if X2_init<0 || isnumeric(X2_init)==0
    set(handles.x2_init,'String',0.1);
    X2_init=str2double(get(handles.x2_init,'String'));
    msgbox('The value entered for X2_init was not valid so the default value was reset','Error','error')
else
end
S3_init=str2double(get(handles.s3_init,'String'));
if S3_init<0 || isnumeric(S3_init)==0
    set(handles.s3_init,'String',0.1);
    S3_init=str2double(get(handles.s3_init,'String'));
    msgbox('The value entered for S3_init was not valid so the default value was reset','Error','error')
else
end
X3_init=str2double(get(handles.x3_init,'String'));
if X3_init<0 || isnumeric(X3_init)==0
    set(handles.x3_init,'String',0.1);
    X3_init=str2double(get(handles.x3_init,'String'));
    msgbox('The value entered for X3_init was not valid so the default value was reset','Error','error')
else
end

abstol=str2double(get(handles.abstol,'String'));
if abstol<0 || isnumeric(abstol)==0
    set(handles.abstol,'String',1e-8);
    abstol=str2double(get(handles.abstol,'String'));
    msgbox('The value entered for Abs. tolerance was not valid so the default value was reset','Error','error')
else
end

reltol=str2double(get(handles.reltol,'String'));
if reltol<0 || isnumeric(reltol)==0
    set(handles.reltol,'String',1e-8);
    reltol=str2double(get(handles.reltol,'String'));
    msgbox('The value entered for Rel. tolerance was not valid so the default value was reset','Error','error')
else
end

solvstr=get(handles.solver,'String');
solvnum=get(handles.solver,'Value');
if solvnum==1
    solvnum=3;
end
solver=solvstr(solvnum,:);

s1_check=get(handles.s1_check,'Value');
x1_check=get(handles.x1_check,'Value');
s2_check=get(handles.s2_check,'Value');
x2_check=get(handles.x2_check,'Value');
s3_check=get(handles.s3_check,'Value');
x3_check=get(handles.x3_check,'Value');

run_script_gui;

function fixed_points_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function fixed_points_s1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function models_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function commensalism_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','4 ODE: Food Chain');
set(handles.commensalism,'Checked','on'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='fc';
axes(handles.motif_image)
imshow('motifs_images/fc.eps');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','off'); set(handles.s3_init_text,'enable','off')
set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','on')
set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-')
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-')
set(handles.fixed_points_s3,'String','-');

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
set(handles.km3_in,'Enable','off');set(handles.y3_in,'Enable','off')
set(handles.kdec3_in,'Enable','off');set(handles.ks3_in,'Enable','off')
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

handles.plotsolution='time'; handles.plotdim='two';

%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma(1-Y_1)f_1X_1 - f_2X_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';

%Functions depend on growth model
switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2}X_2 + S_2}$';
end

handles.params_sim=strvcat('S1in','D','km1','KS1','Y1','km2','KS2','Y2','kdec1','kdec2');
handles.var_names=strvcat('S1','X1','S2','X2');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
fnctext = strvcat(f1tx,f2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
guidata(hObject,handles);

% --------------------------------------------------------------------
function competition_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','3 ODE: Substrate Competition');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','on');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='sc';
axes(handles.motif_image)
imshow('motifs_images/sc.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','off'); set(handles.s3_init_text,'enable','off')
set(handles.s2_init,'enable','off');set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','off')
set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_s3,'String','-');
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
set(handles.km3_in,'Enable','off');set(handles.y3_in,'Enable','off')
set(handles.kdec3_in,'Enable','off');set(handles.ks3_in,'Enable','off')
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1 - f_2X_2$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
eqtx3='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';

switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2}X_2 + S_2}$';
end

handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2');
handles.var_names=strvcat('S1','X1','X2');
eqtext = strvcat(eqtx1,eqtx2,eqtx3);
fnctext = strvcat(f1tx,f2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;

guidata(hObject,handles);

% --------------------------------------------------------------------
function predation_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','5 ODE: Food Chain with Waste Product Inhibition');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','on'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='fci';
axes(handles.motif_image)
imshow('motifs_images/fcpi.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','on'); set(handles.s3_init_text,'enable','on')
set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','on')
set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_s3,'String','-');
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
set(handles.km3_in,'Enable','off');set(handles.y3_in,'Enable','off')
set(handles.kdec3_in,'Enable','off');set(handles.ks3_in,'Enable','off')
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1I_3$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1I_3 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma(1-Y_1)f_1X_1I_3 - f_2X_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
eqtx5='$\frac{dS_3}{dt} = -DS_3 + f_2X_2$';
I3tx = '$I_3 = \frac{1}{1+\frac{S_3}{K_{i,3}}}$';
%Functions depend on growth model
switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2}X_2 + S_2}$';
end

handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2','KI3');
handles.var_names=strvcat('S1','X1','S2','X2','S3');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5);
fnctext = strvcat(f1tx,f2tx,I3tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;

guidata(hObject,handles)

% --------------------------------------------------------------------
function no_interaction_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','4 ODE: No Common Metabolites');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','on');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='ni';
axes(handles.motif_image)
imshow('motifs_images/ni.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','off'); set(handles.s3_init_text,'enable','off')
set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','on')
set(handles.s2in_in,'enable','on'); set(handles.text55,'enable','on')
set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_s3,'String','-');
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
set(handles.km3_in,'Enable','off');set(handles.y3_in,'Enable','off')
set(handles.kdec3_in,'Enable','off');set(handles.ks3_in,'Enable','off')
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = D(S_{2,in} - S_2) - f_2X_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';

switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2}X_2 + S_2}$';
end

handles.params_sim=strvcat('S1in','D','S2in','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2');
handles.var_names=strvcat('S1','X1','S2','X2');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
fnctext = strvcat(f1tx,f2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
guidata(hObject,handles)

% --------------------------------------------------------------------
function cooperation_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','4 ODE: Syntrophy');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','on'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='syn';
axes(handles.motif_image)
imshow('motifs_images/syn.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',0,'enable','off'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','off'); set(handles.s3_init_text,'enable','off')
set(handles.s2_init,'enable','on');set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','on')
set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
set(handles.s3in_in,'enable','off'); set(handles.text48,'enable','off')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_s3,'String','-');
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on')
set(handles.km2_in,'Enable','on','String',35);set(handles.y2_in,'Enable','on','String',0.06)
set(handles.kdec2_in,'Enable','on');set(handles.ks2_in,'Enable','on','String',2.5e-5)

set(handles.text43,'Enable','off');set(handles.text44,'Enable','off')
set(handles.text45,'Enable','off');set(handles.text47,'Enable','off')
set(handles.km3_in,'Enable','off');set(handles.y3_in,'Enable','off')
set(handles.kdec3_in,'Enable','off');set(handles.ks3_in,'Enable','off')
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1I_2$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1I_2 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma(1-Y_1)f_1X_1I_2 - f_2X_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2 - k_{dec,2}X_2$';
I2tx = '$I_2 = \frac{1}{1+\frac{S_2}{K_{i,2}}}$';

switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2}X_2 + S_2}$';
end

handles.params_sim=strvcat('S1in','D','km1','Ks1','Y1','km2','Ks2','Y2','kdec1','kdec2','KI2');
handles.var_names=strvcat('S1','X1','S2','X2');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4);
fnctext = strvcat(f1tx,f2tx,I2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
guidata(hObject,handles)

% --------------------------------------------------------------------
function amensalism_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','5 ODE: Waste Product Inhibition');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','on'); set(handles.threesp,'Checked','off');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='pi';
axes(handles.motif_image)
imshow('motifs_images/wpi.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',1,'enable','on'); set(handles.x3_check,'value',0,'enable','off')
set(handles.s3_init,'enable','on'); set(handles.s3_init_text,'enable','on')
set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','off') 
set(handles.x3_init_text,'enable','off'); set(handles.S2_0,'enable','on')
set(handles.s2in_in,'enable','off'); set(handles.text55,'enable','off')
set(handles.s3in_in,'enable','on'); set(handles.text48,'enable','on')
set(handles.ks32_in,'enable','off','String',1e-6); set(handles.text81,'enable','off')
set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off');
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-');
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-');
set(handles.fixed_points_s3,'String','-');
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','off');set(handles.text15,'Enable','off')
set(handles.text16,'Enable','off');set(handles.text19,'Enable','off')
set(handles.km2_in,'Enable','off');set(handles.y2_in,'Enable','off')
set(handles.kdec2_in,'Enable','off');set(handles.ks2_in,'Enable','off')

set(handles.text43,'Enable','on');set(handles.text44,'Enable','on')
set(handles.text45,'Enable','on');set(handles.text47,'Enable','on')
set(handles.km3_in,'Enable','on','Value',21);set(handles.y3_in,'Enable','on','String',0.04)
set(handles.kdec3_in,'Enable','on');set(handles.ks3_in,'Enable','on','String',1e-3)
set(handles.gamma0,'Enable','on','String',0.43); set(handles.gamma1,'Enable','off')
set(handles.gamma2,'Enable','off'); set(handles.text82,'Enable','off')
set(handles.text83,'Enable','off')

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = -DS_2 + \gamma(1-Y_1)f_1X_1$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_3f_3X_2I_4 - k_{dec,3}X_2$';
eqtx5='$\frac{dS_3}{dt} = D(S_{3,in}-S_3) - f_3X_2I_4$';
I3tx = '$I_4 = \frac{1}{1+\frac{K_{i,2}}{S_2}}$';

switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_1}{K_{S,1}X_1 + S_1}$';
        f2tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3}X_3 + S_3}$';
end


handles.params_sim=strvcat('S1in','D','S3in','km1','Ks1','Y1','km3','Ks3','Y3','kdec1','kdec2','KI2');
handles.var_names=strvcat('S1','X1','S2','X2','S3');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5);
fnctext = strvcat(f1tx,f2tx,I3tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.5,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
text(0.7,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',14)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
guidata(hObject,handles)

% --------------------------------------------------------------------
function threesp_Callback(hObject, eventdata, handles)
set(handles.model_label,'String','6 ODE: 3 species food-web');
set(handles.commensalism,'Checked','off'); set(handles.competition,'Checked','off');
set(handles.predation,'Checked','off'); set(handles.no_interaction,'Checked','off');
set(handles.cooperation,'Checked','off'); set(handles.amensalism,'Checked','off'); set(handles.threesp,'Checked','on');
set(handles.spmatrix,'Visible','off'); axes(handles.spmatrix); cla
handles.motif='ths';
axes(handles.motif_image)
imshow('motifs_images/ths.png');
set(handles.plot_button,'enable','off')
set(handles.plot_button2,'enable','off')
set(handles.s3_check,'value',1,'enable','on'); set(handles.x3_check,'value',1,'enable','on')
set(handles.s3_init,'enable','on'); set(handles.s3_init_text,'enable','on')
set(handles.s2_init,'enable','on'); set(handles.x3_init,'enable','on') 
set(handles.x3_init_text,'enable','on'); set(handles.S2_0,'enable','on')
set(handles.ks2_in,'enable','on','String',0.302); set(handles.ks1_in,'String',5.2e-5)
set(handles.s2in_in,'enable','on','String',0); set(handles.text55,'enable','on')
set(handles.s3in_in,'enable','on','String',0); set(handles.text48,'enable','on')
set(handles.km3_in,'enable','on','String',35); set(handles.text43,'enable','on')
set(handles.y3_in,'enable','on','String',0.06); set(handles.text44,'enable','on')
set(handles.kdec3_in,'enable','on'); set(handles.text45,'enable','on')
set(handles.ks3_in,'enable','on','String',2.5e-5); set(handles.text47,'enable','on')
set(handles.ks32_in,'enable','on','String',1e-6); set(handles.text81,'enable','on')
set(handles.gamma0,'String',1.0769); set(handles.gamma1,'enable','on','String',0.1429) 
set(handles.gamma2,'enable','on','String',0.0769); set(handles.text82,'enable','on')
set(handles.text83,'enable','on');

set(handles.timeplot,'checked','on'); set(handles.phaseplot,'checked','off')
set(handles.fixed_points_s1,'String','-'); set(handles.fixed_points_x1,'String','-')
set(handles.fixed_points_s2,'String','-'); set(handles.fixed_points_x2,'String','-')
set(handles.fixed_points_s3,'String','-'); set(handles.fixed_points_x3,'String','-')
handles.plotsolution='time'; handles.plotdim='two';

set(handles.text14,'Enable','on');set(handles.text15,'Enable','on')
set(handles.text16,'Enable','on');set(handles.text19,'Enable','on'); set(handles.text55,'Enable','on')
set(handles.km2_in,'Enable','on','Value',26);set(handles.y2_in,'Enable','on','String',0.04)
set(handles.kdec2_in,'Enable','on');

axes(handles.trajectoryplot)
colorbar off
%Equations
eqtx1='$\frac{dS_1}{dt} = D(S_{1,in} - S_1) - f_1X_1$';
eqtx2='$\frac{dX_1}{dt} = -DX_1 + Y_1f_1X_1 - k_{dec,1}X_1$';
eqtx3='$\frac{dS_2}{dt} = D(S_{2,in} - S_2) + \gamma_0(1-Y_1)f_1X_1 - f_2X_2I_2$';
eqtx4='$\frac{dX_2}{dt} = -DX_2 + Y_2f_2X_2I_2 - k_{dec,2}X_2$';
eqtx5='$\frac{dS_3}{dt} = D(S_{3,in}-S_3) + \gamma_1(1-Y_2)f_2X_2I_2 - f_3X_3 - \gamma_2f_1X_1$';
eqtx6='$\frac{dX_3}{dt} = -DX_3 + Y_3f_3X_3 - k_{dec,3}X_3$';
I2tx = '$I_2 = \frac{1}{1+\frac{S_3}{K_{i,2}}}$';

switch handles.growthmodel
    case 'Monod'
        f1tx = '$f_1 = \frac{k_{m,1}S_3}{K_{S,3c} + S_3}\frac{S_1}{K_{S,1} + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2} + S_2}$';
        f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3} + S_3}$';
    case 'Contois'
        f1tx = '$f_1 = \frac{k_{m,1}S_3}{K_{S,3c}X_3 + S_3}\frac{S_1}{K_{S,1co}X_1 + S_1}$';
        f2tx = '$f_2 = \frac{k_{m,2}S_2}{K_{S,2co}X_2 + S_2}$';
        f3tx = '$f_3 = \frac{k_{m,3}S_3}{K_{S,3co}X_3 + S_3}$';
end

handles.params_sim=strvcat('S1in','D','S3in','S2in','km1','Ks1','Y1','km2','Ks2','Y2','km3','Ks3','Y3','kdec1','kdec2','kdec3','Ks3c','KI2');
handles.var_names=strvcat('S1','X1','S2','X2','S3','X3');
eqtext = strvcat(eqtx1,eqtx2,eqtx3,eqtx4,eqtx5,eqtx6);
fnctext = strvcat(f1tx,f2tx,f3tx,I2tx);
cla(handles.txax)
axes(handles.txax)

text(0,0.45,eqtext,'interpreter','latex','horiz','left','vert','middle','fontsize',12)
text(0.81,0.5,fnctext,'interpreter','latex','horiz','left','vert','middle','fontsize',12)
handles.eqtx=eqtext;
handles.fnctx=fnctext;
guidata(hObject,handles)


% --------------------------------------------------------------------
function plotoptions_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function solutionsplot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function trajplot_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function twodimplot_Callback(hObject, eventdata, handles)
set(handles.twodimplot,'Checked','on'); set(handles.threedimplot,'Checked','off')
handles.plotdim='two';
set(handles.s1_check,'value',1)
set(handles.x1_check,'value',1)
set(handles.s2_check,'value',0)
set(handles.x2_check,'value',0)
set(handles.s3_check,'value',0)
set(handles.x3_check,'value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function threedimplot_Callback(hObject, eventdata, handles)
set(handles.twodimplot,'Checked','off'); set(handles.threedimplot,'Checked','on')
handles.plotdim='three';
set(handles.s1_check,'value',1)
set(handles.x1_check,'value',1)
set(handles.s2_check,'value',1)
set(handles.x2_check,'value',0)
set(handles.s3_check,'value',0)
set(handles.x3_check,'value',0)

guidata(hObject,handles)

% --------------------------------------------------------------------
function timeplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','on'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','off')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='time';
set(handles.sol_disp,'String','Display:')
set(handles.s1_check,'Value',1,'Enable','on'); set(handles.x1_check,'Value',1,'Enable','on')
set(handles.s2_check,'Value',1,'Enable','on'); set(handles.x2_check,'Value',1,'Enable','on')
set(handles.normeig,'Visible','off','Value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function phaseplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','on'); set(handles.eigplot,'Checked','off')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='phase';
set(handles.sol_disp,'String','Select 2-3 variables')
set(handles.s1_check,'Value',0,'Enable','on'); set(handles.x1_check,'Value',0,'Enable','on')
set(handles.s2_check,'Value',0,'Enable','on'); set(handles.x2_check,'Value',0,'Enable','on')
set(handles.normeig,'Visible','off','Value',0)
guidata(hObject,handles)

% --------------------------------------------------------------------
function subplot_sol_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','off'); set(handles.subplot_sol,'Checked','on')
set(handles.eigplot,'Checked','off')
handles.plotsolution='subp';
set(handles.sol_disp,'String','Display')
set(handles.s1_check,'Value',1,'Enable','on'); set(handles.x1_check,'Value',1,'Enable','on')
set(handles.s2_check,'Value',1,'Enable','on'); set(handles.x2_check,'Value',1,'Enable','on')
set(handles.normeig,'Visible','on','Value',0,'String','Overlay plots')
guidata(hObject,handles)

% --------------------------------------------------------------------
function eigplot_Callback(hObject, eventdata, handles)
set(handles.timeplot,'Checked','off'); set(handles.phaseplot,'Checked','off'); set(handles.eigplot,'Checked','on')
set(handles.subplot_sol,'Checked','off')
handles.plotsolution='eig';
set(handles.sol_disp,'String','Eigenvalues')
set(handles.s1_check,'Value',0,'Enable','off'); set(handles.x1_check,'Value',0,'Enable','off')
set(handles.s2_check,'Value',0,'Enable','off'); set(handles.x2_check,'Value',0,'Enable','off')
set(handles.s3_check,'Value',0,'Enable','off'); set(handles.x3_check,'Value',0,'Enable','off')
set(handles.normeig,'Visible','on','Value',0,'String','Normalise')
guidata(hObject,handles)

% --- Executes on button press in plot_button.
function plot_button_Callback(hObject, eventdata, handles)
which_plot='sol';
plot_results;

% --- Executes on button press in lsanaly.
function lsanaly_Callback(hObject, eventdata, handles)
gvlsa=get(handles.lsanaly,'Value');
if gvlsa==0
    set(handles.lsanaly,'ForegroundColor','r')
elseif gvlsa==1
    set(handles.lsanaly,'ForegroundColor',[0.078, 0.169, 0.549])
end


% --- Executes on button press in routhcrit.
function routhcrit_Callback(hObject, eventdata, handles)
gvrhc=get(handles.routhcrit,'Value');
if gvrhc==0
    set(handles.routhcrit,'ForegroundColor','r')
elseif gvrhc==1
    set(handles.routhcrit,'ForegroundColor',[0.078, 0.169, 0.549])
end
   

% --- Executes on button press in jacobian_but.
function jacobian_but_Callback(hObject, eventdata, handles)
gvjac=get(handles.jacobian_but,'Value');
if gvjac==0
    set(handles.jacobian_but,'ForegroundColor','r')
elseif gvjac==1
    set(handles.jacobian_but,'ForegroundColor',[0.078, 0.169, 0.549])
end

function pert_p_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pert_p_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function sim_opts_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function sim_sp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','on'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor',[0.502 0.502 0.502])
set(handles.param1txt,'Enable','off'); set(handles.param2txt,'Enable','off')
set(handles.simparam1,'Enable','off','String','.'); set(handles.simparam2,'Enable','off','String','.'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.min_p1,'Enable','off','String','','Visible','on'); set(handles.max_p1,'Enable','off','String','','Visible','on'); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.min_p2,'Enable','off','String','','Visible','on'); set(handles.max_p2,'Enable','off','String','','Visible','on');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Enable','off','String','','Visible','on'); set(handles.step_p2,'Enable','off','String','','Visible','on','Style','edit');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Plot of trajector from initial conditions'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off')
handles.simtype='single_p';

guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_mp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','on')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor','r'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Parameter 1'); set(handles.param2txt,'Enable','on','String','Parameter 2'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.simparam1,'Enable','on','String',handles.params_sim); set(handles.simparam2,'Enable','on','String',handles.params_sim); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.min_p1,'Enable','on','String','','Visible','on'); set(handles.max_p1,'Enable','on','String','','Visible','on');
set(handles.min_p2,'Enable','on','String','','Visible','on'); set(handles.max_p2,'Enable','on','String','','Visible','on');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Enable','on','String','','Visible','on'); set(handles.step_p2,'Enable','on','String','','Visible','on','Style','edit');
set(handles.min_txt,'Enable','on','Visible','on','String','Minimum'); set(handles.max_txt,'Enable','on','Visible','on','String','Maximum');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Bifurcation phase plot'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
handles.simtype='multiple_p';
guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_boa_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','on'); set(handles.sim_pp,'Checked','off')
set(handles.simpanel,'ForegroundColor','b'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Variable 1'); set(handles.param2txt,'Enable','on','String','Variable 2'); set(handles.var3_txt,'Visible','off','Enable','off')
set(handles.simparam1,'Enable','on','String',handles.var_names); set(handles.simparam2,'Enable','on','String',handles.var_names); set(handles.simparam3,'Visible','off','Enable','off','String','');
set(handles.simparam3,'Visible','off','Enable','off','String',handles.var_names);
set(handles.min_p1,'Visible','on','Enable','on','String','0'); set(handles.max_p1,'Visible','on','Enable','on','String','1');
set(handles.min_p2,'Visible','off','String',''); set(handles.max_p2,'Visible','off','String','');
set(handles.min_p3,'Visible','off','String',''); set(handles.max_p3,'Visible','off','String','');
set(handles.step_p1,'Visible','on','Enable','on','String','50'); set(handles.step_p2,'Visible','off','String','','Style','edit');
set(handles.min_txt,'Visible','on','Enable','on','String','Lower Limit'); set(handles.max_txt,'Visible','on','Enable','on','String','Upper Limit');
set(handles.use_3v,'Visible','off','Value',0);
set(handles.uipanel6,'Title','Basin of Attraction'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
handles.simtype='boa';

guidata(hObject,handles)

% --------------------------------------------------------------------
function sim_pp_Callback(hObject, eventdata, handles)
set(handles.sim_sp,'Checked','off'); set(handles.sim_mp,'Checked','off')
set(handles.sim_boa,'Checked','off'); set(handles.sim_pp,'Checked','on')
set(handles.simpanel,'ForegroundColor','b'); set(handles.step_txt,'Enable','on');
set(handles.param1txt,'Enable','on','String','Variable 1'); set(handles.param2txt,'Enable','on','String','Variable 2'); set(handles.var3_txt,'Visible','off','Enable','on')
set(handles.simparam1,'Enable','on','String',handles.var_names); set(handles.simparam2,'Enable','on','String',handles.var_names); set(handles.simparam3,'Visible','off','Enable','on','String','');
set(handles.min_p1,'Visible','on','Enable','on','String','0'); set(handles.max_p1,'Visible','on','Enable','on','String','1');
set(handles.min_p2,'Visible','on','Enable','on','String','0'); set(handles.max_p2,'Visible','on','Enable','on','String','1');
set(handles.min_p3,'Visible','off','Enable','on','String',''); set(handles.max_p3,'Visible','off','Enable','on','String','');
set(handles.step_p1,'Visible','on','Enable','on','String','50'); set(handles.step_p2,'Visible','on','Enable','on','String',strvcat('Fixed','Random'),'Style','popupmenu','Value',1);
set(handles.min_txt,'Visible','on','Enable','on','String','Lower IC'); set(handles.max_txt,'Visible','on','Enable','on','String','Upper IC');
set(handles.use_3v,'Visible','on','Value',0);
set(handles.uipanel6,'Title','Phase Portrait'); set(handles.overlay,'enable','off')
set(handles.twodphase,'Enable','off'); set(handles.threedphase,'Enable','off')
set(handles.colormp,'Visible','off')
handles.simtype='pport';

guidata(hObject,handles)


% --------------------------------------------------------------------
function gui_opts_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function gui_reset_Callback(hObject, eventdata, handles)
%reset all of the values
set(handles.km1_in,'String',13);
set(handles.y1_in,'String',0.04);
set(handles.kdec1_in,'String',0.02);
set(handles.km2_in,'String',35);
set(handles.y2_in,'String',0.06);
set(handles.kdec2_in,'String',0.02);
set(handles.s1in_in,'String',5);
set(handles.ki2_in,'String',0.0000035);
set(handles.ks1_in,'String',0.3);
set(handles.ks2_in,'String',0.000025);
set(handles.ks1a_in,'String',3);
set(handles.ks2a_in,'String',2);
set(handles.d_in,'String',0.1);
set(handles.time_in,'String',1000);
set(handles.s1_init,'String',0.1);
set(handles.x1_init,'String',0.1);
set(handles.s2_init,'String',0.1);
set(handles.x2_init,'String',0.1);
set(handles.km3_in,'String',21);
set(handles.y3_in,'String',0.04);
set(handles.kdec3_in,'String',0.02);
set(handles.ks3_in,'String',0.001);
set(handles.ks3a_in,'String',2);
set(handles.gamma0,'String',0.43);
set(handles.gamma1,'String',0.1429);
set(handles.gamma2,'String',0.0769);
set(handles.ks32_in,'String',1e-6);
set(handles.abstol,'String',1e-8);
set(handles.reltol,'String',1e-8);
set(handles.error_val,'String','1e-4');
set(handles.solver,'value',3);
set(handles.pert_p,'value',0.00001);
set(handles.normeig,'Visible','off','Value',0)


% --------------------------------------------------------------------
function gui_restart_Callback(hObject, eventdata, handles)
close(gcbf)
MI_SIM

% --------------------------------------------------------------------
function gui_close_Callback(hObject, eventdata, handles)
close(gcbf)


% --- Executes on selection change in simparam1.
function simparam1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function simparam1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in simparam2.
function simparam2_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function simparam2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p1_Callback(hObject, eventdata, handles)
minp1a = get(handles.min_p1,'String');
minp1 = str2double(minp1a);

% --- Executes during object creation, after setting all properties.
function min_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_p1_Callback(hObject, eventdata, handles)
maxp1a = get(handles.max_p1,'String');
maxp1 = str2double(maxp1a);

% --- Executes during object creation, after setting all properties.
function max_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p2_Callback(hObject, eventdata, handles)
minp2a = get(handles.min_p2,'String');
minp2 = str2double(minp2a);

% --- Executes during object creation, after setting all properties.
function min_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_p2_Callback(hObject, eventdata, handles)
maxp2a = get(handles.max_p2,'String');
maxp2 = str2double(maxp2a);
% --- Executes during object creation, after setting all properties.
function max_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function step_p1_Callback(hObject, eventdata, handles)
stepp1a = get(handles.step_p1,'String');
stepp1 = str2double(stepp1a);

% --- Executes during object creation, after setting all properties.
function step_p1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function step_p2_Callback(hObject, eventdata, handles)
stepp2a = get(handles.step_p2,'String');
stepp2 = str2double(stepp2a);

% --- Executes during object creation, after setting all properties.
function step_p2_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in simparam3.
function simparam3_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function simparam3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function min_p3_Callback(hObject, eventdata, handles)
minp3a = get(handles.min_p3,'String');
minp3 = str2double(minp3a);


% --- Executes during object creation, after setting all properties.
function min_p3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function max_p3_Callback(hObject, eventdata, handles)
maxp3a = get(handles.max_p3,'String');
maxp3 = str2double(maxp3a);

% --- Executes during object creation, after setting all properties.
function max_p3_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in use_3v.
function use_3v_Callback(hObject, eventdata, handles)
on3d=get(handles.use_3v,'Value');

if on3d==0
    set(handles.var3_txt,'Visible','off','Enable','on')
    set(handles.simparam3,'Visible','off','Enable','on','String','');
    set(handles.min_p3,'Visible','off','Enable','on','String',''); set(handles.max_p3,'Visible','off','Enable','on','String','');
elseif on3d==1
    set(handles.var3_txt,'Visible','on','Enable','on')
    set(handles.simparam3,'Visible','on','Enable','on','String',handles.var_names);
    set(handles.min_p3,'Visible','on','Enable','on','String','0'); set(handles.max_p3,'Visible','on','Enable','on','String','1');
end
guidata(hObject,handles)

% --- Executes on button press in parallel.
function parallel_Callback(hObject, eventdata, handles)

% --- Executes on selection change in bifursp.
function bifursp_Callback(hObject, eventdata, handles)
fp=get(handles.bifursp,'Value');
if fp>1
    handles.mpltopt=1; handles.simtype='multiple_p';
    run_script_gui
else
    handles.mpltopt=0;
end
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function bifursp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function undock_Callback(hObject, eventdata, handles)
fig=figure;
ax=axes;
new_h=copyobj(handles.plothandle,ax,'legacy');


% --------------------------------------------------------------------
function Dimensions_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function twodphase_Callback(hObject, eventdata, handles)
if ~isempty(handles.yout_phase)
    reset(gca)
    axes(handles.trajectoryplot)
    str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
    val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
    param1=str1(val1,:); param2=str2(val2,:);

    h=plot(handles.yout_phase,handles.yout_phase2,handles.colr_phase,'linestyle',handles.mrk_phase); hold on
    %Plot initial conditions
    plot(handles.xlns,handles.ylns,'ks','markersize',7,'markerfacecolor','r')
    %Plot steady-states (if available), otherwise plot final values
    plot(handles.yout_phase_end,handles.yout_phase_end2,'ko','markersize',10,'markerfacecolor',handles.colre_phase)
    axis tight
    xlabel(param1); ylabel(param2)
    set(h,'linewidth',2)
    %Save Report
    newfig_j=figure;
    axesobject_j=copyobj(handles.trajectoryplot,newfig_j);
    title('2D Trajectory Phase Plot')
    export_fig temp_fig/ExpandReport -pdf -r600 -append
    close(newfig_j)
end

% --------------------------------------------------------------------
function threedphase_Callback(hObject, eventdata, handles)
if ~isempty(handles.yout_phase)
    reset(gca)
    axes(handles.trajectoryplot)
    str1=get(handles.simparam1,'String'); str2=get(handles.simparam2,'String');
    val1=get(handles.simparam1,'Value'); val2=get(handles.simparam2,'Value');
    param1=str1(val1,:); param2=str2(val2,:);

    h=plot3(handles.yout_phase,handles.yout_phase2,handles.yout_phase3,handles.colr_phase,'linestyle',handles.mrk_phase); hold on
    %Plot initial conditions
    plot3(handles.xlns,handles.ylns,handles.zlns,'ks','markersize',7,'markerfacecolor','r')
    %Plot steady-states (if available), otherwise plot final values
    plot3(handles.yout_phase_end,handles.yout_phase_end2,handles.yout_phase_end3,'ko','markersize',10,'markerfacecolor',handles.colre_phase)
    axis tight
    xlabel(param1); ylabel(param2); zlabel(str1(handles.val3,:));
    set(h,'linewidth',2); grid on
    %Save Report
    newfig_k=figure;
    axesobject_k=copyobj(handles.trajectoryplot,newfig_k);
    title('3D Trajectory Phase Plot')
    export_fig temp_fig/ExpandReport -pdf -r600 -append
    close(newfig_k)
end

% --- Executes on button press in normeig.
function normeig_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function uiundock_ClickedCallback(hObject, eventdata, handles)
fig=figure;
%Check if subplot or normal figure to be copied
val_p=get(handles.subplot_sol,'Checked');
val_pp=get(handles.sim_pp,'Checked');
if strcmp(val_p,'on') && strcmp(val_pp,'off')
    kl=length(handles.plotaxes_sub);
    if kl == 3 || kl == 4
        sa = 2; sb=2;
    elseif kl == 5 || kl == 6
        sa = 3; sb=2;
    end
    for k=1:kl
        spl(k)=subplot(sa,sb,k);
        copyobj(handles.plotaxes_sub{k},spl(k),'legacy');
        ylim=get(gca,'Ylim');
        endv=handles.plotaxes_sub{1}.XData;
        set(gca,'Xlim',[0 endv(end)],'Ylim',[0 ylim(2)]);
        grid on
        legend(handles.hlegend(k).String)

    end
        suplabel('Time (days)');
        suplabel('Concentration (kgCOD m^{-3})','y');
else
    ax=axes;
    new_h=copyobj(handles.plothandle,ax,'legacy');
    grid on
    axis tight
end


% --------------------------------------------------------------------
function report_menu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function gen_rep_Callback(hObject, eventdata, handles)
%Set Check mark%%%%%
%%%%%%%%%%%%%%%%%%%%
%Check for email and server addresses
strem=get(handles.email_add,'String');
strserv=get(handles.server_add,'String');
if strfind(strem,'@') && ~isempty(strserv)
try
   setpref('Internet','E_mail','matthew.wade@newcastle.ac.uk'); 
   setpref('Internet','SMTP_Server',strserv);
   sendmail(strem,'Microbial Ecology System Analysis Report',['Attached is a PDF report of the mathematical analysis of your microbial models created on ',char(datetime)],'temp_fig/ExpandReport.pdf');
   %New name of report
   dte=strcat(datestr(clock,'yyyy-mm-dd-HHMM'));
   newfile=['Reports/ExpandReport_',dte,'.pdf'];
   movefile('temp_fig/ExpandReport.pdf',newfile)
end
elseif strfind(strem,'@')==0 && ~isempty(strserv)
    msgbox('Enter valid e-mail address','Error: No valid e-mail');
elseif isempty(strserv) && strfind(strem,'@')==1
    msgbox('Enter valid server address','Error: No valid e-mail');
else
    msgbox('Enter valid e-mail and server address','Error: No valid e-mail');
end

% --------------------------------------------------------------------
function del_rep_Callback(hObject, eventdata, handles)
%Delete report(s) function here
%Set Check mark%%%%%
%%%%%%%%%%%%%%%%%%%%

function email_add_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function email_add_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function server_add_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function server_add_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in solmeth_select.
function solmeth_select_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function solmeth_select_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in colormp.
function colormp_Callback(hObject, eventdata, handles)
axes(handles.trajectoryplot)
strcmap=get(handles.colormp,'String');
numstrcmap=get(handles.colormp,'Value');
if numstrcmap==1
    numstrcmap==2;
end
cmap = strtrim(strcmap(numstrcmap,:));
colormap(char(cmap));


% --- Executes during object creation, after setting all properties.
function colormp_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
