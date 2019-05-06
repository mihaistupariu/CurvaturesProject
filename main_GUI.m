%% Description 
% This is the GUI for computing/comparing curvatures methods for triangle
% meshes. It calls several functions / scripts. 
% July 2013 - May 2019, Mihai-Sorin Stupariu (all scripts / functions)

function varargout = main_GUI(varargin)


%
% MAIN_GUI MATLAB code for main_GUI.fig
%      MAIN_GUI, by itself, creates a new MAIN_GUI or raises the existing
%      singleton*.
%
%      H = MAIN_GUI returns the handle to a new MAIN_GUI or the handle to
%      the existing singleton*.
%
%      MAIN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GUI.M with the given input arguments.
%
%      MAIN_GUI('Property','Value',...) creates a new MAIN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_GUI

% Last Modified by GUIDE v2.5 06-May-2019 15:02:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @main_GUI_OutputFcn, ...
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


% --- Executes just before main_GUI is made visible.
function main_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_GUI (see VARARGIN)

% Choose default command line output for main_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%% CONSOLE, INITIALIZATIONS & MORE
 
function Console_Callback(hObject, eventdata, handles)
% hObject    handle to Console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Console as text
%        str2double(get(hObject,'String')) returns contents of Console as a double


% --- Executes during object creation, after setting all properties.
function Console_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Console (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Initialize some variables used later
handles.dataType=1;
handles.alphaVector=[];
handles.alphaIndexGraph=0;
handles.writeToFile=0;
handles.writeToFileCorr=0;
handles.nrAlphas=0;
handles.chosenMethods=zeros(2,8+handles.nrAlphas);
handles.chosenMethods(2,8+handles.nrAlphas)=1;
handles.boolAlphaGraph=0;
handles.descriptor='correlation';
handles.ylabel='Correlation coefficient';
handles.methodCAE='_GC';
handles.Az=-159; 
handles.Elev=12;
guidata (hObject, handles);
 
%% THE INPUT DATA BLOCK

%--------------------------------
% Read the file with the data
%--------------------------------

% --- Executes on button press in Button_LoadData.
function Button_LoadData_Callback(hObject, eventdata, handles)
% hObject    handle to Button_LoadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename1,filepath1]=uigetfile({'*.*','All Files'},...
  'Select Data File 1');
  cd(filepath1);
handles.file_data=load([filepath1 filename1]);

set(handles.Console, 'string', 'The input file was read.');
guidata (hObject, handles);

%--------------------------------
% Read the limits (x, y)
%--------------------------------

function Value_xMin_Callback(hObject, eventdata, handles)
% hObject    handle to Value_xMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_xMin as text
%        str2double(get(hObject,'String')) returns contents of Value_xMin as a double

handles.implicit_limits=0;
handles.Value_xMin=str2double(get(hObject,'String'));
set(handles.Console, 'string', ['xMin=' get(hObject,'String')]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_xMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_xMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Value_xMax_Callback(hObject, eventdata, handles)
% hObject    handle to Value_xMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_xMax as text
%        str2double(get(hObject,'String')) returns contents of Value_xMax as a double
handles.Value_xMax=str2double(get(hObject,'String'));
set(handles.Console, 'string', ['xMax=' get(hObject,'String')]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_xMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_xMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Value_yMin_Callback(hObject, eventdata, handles)
% hObject    handle to Value_yMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_yMin as text
%        str2double(get(hObject,'String')) returns contents of Value_yMin as a double
handles.Value_yMin=str2double(get(hObject,'String'));
set(handles.Console, 'string', ['yMin=' get(hObject,'String')]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_yMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_yMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Value_yMax_Callback(hObject, eventdata, handles)
% hObject    handle to Value_yMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_yMax as text
%        str2double(get(hObject,'String')) returns contents of Value_yMax as a double
handles.Value_yMax=str2double(get(hObject,'String'));
set(handles.Console, 'string', ['yMax=' get(hObject,'String')])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_yMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_yMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Choose implicit limits (PC)
%--------------------------------

% --- Executes on button press in Select_Limits_PC.
function Select_Limits_PC_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Limits_PC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Limits_PC
handles.implicit_limits=get(hObject, 'Value');
guidata (hObject, handles);

%% THE PARAMETERS BLOCK
 
%--------------------------------
% 'Levels of resolution' 
%--------------------------------

function Values_Resolutions_Callback(hObject, eventdata, handles)
% hObject    handle to Values_Resolutions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Values_Resolutions as text
%        str2double(get(hObject,'String')) returns contents of Values_Resolutions as a double

levelResol_aux_char=textscan(get(hObject,'String'),'%f');
levelResol_aux=cell2mat(levelResol_aux_char);
handles.levelResol=transpose(sort(levelResol_aux));

set(handles.Console, 'string', ['The levels of resolution are ', mat2str(handles.levelResol)]);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Values_Resolutions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Values_Resolutions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% The alpha values (Integral Approach)
%--------------------------------

function Values_IA_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Values_IA_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Values_IA_alpha as text
%        str2double(get(hObject,'String')) returns contents of Values_IA_alpha as a double
read_alpha_val=get(hObject,'String');

if isempty(read_alpha_val)
    alphaVector_aux=[];
else
    alphaVector_aux_char=textscan(get(hObject,'String'),'%f');
    alphaVector_aux=cell2mat(alphaVector_aux_char);
end
handles.alphaVector=transpose(alphaVector_aux);
[~,handles.nrAlphas]=size(handles.alphaVector);
set(handles.Console, 'string', ['The alpha-values are ', mat2str(handles.alphaVector)]);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Values_IA_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Values_IA_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% THE CURVATURES BLOCK - MAIN COMPUTATIONS

%--------------------------------
% Option to write to file
%--------------------------------

% --- Executes during object creation, after setting all properties.
function WriteToFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WriteToFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes when selected object is changed in WriteToFile.
function WriteToFile_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in WriteToFile 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch (get(eventdata.NewValue, 'Tag'))
    case 'Yes'
        handles.writeToFile=1;
        set(handles.Console, 'string', 'Write curvatures to files');
    case 'No'
        handles.writeToFile=0;
        set(handles.Console, 'string', 'Does not write curvatures to files');
end
guidata (hObject, handles);


%--------------------------------
% Curvature Computations
%--------------------------------

% --- Executes on button press in Button_Curvatures.
function Button_Curvatures_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Curvatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% THE MAIN COMPUTATIONS
    
    pcInput=handles.file_data;
    cellSize=handles.levelResol;
    xMin=handles.Value_xMin;
    xMax=handles.Value_xMax;
    yMin=handles.Value_yMin;
    yMax=handles.Value_yMax;
    alphaVector=handles.alphaVector;
    levelsResolution=handles.levelResol;
    
    if (handles.implicit_limits==1)
        xMin=min(pcInput(:,1));
        xMax=max(pcInput(:,1));
        yMin=min(pcInput(:,2));
        yMax=max(pcInput(:,2));
    end
    
    set(handles.Console, 'string', 'Computing ...... ');
    guidata(hObject, handles);
    
    % Prepare point cloud
    pcAuxiliary=generateAuxPointCloud(pcInput, xMin, xMax, yMin, yMax);
    
    % GenerateTriangleMesh and curvatures at each level of resolution
    [~,nrLevelRes]=size(levelsResolution);
    for ll=1:nrLevelRes
        cellSize=levelsResolution(ll); % set the size of the cell
        [pcFinal, regGrid]= generateFinalPointCloud (pcAuxiliary, cellSize, xMin, xMax, yMin, yMax); % generate the point cloud through interpolation
        
        % Write variables to workspace
        level=num2str(ll);
        assignin('base', ['PointCloud_' level], pcFinal);
        assignin('base', ['Reg_Grid_', level], regGrid);  
        
        % COMPUTE CURVATURES
        [GC,MC,vfConvexHull]=computeCurvatures(pcFinal, alphaVector,cellSize, regGrid);
        
        % Write curvatures to workspace
        assignin('base', ['GC_' level], GC);
        assignin('base', ['MC_' level], MC); 
        assignin('base', ['vfConvexHull_' level], vfConvexHull); 
        
        

        
        % Write to file if required
        if handles.writeToFile==1
            fileID = fopen(['GC' level '.txt'], 'w');
            dlmwrite(['GC' level '.txt'], GC, 'delimiter',' ','precision','% 8.2f')
            fclose(fileID);

            fileID = fopen(['MC' level '.txt'], 'w');
            dlmwrite(['MC' level '.txt'], MC, 'delimiter',' ','precision','% 8.2f')
            fclose(fileID);
        end
        
    end
    
    % Initialize values
    handles.indexDrawCurvature=1;
    handles.typeDrawCurvature=1;
    handles.chosenMethods=zeros(2,8+handles.nrAlphas);
    handles.chosenMethods(2,8+handles.nrAlphas)=1;
    handles.aux=zeros(1,8+handles.nrAlphas);
    set(handles.Console, 'string', 'Curvatures computed!');
    guidata(hObject, handles);
 

    
%% THE CORRELATIONS BLOCK

%--------------------------------
% Select curvatures to be considered
%--------------------------------


% --- Executes on button press in Check_GB1.
function Check_GB1_Callback(hObject, eventdata, handles)
% hObject    handle to Check_GB1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_GB1
value=get(hObject, 'Value');
handles.chosenMethods(1,1)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_GB2.
function Check_GB2_Callback(hObject, eventdata, handles)
% hObject    handle to Check_GB2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_GB2
value=get(hObject, 'Value');
handles.chosenMethods(1,2)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_ET.
function Check_ET_Callback(hObject, eventdata, handles)
% hObject    handle to Check_ET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_ET
value=get(hObject, 'Value');
handles.chosenMethods(1,3)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_TA.
function Check_TA_Callback(hObject, eventdata, handles)
% hObject    handle to Check_TA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_TA
value=get(hObject, 'Value');
handles.chosenMethods(1,4)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_JF.
function Check_JF_Callback(hObject, eventdata, handles)
% hObject    handle to Check_JF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_JF
value=get(hObject, 'Value');
handles.chosenMethods(1,5)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_NC_1R.
function Check_NC_1R_Callback(hObject, eventdata, handles)
% hObject    handle to Check_NC_1R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_NC_1R
value=get(hObject, 'Value');
handles.chosenMethods(1,6)=value;
guidata(hObject, handles);

% --- Executes on button press in Check_NC_2R.
function Check_NC_2R_Callback(hObject, eventdata, handles)
% hObject    handle to Check_NC_2R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_NC_2R
value=get(hObject, 'Value');
handles.chosenMethods(1,7)=value;
guidata(hObject, handles);


%-------------
% The alpha-values(IA)
%-------------
function Check_Select_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Check_Select_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Check_Select_alpha as text
%        str2double(get(hObject,'String')) returns contents of Check_Select_alpha as a double

handles.chosenMethods(1,7+1:7+handles.nrAlphas)=zeros(1,handles.nrAlphas);
read_alpha=get(hObject,'String');
if isempty(read_alpha)
    alphaVector_aux=[];
else
	alphaVector_aux_char=textscan(get(hObject,'String'),'%f');
    alphaVector_aux=cell2mat(alphaVector_aux_char);
end

handles.alphaVectorChosen=transpose(alphaVector_aux);
[~,nrAlphasChosen]=size(handles.alphaVectorChosen);

for ii=1:nrAlphasChosen
    alphaAux=handles.alphaVectorChosen(1,ii); 
    for jj=1:handles.nrAlphas
        if handles.alphaVector(1,jj)==alphaAux
           handles.chosenMethods(1,7+jj)=1;
        end
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Check_Select_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Check_Select_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Check_RG_SS.
function Check_RG_SS_Callback(hObject, eventdata, handles)
% hObject    handle to Check_RG_SS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Check_RG_SS
value=get(hObject, 'Value');
handles.chosenMethods(1,8+handles.nrAlphas)=value;
guidata(hObject, handles);

%--------------------------------
% Select the reference curvature (correlations)
%--------------------------------

% --- Executes during object creation, after setting all properties.
function ButtonGroup_Basis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ButtonGroup_Basis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 
% --- Executes when selected object is changed in ButtonGroup_Basis.
function ButtonGroup_Basis_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ButtonGroup_Basis 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aux(1,:)=zeros(1,8+handles.nrAlphas);
switch (get(eventdata.NewValue, 'Tag'))
    case 'Button_GB1'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,1)=1;
        handles.aux(1,1)=1;
    case 'Button_GB2'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,2)=1;
        handles.aux(1,2)=1;
    case 'Button_ET'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,3)=1;
        handles.aux(1,3)=1;
    case 'Button_TA'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,4)=1;
        handles.aux(1,4)=1;
    case 'Button_JF'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,5)=1;
        handles.aux(1,5)=1;
    case 'Button_NC_1R'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,6)=1;
        handles.aux(1,6)=1;
    case 'Button_NC_2R'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,7)=1;
        handles.aux(1,7)=1;
    case 'Button_RG_SS'
        handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
        handles.chosenMethods(2,8+handles.nrAlphas)=1;
        handles.aux(1,8+handles.nrAlphas)=1;

end
guidata (hObject, handles);

%-------------
% The alpha-values(IA)
%-------------

function Reference_Select_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to Reference_Select_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Reference_Select_alpha as text
%        str2double(get(hObject,'String')) returns contents of Reference_Select_alpha as a double
read_alpha_ref=get(hObject,'String');

handles.chosenMethods(2,:)=handles.aux(1,:);
if isempty(read_alpha_ref)
else 
   alphaReference_char=textscan(get(hObject,'String'),'%f');
   alphaReference=cell2mat(alphaReference_char);
   for jj=1:handles.nrAlphas
       if handles.alphaVector(1,jj)==alphaReference
            handles.chosenMethods(2,1:8+handles.nrAlphas)=zeros(1,8+handles.nrAlphas);
            handles.chosenMethods(2,7+jj)=1;
       end
   end 

    
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Reference_Select_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Reference_Select_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Option to write to file correlations
%--------------------------------

% --- Executes when selected object is changed in WriteToFileCorr.
function WriteToFileCorr_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in WriteToFileCorr 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch (get(eventdata.NewValue, 'Tag'))
    case 'CorrYes'
        handles.writeToFileCorr=1;
        set(handles.Console, 'string', 'Write correlations to files');
    case 'CorrNo'
        handles.writeToFileCorr=0;
        set(handles.Console, 'string', 'Does not write correlations to files');
end
guidata (hObject, handles);


% --- Executes during object creation, after setting all properties.
function WriteToFileCorr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WriteToFileCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%--------------
% Button - compute
%--------------

% --- Executes on button press in Button_Correlations.
function Button_Correlations_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Correlations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 
% Select the reference method
for ii=1:8+handles.nrAlphas
    if handles.chosenMethods(2,ii)==1
        methodReference=ii;
    end  
end

levelsResol=handles.levelResol;
[~,nrLevels]=size(levelsResol);
for ll=1:nrLevels % take each level

    level=num2str(ll);
    GC=evalin('base', ['GC_' level]);
    MC=evalin('base', ['MC_' level]);
    regGrid=evalin('base', ['Reg_Grid_', level]); 
    vfConvexHull=evalin('base', ['vfConvexHull_', level]);
    
    % Set the reference vectors
    vectReference_GC=GC(:,methodReference);
    vectReference_MC=MC(:,methodReference);
    
    % Creates matrices of selected curvatures
    index=0;
    for ii=1:8+handles.nrAlphas
        if handles.chosenMethods(1,ii)==1 || handles.chosenMethods(2,ii)==1
            index=index+1;
            matSelected_GC(:,index)=GC(:,ii);
            matSelected_MC(:,index)=MC(:,ii);
        end
    end
    [correlation_GC(ll,:), absError_GC(ll,:)]=computeCorrelationVector(vectReference_GC, matSelected_GC, vfConvexHull);
    [correlation_MC(ll,:), absError_MC(ll,:)]=computeCorrelationVector(vectReference_MC, matSelected_MC, vfConvexHull);
    clearvars matSelected_GC matSelected_MC
end
% Write correlations to workspace
 assignin('base', 'correlation_GC', correlation_GC);
 assignin('base', 'correlation_MC', correlation_MC);  
 assignin('base', 'absError_GC', absError_GC);
 assignin('base', 'absError_MC', absError_MC); 
 
 % Add correlations to handles
 
 handles.correlation_GC=correlation_GC;
 handles.correlation_MC=correlation_MC;
 handles.absError_GC=absError_GC;
 handles.absError_MC=absError_MC; 
  
  
 % Write to file if required
if handles.writeToFileCorr==1
    filename='correlation_GC.xlsx';
    xlswrite(filename,correlation_GC);
    filename='correlation_MC.xlsx';
    xlswrite(filename,correlation_MC);
    filename='absError_GC.xlsx';
    xlswrite(filename,absError_GC);
    filename='absError_MC.xlsx';
    xlswrite(filename,absError_MC);
end
    set(handles.Console, 'string', 'Correlations computed!');
 
    guidata(hObject, handles);

        


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Button_Correlations.
function Button_Correlations_ButtonDownFcn(hObject, eventdata, handles)



% hObject    handle to Button_Correlations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%% GRAPHICS BLOCK 

%--------------------------------
% SINGLE CURVATURE
%--------------------------------

%--------------------------------
% Create axes
%--------------------------------

% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------
% Select curvature type and method
%--------------------------------

% --- Executes on selection change in PopMenu_CurvatureType.
function PopMenu_CurvatureType_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenu_CurvatureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenu_CurvatureType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenu_CurvatureType
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Gaussian curvature'
        handles.typeDrawCurvature=1;
    case 'Mean curvature'
        handles.typeDrawCurvature=2;
end

if handles.boolAlphaGraph==1;
    handles.indexDrawCurvature=7+handles.alphaIndexGraph;
end
        
if handles.typeDrawCurvature==1;
    handles.vecCol=handles.GC(:,handles.indexDrawCurvature);
end

if handles.typeDrawCurvature==2;
    handles.vecCol=handles.MC(:,handles.indexDrawCurvature);
end
generateGraphicsCurvature (handles.regGrid, handles.vecCol, handles.Az, handles.Elev);
guidata(hObject, handles);        

% --- Executes during object creation, after setting all properties.
function PopMenu_CurvatureType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopMenu_CurvatureType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PopMenu_Method.
function PopMenu_Method_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenu_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenu_Method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenu_Method

handles.boolAlphaGraph=0;
[~,nrAlphas]=size(handles.alphaVector);
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Gauss-Bonnet (I)'
        handles.indexDrawCurvature=1;
    case 'Gauss-Bonnet (II)'
        handles.indexDrawCurvature=2;
    case 'Euler Theorem'
        handles.indexDrawCurvature=3;
    case 'Tensor Approach'
        handles.indexDrawCurvature=4;
    case 'Jet Fitting'
        handles.indexDrawCurvature=5;
    case 'Normal Cycle (1-R)'
        handles.indexDrawCurvature=6;
    case 'Normal Cycle (2-R)'
        handles.indexDrawCurvature=7;
    case 'Integral Approach'
        handles.boolAlphaGraph=1;
        handles.indexDrawCurvature=7+handles.alphaIndexGraph;
    case 'Grid / Smooth'
        handles.indexDrawCurvature=7+nrAlphas+1;
end
guidata(hObject, handles);

if handles.typeDrawCurvature==1;
    handles.vecCol=handles.GC(:,handles.indexDrawCurvature);
end

if handles.typeDrawCurvature==2;
    handles.vecCol=handles.MC(:,handles.indexDrawCurvature);
end
generateGraphicsCurvature (handles.regGrid, handles.vecCol, handles.Az, handles.Elev);
guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function PopMenu_Method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopMenu_Method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Select alpha-values
%--------------------------------

function Value_Alpha_Graphics_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Alpha_Graphics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Alpha_Graphics as text
%        str2double(get(hObject,'String')) returns contents of Value_Alpha_Graphics as a double
alphaValueGraph_char=textscan(get(hObject,'String'),'%f');
alphaValueGraph=cell2mat(alphaValueGraph_char);
alphaVector=handles.alphaVector;
[~,nrAlphas]=size(alphaVector);
for ii=1:nrAlphas
    if alphaVector(1,ii)==alphaValueGraph
        handles.alphaIndexGraph=ii;
        handles.boolAlphaGraph=1;
    end
end
 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Value_Alpha_Graphics_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Alpha_Graphics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Select level of resolution
%--------------------------------

function Value_Resolution_Graphics_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Resolution_Graphics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Resolution_Graphics as text
%        str2double(get(hObject,'String')) returns contents of Value_Resolution_Graphics as a double
levelValueGraph_char=textscan(get(hObject,'String'),'%f');
levelValueGraph=cell2mat(levelValueGraph_char);
levelsResol=handles.levelResol;
[~,nrLevels]=size(levelsResol);
for ii=1:nrLevels
    if levelsResol(1,ii)==levelValueGraph
        ll=ii;
    end
end
 
level=num2str(ll);
handles.regGrid=evalin('base', ['Reg_Grid_' level]);
handles.GC=evalin('base', ['GC_' level]);
handles.MC=evalin('base', ['MC_' level]);


guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function Value_Resolution_Graphics_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Resolution_Graphics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% GRAPHICS FOR CORRELATIONS / ABSOLUTE ERROR 
%--------------------------------

%--------------------------------
% Select curvature type
%--------------------------------

% --- Executes on selection change in PopMenu_CorrAE_Curvature_Type.
function PopMenu_CorrAE_Curvature_Type_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenu_CorrAE_Curvature_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenu_CorrAE_Curvature_Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenu_CorrAE_Curvature_Type
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Gaussian curvature'
        handles.methodCAE='_GC';
    case 'Mean curvature'
        handles.methodCAE='_MC';
end

nume=['handles.' handles.descriptor handles.methodCAE];
variabila=eval(nume);
generateGraphicsCorrelations(variabila,handles.levelResol, handles.chosenMethods, handles.alphaVector);
ylabel(handles.ylabel);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PopMenu_CorrAE_Curvature_Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopMenu_CorrAE_Curvature_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select descriptor
%--------------------------------

% --- Executes on selection change in PopMenu_Corr_AE.
function PopMenu_Corr_AE_Callback(hObject, eventdata, handles)
% hObject    handle to PopMenu_Corr_AE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PopMenu_Corr_AE contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PopMenu_Corr_AE
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Correlations'
        handles.descriptor='correlation';
        handles.ylabel='Correlation coefficient';
    case 'Absolute Errors'
        handles.descriptor='absError';
        handles.ylabel='Absolute error';
end
nume=['handles.' handles.descriptor handles.methodCAE];
variabila=eval(nume);
generateGraphicsCorrelations(variabila,handles.levelResol, handles.chosenMethods, handles.alphaVector);
ylabel(handles.ylabel);
guidata(hObject, handles);
 

% --- Executes during object creation, after setting all properties.
function PopMenu_Corr_AE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PopMenu_Corr_AE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 

%--------------------------------
% Parameters for 3D view
%--------------------------------

function Value_Az_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Az as text
%        str2double(get(hObject,'String')) returns contents of Value_Az as a double
read_az=get(hObject,'String');
if isempty(read_az)
    handles.Az=-159;
else
    az_aux_char=textscan(get(hObject,'String'),'%f');
    handles.Az=cell2mat(az_aux_char);
end
generateGraphicsCurvature (handles.regGrid, handles.vecCol, handles.Az, handles.Elev);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_Az_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Az (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Value_Elev_Callback(hObject, eventdata, handles)
% hObject    handle to Value_Elev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Value_Elev as text
%        str2double(get(hObject,'String')) returns contents of Value_Elev as a double
read_elev=get(hObject,'String');
if isempty(read_elev)
    handles.Elev=12;
else
    elev_aux_char=textscan(get(hObject,'String'),'%f');
    handles.Elev=cell2mat(elev_aux_char);
end
generateGraphicsCurvature (handles.regGrid, handles.vecCol, handles.Az, handles.Elev);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Value_Elev_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Value_Elev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
