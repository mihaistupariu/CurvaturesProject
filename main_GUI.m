%% Description 

% This is the GUI for computing/comparing curvatures methods for triangle
% meshes. It calls several functions / scripts. 
% July 2013 - February 2020, Mihai-Sorin Stupariu (all scripts / functions)

%% 
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

% Last Modified by GUIDE v2.5 01-Jul-2020 12:28:03

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
% Noise
handles.boolNoise=0;
handles.levelNoise_x=0;
handles.levelNoise_y=0;
handles.levelNoise_z=0;
handles.levelOfNoise=1;
handles.noiseFromStructure=0;
handles.writeToFile=0;
handles.writeToFileCorr=0;
handles.nrAlphas=0;
handles.chosenMethods=zeros(2,8+handles.nrAlphas);
handles.chosenMethods(2,8+handles.nrAlphas)=1;
handles.boolAlphaGraph=0;
% Graphics single curvature
handles.typeDrawCurvature=1;
handles.indexDrawCurvature=1;
handles.Az=-126; 
handles.Elev=60;
handles.shadingStyle=1;
handles.drawColorbar=0;
handles.drawFigureCurvature=0;
% Graphics correlations
handles.methodCAE='_GC';
handles.descriptor='correlation';
handles.ylabel='Correlation coefficient';
handles.drawFigureCorr_AE=0;
guidata (hObject, handles);
warning('off', 'all');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THE INPUT DATA BLOCK
 
%--------------------------------
% Read the file with the data
%--------------------------------

% --- Executes during object creation, after setting all properties.
function Button_Group_Data_Type_CreateFcn(hObject, eventdata, handles)
function Button_Data_PC_CreateFcn(hObject, eventdata, handles)
function Button_Data_PC_Callback(hObject, eventdata, handles)
function Button_Data_SS_Callback(hObject, eventdata, handles)
function Button_Group_Data_Type_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in Button_Group_Data_Type 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch (get(eventdata.NewValue, 'Tag'))
    case 'Button_Data_PC'
        handles.dataType=1;
    case 'Button_Data_SS'
        handles.dataType=2;
end
guidata (hObject, handles);

function Button_Load_Data_CreateFcn(hObject, eventdata, handles)

% --- Executes on button press in Button_Load_Data.
function Button_Load_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Load_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch (handles.dataType)
    case 1
        set(handles.Console, 'string', 'The point cloud file was read.');
        [filename1,filepath1]=uigetfile({'*.*','All Files'},...
        'Select Data File 1');
        cd(filepath1);      
        handles.file_data=load([filepath1 filename1]);
    case 2
        set(handles.Console, 'string', 'The file with PD functions was read.');
        [filename1,filepath1]=uigetfile({'*.*','All Files'},...
        'Select Data File 1');
        cd(filepath1);      
        handles.file_functions=fopen([filepath1 filename1]);
end
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
aux_levelResol=transpose(sort(levelResol_aux));
if (handles.dataType==2)
    handles.levelResol=flip(aux_levelResol);
else
    handles.levelResol=aux_levelResol;
end

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

%% NOISE

%--------------------------------
% x_Noise
%--------------------------------
function Values_Noise_x_Callback(hObject, eventdata, handles)
if isempty(get(hObject,'String'))
    levelNoise_x=0;
else 
    levelNoise_x_aux_char=textscan(get(hObject,'String'),'%f');
    levelNoise_x_aux=cell2mat(levelNoise_x_aux_char);
    levelNoise_x=flip(transpose(sort(levelNoise_x_aux)));
    handles.boolNoise=1;
end 
handles.levelNoise_x=levelNoise_x;
handles.noiseFromStructure=0;
set(handles.Console, 'string', ['The values for x-amplitude are ', mat2str(handles.levelNoise_x)]);
guidata(hObject, handles);

function Values_Noise_x_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% y_Noise
%--------------------------------
function Values_Noise_y_Callback(hObject, eventdata, handles)
if isempty(get(hObject,'String'))
    levelNoise_y=0;
else 
    levelNoise_y_aux_char=textscan(get(hObject,'String'),'%f');
    levelNoise_y_aux=cell2mat(levelNoise_y_aux_char);
    levelNoise_y=flip(transpose(sort(levelNoise_y_aux)));
    handles.boolNoise=1;
end   
handles.levelNoise_y=levelNoise_y;
handles.noiseFromStructure=0;
set(handles.Console, 'string', ['The values for y-amplitude are ', mat2str(handles.levelNoise_y)]);
guidata(hObject, handles);

function Values_Noise_y_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% z-Noise
%--------------------------------
function Values_Noise_z_Callback(hObject, eventdata, handles)
if isempty(get(hObject,'String'))
    levelNoise_z=0;
else 
    levelNoise_z_aux_char=textscan(get(hObject,'String'),'%f');
    levelNoise_z_aux=cell2mat(levelNoise_z_aux_char);
    levelNoise_z=flip(transpose(sort(levelNoise_z_aux)));
    handles.boolNoise=1;
end
handles.levelNoise_z=levelNoise_z;
handles.noiseFromStructure=0;
set(handles.Console, 'string', ['The values for z-amplitude are ', mat2str(handles.levelNoise_z)]);
guidata(hObject, handles);


function Values_Noise_z_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Noise read from workspace as a structure
%--------------------------------
function Button_Load_Noise_Callback(hObject, eventdata, handles)
handles.noiseStructure=evalin('base', 'noiseStructure');
handles.boolNoise=1;
handles.noiseFromStructure=1;
set(handles.Console, 'string', 'Noise structure was read.');
guidata(hObject, handles);

function Button_Load_Noise_CreateFcn(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THE MAIN COMPUTATIONS

%% 1. The case of a point cloud as input 
if (handles.dataType==1)   
    % Write console message
    set(handles.Console, 'string', 'Computing ...... ');
    guidata(hObject, handles);
    % Read data
    pcInput=handles.file_data;
    % Limits
    xMin=handles.Value_xMin;
    xMax=handles.Value_xMax;
    yMin=handles.Value_yMin;
    yMax=handles.Value_yMax;
    if (handles.implicit_limits==1)
        xMin=min(pcInput(:,1));
        xMax=max(pcInput(:,1));
        yMin=min(pcInput(:,2));
        yMax=max(pcInput(:,2));
    end
    alphaVector=handles.alphaVector;
    % Levels of Resolution
    levelsResolution=handles.levelResol;
    

    
    % Prepare the point cloud
    pcAuxiliary=generateAuxPointCloud(pcInput, xMin, xMax, yMin, yMax);
    
    % GenerateTriangleMesh and compute curvatures for each level of resolution
    [~,nrLevelRes]=size(levelsResolution);
    for ll=1:nrLevelRes
        cellSize=levelsResolution(ll); % set the size of the cell
        [pcFinal, regGrid]= generateFinalPointCloud (pcAuxiliary, cellSize, xMin, xMax, yMin, yMax); % generate the point cloud through interpolation   
        
        % Write variables to workspace
        level=num2str(ll);
        assignin('base', ['PointCloud_' level], pcFinal);
        assignin('base', ['Reg_Grid_', level], regGrid);  
        
        % COMPUTE CURVATURES
        [GC,MC,vfConvexHull]=computeCurvatures_PC(pcFinal, alphaVector,cellSize, regGrid);
        
        % Add curvatures/grid/convex hull to structure 
        label=['level' level];
        structureGC.(matlab.lang.makeValidName(label))=GC;
        structureMC.(matlab.lang.makeValidName(label))=MC;
        structureRegGrid.(matlab.lang.makeValidName(label))=regGrid;
        structureVfConvexHull.(matlab.lang.makeValidName(label))=vfConvexHull;
        
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
    
    % Initialize values for subsequent steps
    handles.GC=structureGC;
    handles.MC=structureMC;
    handles.RegGrid=structureRegGrid;
    handles.VfConvexHull=structureVfConvexHull;
    handles.chosenMethods=zeros(2,8+handles.nrAlphas);
    handles.chosenMethods(2,8+handles.nrAlphas)=1;
    handles.aux=zeros(1,8+handles.nrAlphas);
    set(handles.Console, 'string', 'Curvatures computed!');
    guidata(hObject, handles);
end
 
%% 2. The case of a smooth surface 
if (handles.dataType==2)   
 
    % Limits
    xMin=handles.Value_xMin;
    xMax=handles.Value_xMax;
    yMin=handles.Value_yMin;
    yMax=handles.Value_yMax;
    alphaVector=handles.alphaVector;
   % Write console message
    set(handles.Console, 'string', 'Computing ...... ');
    guidata(hObject, handles);
    
    % Levels of Resolution
    levelsResolution=handles.levelResol;
    
    % Write console message
    set(handles.Console, 'string', 'Computing ...... ');
    guidata(hObject, handles);
    
    % Read levels of resolution
    [~,nrLevelRes]=size(levelsResolution);
    
        [~, nc_x]=size(handles.levelNoise_x);
        [~, nc_y]=size(handles.levelNoise_y);
        [~, nc_z]=size(handles.levelNoise_z);
        if nc_x==1 && nc_y==1 && nc_z==1 && handles.levelNoise_x(1,1)==0 &&  handles.levelNoise_y(1,1)==0 &&  handles.levelNoise_z(1,1)==0  && handles.noiseFromStructure==0
            handles.boolNoise=0;
        else
            handles.boolNoise=1;
        end
    
    if handles.boolNoise==0
      % Case 1. Without noise. GenerateTriangleMesh and compute curvatures for each level of resolution
        for ll=1:nrLevelRes
            cellSize=levelsResolution(ll); % set the size of the cell
            [X,Y]=meshgrid(xMin:cellSize:xMax,yMin:cellSize:yMax);  % generate 2D mesh
            fid=handles.file_functions;
            frewind(fid);
            line1=fgetl(fid);
            func1=str2func(line1);
            regGrid=func1(X,Y); 
            pcFinal=[];
            pcFinal(:,1)=makeVector(X);
            pcFinal(:,2)=makeVector(Y);
            pcFinal(:,3)=makeVector(regGrid);
            % Write variables to workspace
            level=num2str(ll);
            assignin('base', ['PointCloud_' level], pcFinal);
            assignin('base', ['Reg_Grid_', level], regGrid); 

            % COMPUTE CURVATURES
            [GC,MC,vfConvexHull]=computeCurvatures_SS(pcFinal, alphaVector,cellSize, handles.file_functions);

            % Add curvatures/grid/convex hull to structure 
            label=['level' level];
            structureGC.(matlab.lang.makeValidName(label))=GC;
            structureMC.(matlab.lang.makeValidName(label))=MC;
            structureRegGrid.(matlab.lang.makeValidName(label))=regGrid;
            structureVfConvexHull.(matlab.lang.makeValidName(label))=vfConvexHull;

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
    else
        % Case 2. With noise
        if handles.noiseFromStructure==0
            cellSize=levelsResolution(1); % set the size of the cell; take only the one for the first level of resolution
            [Xorig,Yorig]=meshgrid(xMin:cellSize:xMax,yMin:cellSize:yMax);  % generate 2D mesh, it will be altered for each level of noise
            % GENERATE UNIFORM LEVELS OF AMPLITUDE FOR NOISE
            % Find the largest vector of amplitudes between x, y, z and extend
            % the others with 0 ==> vectors of same size for x, y, z
            [~, nc_x]=size(handles.levelNoise_x);
            [~, nc_y]=size(handles.levelNoise_y);
            [~, nc_z]=size(handles.levelNoise_z);
            maxNoiseLevels=max(nc_x, max( nc_y, nc_z));
            compl_x=zeros(1, maxNoiseLevels-nc_x);
            levelNoise_x=[handles.levelNoise_x compl_x];
            compl_y=zeros(1, maxNoiseLevels-nc_y);
            levelNoise_y=[handles.levelNoise_y compl_y];
            compl_z=zeros(1, maxNoiseLevels-nc_z);
            levelNoise_z=[handles.levelNoise_z compl_z];
            % As handles for subsequent steps
            handles.levelNoise_x=levelNoise_x;
            handles.levelNoise_y=levelNoise_y;
            handles.levelNoise_z=levelNoise_z;
            handles.levelNoise_xyz=[levelNoise_x;levelNoise_y;levelNoise_z];
            %% For each level of amplitude 
            for noiseLevel=1:maxNoiseLevels
                    % GENERATE THE 'NOISY' POINT CLOUD
                    % Prepare for generating random matrices
                    zClock=clock; new_gen=floor(zClock(1,6)); rng(new_gen);
                    [rows_x,cols_x]=size(Xorig);
                    % Generate x,y-noise by perturbation of cellSize
                    matNoise_X=zeros(rows_x,cols_x);
                    matNoise_Y=zeros(rows_x,cols_x);
                    matNoise_X(2:rows_x-1,2:cols_x-1)=cellSize*levelNoise_x(noiseLevel)*(2*rand(rows_x-2,cols_x-2)-1);
                    matNoise_Y(2:rows_x-1,2:cols_x-1)=cellSize*levelNoise_y(noiseLevel)*(2*rand(rows_x-2,cols_x-2)-1);
                    % Add noise to original X Y
                    X=Xorig+matNoise_X;
                    Y=Yorig+matNoise_Y;
                    % Prepare generation of z-values
                    fid=handles.file_functions;
                    frewind(fid);
                    line1=fgetl(fid);
                    func1=str2func(line1);
                    regGridOrig=func1(X,Y); 
                    % Generate z-noise by perturbation of the z-value
                    matNoise_Z=zeros(rows_x,cols_x);
                    matNoise_Z(2:rows_x-1,2:cols_x-1)=levelNoise_z(noiseLevel)*(2*rand(rows_x-2,cols_x-2)-1);
                    matNoise_Z=regGridOrig.*matNoise_Z;
                    regGrid=regGridOrig+matNoise_Z;
                    pcFinal(:,1)=makeVector(X);
                    pcFinal(:,2)=makeVector(Y);
                    pcFinal(:,3)=makeVector(regGrid);

                    % Write variables to workspace
                    strNoiseLevel=num2str(noiseLevel);
                    assignin('base', ['PointCloud_' strNoiseLevel], pcFinal);
                    assignin('base', ['Reg_Grid_', strNoiseLevel], regGrid); 

                    % COMPUTE CURVATURES
                    [GC,MC,vfConvexHull]=computeCurvatures_SS(pcFinal, alphaVector,cellSize, handles.file_functions);

                     % Add curvatures/grid/convex hull to structure 
                    label=['levelNoise' strNoiseLevel];
                    structureGC.(matlab.lang.makeValidName(label))=GC;
                    structureMC.(matlab.lang.makeValidName(label))=MC;
                    structureRegGrid.(matlab.lang.makeValidName(label))=regGrid;
                    structureVfConvexHull.(matlab.lang.makeValidName(label))=vfConvexHull;

                    % Write curvatures to workspace
                    assignin('base', ['GC_' strNoiseLevel], GC);
                    assignin('base', ['MC_' strNoiseLevel], MC); 
                    assignin('base', ['vfConvexHull_' strNoiseLevel], vfConvexHull);


                     % Write to file if required
                    if handles.writeToFile==1
                        fileID = fopen(['GC' strNoiseLevel '.txt'], 'w');
                        dlmwrite(['GC' strNoiseLevel '.txt'], GC, 'delimiter',' ','precision','% 8.2f')
                        fclose(fileID);

                        fileID = fopen(['MC' strNoiseLevel '.txt'], 'w');
                        dlmwrite(['MC' strNoiseLevel '.txt'], MC, 'delimiter',' ','precision','% 8.2f')
                        fclose(fileID);
                    end  
                    %{
                    assignin('base', ['X_' strNoiseLevel],matNoise_X);
                    assignin('base', ['Y_' strNoiseLevel], matNoise_Y);
                    assignin('base', ['regGridOrig_' strNoiseLevel], matNoise_Z);
                    assignin('base', ['regGrid_' strNoiseLevel], regGrid);
                    %}

            end 
        end
        %% If matrix of noise provided
        if handles.noiseFromStructure==1
            noiseStructure=handles.noiseStructure; % extract the noise structure
            maxNoiseLevels=noiseStructure.size; % read number of levels of noise
            cellSize=levelsResolution(1); % set the size of the cell; take only the one for the first level of resolution
            [Xorig,Yorig]=meshgrid(xMin:cellSize:xMax,yMin:cellSize:yMax);  % generate 2D mesh, it will be altered for each level of noise
            %% For each level of amplitude 
            for noiseLevel=1:maxNoiseLevels
                    % Extract the noise
                    matNoise_X=noiseStructure.lev(noiseLevel).x;
                    matNoise_Y=noiseStructure.lev(noiseLevel).y;
                    % Add noise to original X Y
                    X=Xorig+matNoise_X;
                    Y=Yorig+matNoise_Y;
                    % Prepare generation of z-values
                    fid=handles.file_functions;
                    frewind(fid);
                    line1=fgetl(fid);
                    func1=str2func(line1);
                    regGridOrig=func1(X,Y); 
                    % Generate z-noise by perturbation of the z-value
                    matNoise_Z=noiseStructure.lev(noiseLevel).z;
                    matNoise_Z=regGridOrig.*matNoise_Z;
                    regGrid=regGridOrig+matNoise_Z;
                    pcFinal(:,1)=makeVector(X);
                    pcFinal(:,2)=makeVector(Y);
                    pcFinal(:,3)=makeVector(regGrid);

                    % Write variables to workspace
                    strNoiseLevel=num2str(noiseLevel);
                    assignin('base', ['PointCloud_' strNoiseLevel], pcFinal);
                    assignin('base', ['Reg_Grid_', strNoiseLevel], regGrid); 

                    % COMPUTE CURVATURES
                    [GC,MC,vfConvexHull]=computeCurvatures_SS(pcFinal, alphaVector,cellSize, handles.file_functions);

                     % Add curvatures/grid/convex hull to structure 
                    label=['levelNoise' strNoiseLevel];
                    structureGC.(matlab.lang.makeValidName(label))=GC;
                    structureMC.(matlab.lang.makeValidName(label))=MC;
                    structureRegGrid.(matlab.lang.makeValidName(label))=regGrid;
                    structureVfConvexHull.(matlab.lang.makeValidName(label))=vfConvexHull;

                    % Write curvatures to workspace
                    assignin('base', ['GC_' strNoiseLevel], GC);
                    assignin('base', ['MC_' strNoiseLevel], MC); 
                    assignin('base', ['vfConvexHull_' strNoiseLevel], vfConvexHull);


                     % Write to file if required
                    if handles.writeToFile==1
                        fileID = fopen(['GC' strNoiseLevel '.txt'], 'w');
                        dlmwrite(['GC' strNoiseLevel '.txt'], GC, 'delimiter',' ','precision','% 8.2f')
                        fclose(fileID);

                        fileID = fopen(['MC' strNoiseLevel '.txt'], 'w');
                        dlmwrite(['MC' strNoiseLevel '.txt'], MC, 'delimiter',' ','precision','% 8.2f')
                        fclose(fileID);
                    end  
                    %{
                    assignin('base', ['X_' strNoiseLevel],matNoise_X);
                    assignin('base', ['Y_' strNoiseLevel], matNoise_Y);
                    assignin('base', ['regGridOrig_' strNoiseLevel], matNoise_Z);
                    assignin('base', ['regGrid_' strNoiseLevel], regGrid);
                    %}

            end 
        end
    end
 
    % Initialize values for subsequent steps
    handles.GC=structureGC;
    handles.MC=structureMC;
    handles.RegGrid=structureRegGrid;
    handles.VfConvexHull=structureVfConvexHull;
    handles.chosenMethods=zeros(2,8+handles.nrAlphas);
    handles.chosenMethods(2,8+handles.nrAlphas)=1;
    handles.aux=zeros(1,8+handles.nrAlphas);
    set(handles.Console, 'string', 'Curvatures computed!');
    guidata(hObject, handles);
 
end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% THE CORRELATIONS BLOCK

%--------------------------------
% Select curvatures to be considered
%--------------------------------
function Check_GB1_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,1)=value;
guidata(hObject, handles);
function Check_GB2_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,2)=value;
guidata(hObject, handles);
function Check_ET_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,3)=value;
guidata(hObject, handles);
function Check_TA_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,4)=value;
guidata(hObject, handles);
function Check_JF_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,5)=value;
guidata(hObject, handles);
function Check_NC_1R_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,6)=value;
guidata(hObject, handles);
function Check_NC_2R_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,7)=value;
guidata(hObject, handles);

%-------------
% The alpha-values(IA)
%-------------
function Check_Select_alpha_Callback(hObject, eventdata, handles)
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

function Check_Select_alpha_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Check_RG_SS_Callback(hObject, eventdata, handles)
value=get(hObject, 'Value');
handles.chosenMethods(1,8+handles.nrAlphas)=value;
guidata(hObject, handles);

%--------------------------------
% Select the reference curvature (correlations)
%--------------------------------
function ButtonGroup_Basis_CreateFcn(hObject, eventdata, handles)
function ButtonGroup_Basis_SelectionChangedFcn(hObject, eventdata, handles)
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


function Reference_Select_alpha_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%--------------------------------
% Option to write to file correlations
%--------------------------------
function WriteToFileCorr_SelectionChangedFcn(hObject, eventdata, handles)
switch (get(eventdata.NewValue, 'Tag'))
    case 'CorrYes'
        handles.writeToFileCorr=1;
        set(handles.Console, 'string', 'Write correlations to files');
    case 'CorrNo'
        handles.writeToFileCorr=0;
        set(handles.Console, 'string', 'Does not write correlations to files');
end
guidata (hObject, handles);
function WriteToFileCorr_CreateFcn(hObject, eventdata, handles)

%--------------
% Button - compute
%--------------
function Button_Correlations_Callback(hObject, eventdata, handles)
for ii=1:8+handles.nrAlphas
    if handles.chosenMethods(2,ii)==1
        methodReference=ii;
    end  
end

structureGC=handles.GC;
structureMC=handles.MC;
structureRegGrid=handles.RegGrid;
structureVfConvexHull=handles.VfConvexHull;
if handles.boolNoise==0
    % if no noise, the levels are the levels of resolution
    levelsResol=handles.levelResol;
    [~,nrLevels]=size(levelsResol);
    for ll=1:nrLevels % take each level
        level_char=['level' num2str(ll)];
        GC=structureGC.(matlab.lang.makeValidName(level_char));
        MC=structureMC.(matlab.lang.makeValidName(level_char));
        regGrid=structureRegGrid.(matlab.lang.makeValidName(level_char));
        vfConvexHull=structureVfConvexHull.(matlab.lang.makeValidName(level_char));

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
else
    if handles.noiseFromStructure==0 % noise from amplitude
        % if noise, the levels correspond to the noise amplitudes
        [~,totNoiseLevels]=size(handles.levelNoise_x);
        for ll=1:totNoiseLevels % take each level
            noise_char=['levelNoise' num2str(ll)];
            GC=structureGC.(matlab.lang.makeValidName(noise_char));
            MC=structureMC.(matlab.lang.makeValidName(noise_char));
            regGrid=structureRegGrid.(matlab.lang.makeValidName(noise_char));
            vfConvexHull=structureVfConvexHull.(matlab.lang.makeValidName(noise_char));

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
    else   % if noise from structure
        % ithe levels correspond to the noise amplitudes
        totNoiseLevels=handles.noiseStructure.size;
        for ll=1:totNoiseLevels % take each level
            noise_char=['levelNoise' num2str(ll)];
            GC=structureGC.(matlab.lang.makeValidName(noise_char));
            MC=structureMC.(matlab.lang.makeValidName(noise_char));
            regGrid=structureRegGrid.(matlab.lang.makeValidName(noise_char));
            vfConvexHull=structureVfConvexHull.(matlab.lang.makeValidName(noise_char));

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
 
           
    end
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
function Button_Correlations_ButtonDownFcn(hObject, eventdata, handles)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRAPHICS BLOCK 

%--------------------------------
% SINGLE CURVATURE
%--------------------------------

%--------------------------------
% Create axes
%--------------------------------
function axes2_ButtonDownFcn(hObject, eventdata, handles)

%--------------------------------
% Select curvature type
%--------------------------------
function PopMenu_CurvatureType_Callback(hObject, eventdata, handles)
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Gaussian curvature'
        handles.typeDrawCurvature=1;
    case 'Mean curvature'
        handles.typeDrawCurvature=2;
end       
guidata(hObject, handles);        
%
function PopMenu_CurvatureType_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select curvature method
%--------------------------------
function PopMenu_Method_Callback(hObject, eventdata, handles)
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
%
function PopMenu_Method_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select alpha-values
%--------------------------------
function Value_Alpha_Graphics_Callback(hObject, eventdata, handles)
read_alpha_graph=get(hObject,'String');
if isempty (read_alpha_graph)   
else
    alphaValueGraph_char=textscan(read_alpha_graph,'%f');
    alphaValueGraph=cell2mat(alphaValueGraph_char);
    alphaVector=handles.alphaVector;
    [~,nrAlphas]=size(alphaVector);
    for ii=1:nrAlphas
        if alphaVector(1,ii)==alphaValueGraph
            handles.alphaIndexGraph=ii;
            handles.boolAlphaGraph=1;
        end
    end
end
guidata(hObject, handles);
%
function Value_Alpha_Graphics_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select level of resolution
%--------------------------------
function Value_Resolution_Graphics_Callback(hObject, eventdata, handles)
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
handles.level_resol=level;
guidata(hObject, handles);
%
function Value_Resolution_Graphics_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select level of noise
%--------------------------------
function Value_Noise_Graphics_Callback(hObject, eventdata, handles)
read_levelOfNoise=get(hObject,'String');
if isempty(read_levelOfNoise)
    handles.levelOfNoise=1;
else
    levelOfNoise_aux_char=textscan(get(hObject,'String'),'%f');
    handles.levelOfNoise=cell2mat(levelOfNoise_aux_char);
end
guidata(hObject, handles);
function Value_Noise_Graphics_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%--------------------------------
% Parameters for 3D view
%--------------------------------
function Value_Az_Callback(hObject, eventdata, handles)
read_az=get(hObject,'String');
if isempty(read_az)
    handles.Az=-126;
else
    az_aux_char=textscan(get(hObject,'String'),'%f');
    handles.Az=cell2mat(az_aux_char);
end
guidata(hObject, handles);
%
function Value_Az_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
function Value_Elev_Callback(hObject, eventdata, handles)
read_elev=get(hObject,'String');
if isempty(read_elev)
    handles.Elev=60;
else
    elev_aux_char=textscan(get(hObject,'String'),'%f');
    handles.Elev=cell2mat(elev_aux_char);
end
guidata(hObject, handles);
%
function Value_Elev_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select shading style
%--------------------------------
function Popup_Menu_Shading_Callback(hObject, eventdata, handles)
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Faceted'
        handles.shadingStyle=1;
    case 'Flat'
        handles.shadingStyle=2;
    case 'Interp'
        handles.shadingStyle=3;
end
guidata(hObject, handles);  
%
function Popup_Menu_Shading_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Colorbar
%--------------------------------
function DrawColorbarButton_CreateFcn(hObject, eventdata, handles)
%
function DrawColorbarButton_SelectionChangedFcn(hObject, eventdata, handles)
switch (get(eventdata.NewValue, 'Tag'))
    case 'Yes'
        handles.drawColorbar=1;
    case 'No'
        handles.drawColorbar=0;
end
guidata (hObject, handles);

%--------------------------------
% Figure separately
%--------------------------------
function GenerateFigureCurvature_CreateFcn(hObject, eventdata, handles)
%
function GenerateFigureCurvature_SelectionChangedFcn(hObject, eventdata, handles)
switch (get(eventdata.NewValue, 'Tag'))
    case 'Yes'
        handles.drawFigureCurvature=1;
    case 'No'
        handles.drawFigureCurvature=0;
end
guidata (hObject, handles);

%--------------------------------
% Draw button
%--------------------------------
function Button_Draw_Curvature_CreateFcn(hObject, eventdata, handles)
%
function Button_Draw_Curvature_Callback(hObject, eventdata, handles)
generateGraphicsCurvature (handles.RegGrid, handles.GC, handles.MC, handles.typeDrawCurvature,  handles.indexDrawCurvature, ...
    handles.level_resol, handles.boolNoise, handles.levelOfNoise, handles.Az, handles.Elev, handles.shadingStyle, handles.drawColorbar,  handles.drawFigureCurvature);
guidata(hObject, handles);   

%% 
%--------------------------------
% GRAPHICS FOR CORRELATIONS / ABSOLUTE ERROR 
%--------------------------------

%--------------------------------
% Select curvature type
%--------------------------------
function PopMenu_CorrAE_Curvature_Type_Callback(hObject, eventdata, handles)
str=cellstr(get(hObject, 'String'));
val=get(hObject, 'Value');
switch str{val}
    case 'Gaussian curvature'
        handles.methodCAE='_GC';
    case 'Mean curvature'
        handles.methodCAE='_MC';
end
guidata(hObject, handles);
function PopMenu_CorrAE_Curvature_Type_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Select descriptor
%--------------------------------
function PopMenu_Corr_AE_Callback(hObject, eventdata, handles)
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
 
guidata(hObject, handles);
function PopMenu_Corr_AE_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------
% Generate figure
%--------------------------------
function Generate_Figure_Corr_AE_CreateFcn(hObject, eventdata, handles)
function Generate_Figure_Corr_AE_SelectionChangedFcn(hObject, eventdata, handles)
switch (get(eventdata.NewValue, 'Tag'))
    case 'Yes'
        handles.drawFigureCorr_AE=1;
    case 'No'
        handles.drawFigureCorr_AE=0;
end
guidata (hObject, handles);

%--------------------------------
% Draw button
%--------------------------------
function Button_Draw_Corr_AE_CreateFcn(hObject, eventdata, handles)
function Button_Draw_Corr_AE_Callback(hObject, eventdata, handles)
nume=['handles.' handles.descriptor handles.methodCAE];
variabila=eval(nume);

if handles.boolNoise==0
    generateGraphicsCorrelations(handles.boolNoise, variabila, handles.levelResol, handles.chosenMethods, handles.alphaVector,handles.drawFigureCorr_AE);    
else 
    if handles.noiseFromStructure==0
        % Add the (x,y,z)-amplitudes as characters
        [~,nrLevelsNoiseCorrGr]=size(handles.levelNoise_x);
        for ii=1:nrLevelsNoiseCorrGr
           strLevelNoise(1,ii)={['(' num2str(handles.levelNoise_x(1,ii)) ',' num2str(handles.levelNoise_y(1,ii)) ',' num2str(handles.levelNoise_z(1,ii))  ')']};
        end
    else
        % The i-th level corresponds to i
        for ii=1:handles.noiseStructure.size
           strLevelNoise(1,ii)={[ num2str(ii) ]};
        end
    end
    generateGraphicsCorrelations(handles.boolNoise, variabila, strLevelNoise, handles.chosenMethods, handles.alphaVector,handles.drawFigureCorr_AE); 
end

ylabel(handles.ylabel); 

assignin('base', 'AAMEMORY', handles);
guidata (hObject, handles);
%% VARIA
function Text_Basis_CreateFcn(hObject, eventdata, handles)
function popupmenu10_Callback(hObject, eventdata, handles)
function popupmenu10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
