function varargout = RICMPostProcessing(varargin)
% RICMPOSTPROCESSING MATLAB code for RICMPostProcessing.fig
%      RICMPOSTPROCESSING, by itself, creates a new RICMPOSTPROCESSING or raises the existing
%      singleton*.
%
%      H = RICMPOSTPROCESSING returns the handle to a new RICMPOSTPROCESSING or the handle to
%      the existing singleton*.
%
%      RICMPOSTPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in. RICMPOSTPROCESSING.M with the given input arguments.
%
%      RICMPOSTPROCESSING('Property','Value',...) creates a new RICMPOSTPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RICMPostProcessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RICMPostProcessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RICMPostProcessing

% Last Modified by GUIDE v2.5 07-Jul-2019 16:13:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RICMPostProcessing_OpeningFcn, ...
                   'gui_OutputFcn',  @RICMPostProcessing_OutputFcn, ...
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


% --- Executes just before RICMPostProcessing is made visible.
function RICMPostProcessing_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RICMPostProcessing (see VARARGIN)

% Choose default command line output for RICMPostProcessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RICMPostProcessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);
[file, path] = uigetfile('*.mat', 'Select MATLAB file that contains CellPrints variable');



openStruct = open(fullfile(path,file));
try CellPrints = openStruct.CellPrints;
catch CellPrints = openStruct.cellPrints;
end

handles.CellPrintBi=CellPrints>0;%convert the map to binary

% handles.currentImage = double(handles.CellPrintBi(:,:,1));
% axes(handles.origVideo)
% imagesc(handles.currentImage);
% colormap gray;
handles.frameNum = 1;
handles.sz=size(handles.CellPrintBi);
handles.l = handles.sz(3);

set(handles.videoSlider,'SliderStep',[1/(handles.l-1) 10/(handles.l-1)]);
handles.origStartingFrameNum = 1;

%%load the time and area file (watch out for the old files that have frames instead of time)
[file,path] = uigetfile('*.txt','Pick time/area file');
CellBatch=file(1:end-4);
set(handles.CellBatch, 'String', CellBatch)
fid = fopen(fullfile(path,file),'r');
timeArea = textscan(fid,'%f %f');
time = timeArea{1}';
handles.refArea= timeArea{2}';
if sum(diff(time) == 1.0) > handles.l/5
    uiwait(msgbox('It looks like the area data has the frames instead of the times. Please select the time/pressures txt file to read in the time values'));
    frames = time;
    [file, path] = uigetfile('*.txt', 'Select times/pressure file');
    fid = fopen(fullfile(path,file),'r');
    timepressure = textscan(fid,'%s %s %s %f %s %s %f','Delimiter','\t');
    allTime = timepressure{4}'; % for all frames from video
    time = allTime(frames); % for only selected frames
end
handles.rawTime = time;
handles.time = time - time(1); 
handles.PxScale = [];
searchDir = dir('*.txt');
txtFilesNames = {searchDir.name};
foundAnalysisInfo = false;
for i = 1:length(txtFilesNames)
    idx = strfind(txtFilesNames{i}, 'analysisInfo');
    if ~isempty(idx)
        if ~foundAnalysisInfo
            foundAnalysisInfo = true;
        else
            uiwait(msgbox('Warning: found multiple analysis info files, using first to find scale value'));
            break
        end
        fid = fopen(txtFilesNames{i},'r');
        c = textscan(fid,'%s','Delimiter','\n');
        textCell = c{1};
        for i = 1:length(textCell)
            curLine = textCell{i};
            scaleIdx = strfind(curLine, 'Scale value: ');
            if ~isempty(scaleIdx)
                scaleTxt = curLine(length('Scale value: ')+1:end);
                unitsIdx = strfind(scaleTxt, 'microns/pixel');
                if ~isempty(unitsIdx)
                    handles.PxScale = str2num(scaleTxt(1:unitsIdx-2));
                    handles.PxScaleUnit = 'µm';
                end
            end
        end
    end
end

while isempty(handles.PxScale)
    choice = questdlg('How do you wish to enter the microns/pixel ratio?','Microns/pixel ratio?','Find existing scale value'...
        ,'Find new scale value','Enter by hand','Find new scale value');
    switch choice
        case 'Find existing scale value'
            [file,path] = uigetfile({'*.mat'},'Choose scale value');
            if file ~= 0
                load(fullfile(path,file), 'scale');
                handles.PxScale = scale;
                handles.PxScaleUnit = 'µm';
            end
        case 'Find new scale value'
            scale = CreateScale; %Use CreateScale function to load image and measure scale
            handles.PxScale = scale;
            [file,path] = uiputfile({'*.mat'},'Save scale value for future use');
            save(fullfile(path,file),'scale')
            handles.PxScaleUnit = 'µm';
        case 'Enter by hand'
            scaleCell = inputdlg({'Enter microns per pixel scale value', 'Enter pixel scale unit (e.g. um or px)'});
            handles.PxScale = str2double(scaleCell{1});
            unitText = scaleCell{2};
            if strcmp(unitText,'um')
                handles.PxScaleUnit = 'µm';
            else
                handles.PxScaleUnit = unitText;
            end
    end
end
% load area data using prints and scale val
% for i = 1:handles.l
%     binIm = handles.CellPrintBi(:,:,i);
%     stats = regionprops(~binIm, 'Area', 'PixelIdxList');
%     areas = [stats(:).Area];
%     handles.Area(i) = max(areas) * handles.PxScale^2;
% end


%set the steps for min and max area sliders to calculate roundness and ellipticity
set(handles.minAreaSlider,'SliderStep',[1/(handles.l-1) 10/(handles.l-1)]);
set(handles.maxAreaSlider,'SliderStep',[1/(handles.l-1) 10/(handles.l-1)]);

handles.minWindowFrame=[];
handles.maxWindowFrame=[];

%initilize for the Calculate
handles.frameNum=1;

guidata(hObject,handles) %make the hObject and handles global variables

% --- Outputs from this function are returned to the command line.
function varargout = RICMPostProcessing_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function videoSlider_Callback(hObject, ~, handles)
% hObject    handle to videoSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

position = get(hObject,'Value');%get the position of the slider

handles.frameNum = round(position * (handles.l-1))+1;
text = sprintf('Frame#%d (%.3fs)', handles.frameNum+handles.origStartingFrameNum - 1, handles.time(handles.frameNum));
set(handles.frameText, 'String', text)

handles.currentImageOverlay = handles.CellPrintRGB(:,:,:,handles.frameNum);
axes(handles.origVideo)
axis square
imshow(handles.currentImageOverlay)

%plot demo
axes(handles.CapturedDemo)
plot(handles.x{handles.frameNum,1},handles.sz(1)-handles.y{handles.frameNum,1},'ro',...
    handles.xIn{handles.frameNum,1},handles.sz(1)-handles.yIn{handles.frameNum,1},'b-',...
    handles.xOut{handles.frameNum,1},handles.sz(1)-handles.yOut{handles.frameNum,1},'b-')
grid on
axis equal
pbaspect([1 1 1])

%add centroid
hold on
plot(handles.centroids(handles.frameNum,1), handles.sz(1)-handles.centroids(handles.frameNum,2), '*')
hold off
          
%add displacement cursor
axes(handles.DisplacementOverTime)
hold on
delete(handles.displacementCursor)
handles.displacementCursor=plot(handles.time(handles.frameNum),handles.displacement(handles.frameNum),'rx');
hold off

%area cursor
axes(handles.AreaOverTime)
hold on
delete(handles.areaCursor)
handles.areaCursor=plot(handles.time(handles.frameNum),handles.Area(handles.frameNum),'rx');
hold off

%roundness cursor
axes(handles.RoundnessOverTime)
hold on
delete(handles.roundnessCursor)
handles.roundnessCursor=plot(handles.time(handles.frameNum),handles.Roundness(handles.frameNum),'rx');
hold off

%update the gui variables
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function videoSlider_CreateFcn(hObject, ~, ~)
% hObject    handle to videoSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function minAreaSlider_Callback(hObject, ~, handles)
% hObject    handle to minAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
position = get(hObject,'Value');%get the position of the slider

%get the start point of roundness, based on minAreaSlider
handles.minWindowFrame = round(position * (handles.l-1))+1;
handles.minArea=handles.Area(handles.minWindowFrame);
handles = updateReadout(handles);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function minAreaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over minAreaSlider.
function minAreaSlider_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to minAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function maxAreaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to maxAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
position = get(hObject,'Value');%get the position of the slider

%get the start point of roundness, based on maxAreaSlider
handles.maxWindowFrame = round(position * (handles.l-1))+1;
handles.maxArea=handles.Area(handles.maxWindowFrame);
handles = updateReadout(handles);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function maxAreaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxAreaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%need the cell#, is it frame or time, min and max area, and relative times,
%average roundness and displacement
inputCell = inputdlg({'Site', 'Trial', 'Cell'}, 'Input cell identifiers');
site = str2double(inputCell{1});
trial = str2double(inputCell{2});
cell = str2double(inputCell{3});
% first, save roundness, ellip vs time data'
[file, path] = uiputfile('*.txt', 'Choose place to save roundness and displacement vs time data');
roundnessFile = fullfile(path,file);
fileID = fopen(roundnessFile, 'w');
fprintf(fileID, '%.4f %.4f %.4f %.4f %.4f\r\n', [handles.rawTime', handles.Roundness, handles.displacement, handles.correctedCentroids(:,1)*handles.PxScale, handles.correctedCentroids(:,2)*handles.PxScale]');
fclose(fileID);
saveas(gcf,fullfile(path,horzcat(file(1:end-4),'.png')));
choice = questdlg('New analysis with max over spreading window or old analysis with multiple averages?', 'Which analysis?', 'New', 'Old', 'New');
switch choice
    case 'New'
        Results=[site, trial, cell, handles.spreadingSpeedData, handles.maxAreaData,...
            handles.spreadingTimeFrame...
            handles.timeFrameMaxRoundness, handles.maxRoundness,...
            handles.timeFrameMaxDisplacement, handles.maxDisplacement,...
            handles.timeFrameMaxOverallDispl, handles.maxOverallDispl, handles.rawTime(end)];
        choice = questdlg('Do you wish to append to a results matrix or start a new file?', 'Append or new?', 'Append', 'Start new file', 'Append');
        switch choice
            case 'Append'
                [file,path] = uigetfile('*.txt', 'Choose results file to append to');
                fileID = fopen(fullfile(path,file),'a');
                fprintf(fileID,'%d %d %d %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\r\n',Results);
                fclose(fileID);
            case 'Start new file'
                [file,path]=uiputfile('*.txt','Save Results');
                fileID = fopen(fullfile(path,file),'w');
                fprintf(fileID,'%d %d %d %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\r\n',Results);
                fclose(fileID);
        end
    case 'Old'
        % next, save summary of this cell with option to append to existing txt
        choice = questdlg('Use automated window for roundness averages?', 'Use automated window?', 'Yes', 'No, use current instead', 'Yes');
        useCurWindow = true;
        switch choice
            case 'Yes'
                if ~isfield(handles, 'maxSlopeIdx')
                    uiwait(msgbox('Automated window was not generated. Using current window selection'))
                else
                    useCurWindow = false;
                end
        end
        if useCurWindow
            centralAvgTime = (handles.time(handles.minWindowFrame) + handles.time(handles.maxWindowFrame)) / 2;
            [~,centralFrame] = min(abs(handles.time - centralAvgTime));
            centralTime = handles.time(centralFrame);
        else
            centralFrame = handles.maxSlopeIdx;
            centralTime = handles.maxSlopeTime;
        end
        [err30Before,frame30Before] = min(abs(handles.time - (centralTime-30)));
        [err15Before,frame15Before] = min(abs(handles.time - (centralTime-15)));
        [err15After,frame15After] = min(abs(handles.time - (centralTime+15)));
        [err30After,frame30After] = min(abs(handles.time - (centralTime+30)));
        range1 = frame30Before:centralFrame;
        timeRange1 = handles.time(centralFrame) - handles.time(frame30Before);
        range2 = frame15Before:frame15After;
        timeRange2 = handles.time(frame15After) - handles.time(frame15Before);
        range3 = centralFrame:frame30After;
        timeRange3 = handles.time(frame30After) - handles.time(centralFrame);
        if err30Before > 5
            uiwait(msgbox('First window could not be calculated'));
            timeFrame1 = [0, 0];
            roundness1 = 0;
        else
            roundness1 = trapz(handles.time(range1), handles.Roundness(range1)) / timeRange1;
            timeFrame1 = [handles.rawTime(frame30Before), handles.rawTime(centralFrame)];
        end
        if err15Before > 5 || err15After > 5
            uiwait((msgbox('Middle window could not be calculated')));
            timeFrame2 = [0, 0];
            roundness2 = 0;
        else
            roundness2 = trapz(handles.time(range2), handles.Roundness(range2)) / timeRange2;
            timeFrame2 = [handles.rawTime(frame15Before), handles.rawTime(frame15After)];
        end
        if err30After > 5
            uiwait(msgbox('Third window could not be calculated'));
            timeFrame3 = [0, 0];
            roundness3 = 0;
        else
            roundness3 = trapz(handles.time(range3), handles.Roundness(range3)) / timeRange3;
            timeFrame3 = [handles.rawTime(centralFrame), handles.rawTime(frame30After)];
        end
        Results=[site, trial, cell, handles.spreadingSpeedData, handles.maxAreaData,...
            timeFrame1, roundness1, timeFrame2, roundness2, timeFrame3, roundness3,...
            handles.timeFrameMaxRoundness, handles.maxRoundness,...
            handles.timeFrameMaxDisplacement, handles.maxDisplacement];
        choice = questdlg('Do you wish to append to a results matrix or start a new file?', 'Append or new?', 'Append', 'Start new file', 'Append');
        switch choice
            case 'Append'
                [file,path] = uigetfile('*.txt', 'Choose results file to append to');
                fileID = fopen(fullfile(path,file),'a');
                fprintf(fileID,'%d %d %d %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\r\n',Results);
                fclose(fileID);
            case 'Start new file'
                [file,path]=uiputfile('*.txt','Save Results');
                fileID = fopen(fullfile(path,file),'w');
                fprintf(fileID,'%d %d %d %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\r\n',Results);
                fclose(fileID);
        end
end
guidata(hObject, handles)


% --- Executes on button press in Calculate.
function Calculate_Callback(hObject, eventdata, handles)
% hObject    handle to Calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%label the frame and time, *watch out for some earlier data files that stored frame numbers instead of the time

text = sprintf('Frame#%d (%.3fs)', handles.frameNum+handles.origStartingFrameNum - 1, handles.time(handles.frameNum));
set(handles.frameText, 'String', text)
errorLogic = false;
errorVals = zeros(1,handles.l);
for i = 1:handles.l
    curArea=bwarea(~handles.CellPrintBi(:,:,i))*(handles.PxScale^2);
    checkError = abs(handles.refArea(i) - curArea);
    errorVals(i) = checkError;
    if checkError > 5
        uiwait(errordlg('Area data and cell prints do not line up, please edit accordingly'))
        errorLogic = true;
        break
    end
end
if errorLogic
    close all
    return
end
%get the contour
frames=handles.l;
handles.sz=size(handles.CellPrintBi);
rgb_stack=zeros(handles.sz(1),handles.sz(2),3,handles.sz(3));
handles.OL=zeros(handles.sz(1),handles.sz(2),frames);
for f=1:frames
    rgb = cat(3, handles.CellPrintBi(:,:,f), handles.CellPrintBi(:,:,f), handles.CellPrintBi(:,:,f));%to be overlayed, converted to RGB
    handles.OL(:,:,f) = bwperim(~handles.CellPrintBi(:,:,f));%calculate the outline, ~is to flip the color for bwperim
    rgb_stack(:,:,:,f) = imoverlay(rgb,handles.OL(:,:,f),[1,0,0]);%overlay the outline
end

%plot the image with the contour
handles.CellPrintRGB = rgb_stack;
handles.currentImageOverlay = handles.CellPrintRGB(:,:,:,handles.frameNum);
axes(handles.origVideo)
axis square
imshow(handles.currentImageOverlay)


%calculate the roundness
CIn=zeros(frames,2);RIn=zeros(frames,1);COut=zeros(frames,2);ROut=zeros(frames,1);
handles.x=cell(frames,1);handles.y=cell(frames,1);
handles.xIn=cell(frames,1);handles.yIn=cell(frames,1);
handles.xOut=cell(frames,1);handles.yOut=cell(frames,1);
for f=1:frames
    [handles.y{f,1}, handles.x{f,1}] = find(handles.OL(:,:,f));
    [RIn(f,1),CIn(f,1),CIn(f,2)]=max_inscribed_circle(handles.x{f,1},handles.y{f,1},handles.OL(:,:,f));
    theta = linspace(0,2*pi,100);
    handles.xIn{f,1} = CIn(f,1) + RIn(f,1)*cos(theta);%find the inscribed circle
    handles.yIn{f,1} = CIn(f,2) + RIn(f,1)*sin(theta);
   
    %find the circumscribed circles 
    [COut(f,:),ROut(f,1)]=minboundcircle(handles.x{f,1},handles.y{f,1},false);
    theta = linspace(0,2*pi,100);
    handles.xOut{f,1} = COut(f,1) + ROut(f,1)*cos(theta);
    handles.yOut{f,1} = COut(f,2) + ROut(f,1)*sin(theta);
end
        
%Calculate the centroid of each contour
handles.centroids = zeros(frames,2);
for i = 1:frames
    s = regionprops(~handles.CellPrintBi(:,:,i));
    numAreas = length(s);
    if numAreas > 1
        testAreas = zeros(numAreas,1);
        for j = 1:numAreas
            testAreas(j) = s(j).Area;
        end
        [~, idx] = max(testAreas);
        curCentroid = s(idx).Centroid;
    else
        curCentroid = s.Centroid;
    end
    handles.centroids(i,:) = curCentroid;
end
% subtract stage movement
changeDisplacement = zeros(handles.l-1,0);
for i = 2:handles.l
    curX = handles.centroids(i,1);
    prevX = handles.centroids(i-1,1);
    curY = handles.centroids(i,2);
    prevY = handles.centroids(i-1,2);
    changeDisplacement(i-1) = sqrt((curX-prevX)^2 + (curY-prevY)^2) * handles.PxScale;
end
idx = find(abs(changeDisplacement)./diff(handles.time) > 0.5);
handles.correctedCentroids = handles.centroids;
for i = 1:length(idx)
    changeX = handles.correctedCentroids(idx(i)+1,1)-handles.correctedCentroids(idx(i),1);
    changeY = handles.correctedCentroids(idx(i)+1,2)-handles.correctedCentroids(idx(i),2);
    handles.correctedCentroids(idx(i)+1:end,:) = handles.correctedCentroids(idx(i)+1:end,:) - [changeX, changeY];
end
handles.zeroDisplacementFrame = 1;
initX = mean(handles.correctedCentroids(1:3,1));
initY = mean(handles.correctedCentroids(1:3,2));
handles.displacement = sqrt((handles.correctedCentroids(:,1) - initX).^2 + ...
    (handles.correctedCentroids(:,2) - initY).^2) * handles.PxScale;
uiwait(msgbox(sprintf('Number of stage movements detected: %d', length(idx))))
axes(handles.CapturedDemo)
plot(handles.x{handles.frameNum,1},handles.sz(1)-handles.y{handles.frameNum,1},'ro',...
    handles.xIn{handles.frameNum,1},handles.sz(1)-handles.yIn{handles.frameNum,1},'b-',...
    handles.xOut{handles.frameNum,1},handles.sz(1)-handles.yOut{handles.frameNum,1},'b-')
grid on
axis equal
pbaspect([1 1 1])

%add centroid
hold on
plot(handles.centroids(handles.frameNum,1), handles.sz(1)-handles.centroids(handles.frameNum,2), '*')
hold off
    
%plot the ellipticity, area, roundness, displacement
handles.AreaPlotted=[];
handles.displacementPlotted=[];
handles.RoundnessPlotted=[];
if isempty(handles.displacementPlotted)||isempty(handles.AreaPlotted)||isemtpy(handles.RoundnessPlotted)
    handles.Area=zeros(frames,1);
    handles.Roundness=zeros(frames,1);
    for f=1:frames
%         handles.Area(f,1)=bwarea(~handles.CellPrintBi(:,:,f))*(handles.PxScale^2);
        handles.Area(f,1) = handles.refArea(f);
%         handles.Peri(f,1)=periLength(~handles.CellPrintBi(:,:,f))*(handles.PxScale);%in um
%         handles.displacement(f,1)=4*pi*handles.Area(f,1)/(handles.Peri(f,1)^2);
        handles.Roundness(f,1)=RIn(f,1)/ROut(f,1);
    end
    %plot the displacement
    axes(handles.DisplacementOverTime)
    plot(handles.time,handles.displacement)
    grid on
    xlabel('Time [s]');
    ylabel(sprintf('Displacement'));
    handles.displacementPlotted=1;
    
    %plot the area
    axes(handles.AreaOverTime)
    plot(handles.time, handles.Area, 'LineWidth', 2, 'Color', [.6 .6 .6])
    grid on
    xlabel('Time [s]');
    ylabel(sprintf('Area [%s^2]',handles.PxScaleUnit));
    handles.AreaPlotted=1;
    
%     %plot the peri
%     axes(handles.PeriOverTime)
%     plot(handles.time, handles.Peri)
%     xlabel('Time [s]');
%     ylabel(sprintf('Peri [%s]',handles.PxScaleUnit));
%     handles.PeriPlotted=1;
    
   %plot the roundness
   axes(handles.RoundnessOverTime)
   plot(handles.time,handles.Roundness)
   grid on
   xlabel('Time [s]');
   ylabel(sprintf('Roundness'));
   handles.RoundnessPlotted=1;
end
%add displacement cursor
axes(handles.DisplacementOverTime)
hold on
handles.displacementCursor=plot(handles.time(handles.frameNum),handles.displacement(handles.frameNum),'x');
hold off 
%area cursor
axes(handles.AreaOverTime)
hold on
handles.areaCursor=plot(handles.time(handles.frameNum),handles.Area(handles.frameNum),'x');

%roundness cursor
axes(handles.RoundnessOverTime)
hold on
handles.roundnessCursor=plot(handles.time(handles.frameNum),handles.Roundness(handles.frameNum),'x');
hold off
%initialize the minAreaSlider and maxAreaSlider to later calculate roundness and ellipticity
[~, handles.minWindowFrame] = min(abs(handles.Area-50));
[~, handles.maxWindowFrame] = min(abs(handles.Area-100));
if handles.minWindowFrame == handles.maxWindowFrame
    handles.maxWindowFrame = handles.minWindowFrame + 5;
end
set(handles.minAreaSlider, 'Value', (handles.minWindowFrame-1)/handles.l);
set(handles.maxAreaSlider, 'Value', (handles.maxWindowFrame-1)/handles.l);
%display area and time
handles.minArea=handles.Area(handles.minWindowFrame);
handles.maxArea=handles.Area(handles.maxWindowFrame);
handles = updateReadout(handles);
% find max roundness and max displacement
windowSize = 5;
roundnessFilt = zeros(handles.l, 1);
displacementFilt = zeros(handles.l, 1);
halfWindow = (windowSize-1)/2;
windowIdx = -halfWindow:halfWindow;
for i = halfWindow+1:handles.l-halfWindow
    numIntRoundness = trapz(handles.time(windowIdx + i), handles.Roundness(windowIdx + i));
    numIntDisplacement = trapz(handles.time(windowIdx + i), handles.displacement(windowIdx + i));
    timeInt = handles.time(i+halfWindow) - handles.time(i-halfWindow);
    roundnessFilt(i) = numIntRoundness / timeInt;
    displacementFilt(i) = numIntDisplacement / timeInt;
end
[handles.maxRoundness, maxRoundIdx] = max(roundnessFilt);
[handles.maxDisplacement, maxDisplacementIdx] = max(displacementFilt);
handles.timeFrameMaxRoundness = handles.rawTime([maxRoundIdx-halfWindow,maxRoundIdx+halfWindow]);
handles.timeFrameMaxDisplacement = handles.rawTime([maxDisplacementIdx-halfWindow,maxDisplacementIdx+halfWindow]);
roundText = sprintf('Maximum Roundness = %.3f', handles.maxRoundness);
set(handles.MaxRoundnessText, 'String', roundText)
displacementText = sprintf('Maximum Displacement = %.3f', handles.maxDisplacement);
set(handles.MaxDisplacementText, 'String', displacementText)
guidata(hObject,handles)


% --- Executes on button press in linearFit.
function linearFit_Callback(hObject, eventdata, handles)
% hObject    handle to linearFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fitTime = handles.time(handles.minWindowFrame:handles.maxWindowFrame);
fitArea = handles.Area(handles.minWindowFrame:handles.maxWindowFrame);
p = polyfit(fitTime', fitArea, 1);
slope = p(1);
handles.spreadingSpeed = slope;
axes(handles.AreaOverTime)
cla
plot(handles.time, handles.Area, 'LineWidth', 2, 'Color', [.6 .6 .6])
grid on
xlabel('Time [s]');
ylabel(sprintf('Area [%s^2]',handles.PxScaleUnit));
handles.AreaPlotted=1;
hold on
plot(fitTime, polyval(p,fitTime))
handles = updateReadout(handles);
set(handles.speedText, 'String', sprintf('Spreading Speed = %.3f um^2/s', slope));
handles.spreadingSpeedData = [handles.rawTime(handles.minWindowFrame), handles.rawTime(handles.maxWindowFrame), 1, slope]; % '1' indicates linear fit
guidata(hObject, handles)

% --- Executes on button press in 90oidalFit.
function sigmoidalFit_Callback(hObject, eventdata, handles)
% hObject    handle to sigmoidalFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fitTime = handles.time(handles.minWindowFrame:handles.maxWindowFrame);
fitArea = handles.Area(handles.minWindowFrame:handles.maxWindowFrame);
sigmType = fittype('a + b ./ (1 + exp(-(t-t0)./dt))', 'independent', 't', 'dependent', 'area'); 
sigmFit = fit(fitTime', fitArea, sigmType, 'StartPoint', [min(fitArea), max(fitArea)-min(fitArea), 10, mean(fitTime)]);
slope = (sigmFit.b) / (4*sigmFit.dt);
handles.spreadingSpeed = slope;
axes(handles.AreaOverTime)
cla
plot(handles.time, handles.Area, 'LineWidth', 2, 'Color', [.6 .6 .6])
grid on
xlabel('Time [s]');
ylabel(sprintf('Area [%s^2]',handles.PxScaleUnit));
handles.AreaPlotted=1;
hold on
plot(fitTime, sigmFit(fitTime))
handles = updateReadout(handles);
set(handles.speedText, 'String', sprintf('Spreading Speed = %.3f um^2/s', slope));
handles.spreadingSpeedData = [handles.rawTime(handles.minWindowFrame), handles.rawTime(handles.maxWindowFrame), 2, slope]; % '2' indicates sigmoidal fit
guidata(hObject, handles)

% --- Executes on button press in maxArea.
function maxArea_Callback(hObject, eventdata, handles)
% hObject    handle to maxArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
windowTime = handles.time(handles.minWindowFrame:handles.maxWindowFrame);
windowArea = handles.Area(handles.minWindowFrame:handles.maxWindowFrame);
handles.maxAreaAvg = trapz(windowTime,windowArea) / (windowTime(end)-windowTime(1));
% handles.maxAreaAvg = mean(handles.Area(handles.minWindowFrame:handles.maxWindowFrame));
set(handles.maxAreaText, 'String', sprintf('Maximum Area = %.3f um^2', handles.maxAreaAvg));
handles.maxAreaData = [handles.rawTime(handles.minWindowFrame), handles.rawTime(handles.maxWindowFrame), handles.maxAreaAvg];
guidata(hObject, handles)


% --- Executes on button press in noSlope.
function noSlope_Callback(hObject, eventdata, handles)
% hObject    handle to noSlope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.spreadingSpeedData = [0, 0, 0, 0]; % 0's indicate no slope measured
set(handles.speedText, 'String', 'No spreading speed data');
guidata(hObject,handles)

% --- Executes on button press in noMaxArea.
function noMaxArea_Callback(hObject, eventdata, handles)
% hObject    handle to noMaxArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maxAreaData = [0, 0, 0]; % 0's indicate no max area measured
set(handles.maxAreaText, 'String', 'No max area data');
guidata(hObject,handles)


% --- Executes on button press in automaticWindow.
function automaticWindow_Callback(hObject, eventdata, handles)
% hObject    handle to automaticWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% first, apply sliding window filter to find max slope in area vs time
choice = questdlg('Only include frames after current selection for window?', 'Only later frames?', 'Yes', 'No, include all', 'Yes');
switch choice
    case 'Yes'
        startFrame = handles.frameNum+1;
        % recalculate max roundness
        windowSize = 5;
        roundnessFilt = zeros(handles.l, 1);
        halfWindow = (windowSize-1)/2;
        windowIdx = -halfWindow:halfWindow;
        for i = startFrame+halfWindow:handles.l-halfWindow
            numIntegral = trapz(handles.time(windowIdx + i), handles.Roundness(windowIdx + i));
            timeInt = handles.time(i+halfWindow) - handles.time(i-halfWindow);
            roundnessFilt(i) = numIntegral / timeInt;
        end
        [handles.maxRoundness, maxRoundIdx] = max(roundnessFilt);
        handles.timeFrameMax = handles.rawTime([maxRoundIdx-halfWindow,maxRoundIdx+halfWindow]);
        roundText = sprintf('Maximum Roundness = %.3f', handles.maxRoundness);
        set(handles.MaxRoundnessText, 'String', roundText)
    case 'No, include all'
        startFrame = 1;
end
AreaAveraged = handles.Area(startFrame:end);
for i = 2:length(AreaAveraged)-1
    AreaAveraged(i) = mean(handles.Area(i+startFrame-2:i+startFrame));
end
% use centered finite difference
AreaSlopeVals = zeros(length(AreaAveraged),1);
AreaSlopeVals(1) = (AreaAveraged(2)-AreaAveraged(1)) / (handles.time(startFrame+1)-handles.time(startFrame));
AreaSlopeVals(end) = (AreaAveraged(end)-AreaAveraged(end-1)) / (handles.time(end) - handles.time(end-1));
for i = 2:length(AreaAveraged)-1
    AreaSlopeVals(i) = (AreaAveraged(i+1)-AreaAveraged(i-1)) / (handles.time(startFrame+i)-handles.time(startFrame+i-2));
end
[~,maxIdx] = max(AreaSlopeVals);
handles.maxSlopeIdx = maxIdx + startFrame - 1;
% find range from ~15 s before max slope to ~15 s after
handles.maxSlopeTime = handles.time(handles.maxSlopeIdx);
[~,handles.minWindowFrame] = min(abs(handles.time - (handles.maxSlopeTime-15)));
[~,handles.maxWindowFrame] = min(abs(handles.time - (handles.maxSlopeTime+15)));
set(handles.minAreaSlider, 'Value', (handles.minWindowFrame-1)/handles.l)
set(handles.maxAreaSlider, 'Value', (handles.maxWindowFrame-1)/handles.l)
handles = updateReadout(handles);
guidata(hObject,handles)

% --- Executes on button press in correctStageMovement.
function correctStageMovement_Callback(hObject, eventdata, handles)
% hObject    handle to correctStageMovement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Is the displacement jump before or after this point?', 'When is the jump?', 'Before', 'After', 'Cancel', 'After');
switch choice
    case 'Before'
        indices = handles.frameNum-1:handles.frameNum;
    case 'After'
        indices = handles.frameNum:handles.frameNum+1;
    otherwise
        indices = NaN;
end
if ~isnan(indices)
    changeX = handles.correctedCentroids(indices(2),1)-handles.correctedCentroids(indices(1),1);
    changeY = handles.correctedCentroids(indices(2),2)-handles.correctedCentroids(indices(1),2);
    handles.correctedCentroids(indices(2):end,:) = handles.correctedCentroids(indices(2):end,:) - [changeX, changeY];
    if handles.zeroDisplacementFrame > 1
        initX = mean(handles.correctedCentroids(handles.zeroDisplacementFrame-1:handles.zeroDisplacementFrame+1,1));
        initY = mean(handles.correctedCentroids(handles.zeroDisplacementFrame-1:handles.zeroDisplacementFrame+1,2));
    else
        initX = mean(handles.correctedCentroids(1:3,1));
        initY = mean(handles.correctedCentroids(1:3,2));
    end
    handles.displacement = sqrt((handles.correctedCentroids(:,1) - initX).^2 + ...
        (handles.correctedCentroids(:,2) - initY).^2) * handles.PxScale;
    axes(handles.DisplacementOverTime)
    plot(handles.time,handles.displacement)
    grid on
    xlabel('Time [s]');
    ylabel(sprintf('Displacement'));
    handles.displacementPlotted=1;
    %add displacement cursor
    axes(handles.DisplacementOverTime)
    hold on
    handles.displacementCursor=plot(handles.time(handles.frameNum),handles.displacement(handles.frameNum),'x');
    hold off
    handles = updateReadout(handles);
    % find new max displacement
    windowSize = 5;
    displacementFilt = zeros(handles.l, 1);
    halfWindow = (windowSize-1)/2;
    windowIdx = -halfWindow:halfWindow;
    for i = halfWindow+1:handles.l-halfWindow
        numIntDisplacement = trapz(handles.time(windowIdx + i), handles.displacement(windowIdx + i));
        timeInt = handles.time(i+halfWindow) - handles.time(i-halfWindow);
        displacementFilt(i) = numIntDisplacement / timeInt;
    end
    [handles.maxDisplacement, maxDisplacementIdx] = max(displacementFilt);
    handles.timeFrameMaxDisplacement = handles.rawTime([maxDisplacementIdx-halfWindow,maxDisplacementIdx+halfWindow]);
    displacementText = sprintf('Maximum Displacement = %.3f', handles.maxDisplacement);
    set(handles.MaxDisplacementText, 'String', displacementText)
end
guidata(hObject,handles)

% --- Executes on button press in measureMaxes.
function measureMaxes_Callback(hObject, eventdata, handles)
% hObject    handle to measureMaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Use current window or select automated?', 'Window selection', 'Use current', 'Select automated', 'Select automated');
switch choice
    case 'Select automated'
        choice = questdlg('Only include frames after current selection for window?', 'Only later frames?', 'Yes', 'No, include all', 'Yes');
        switch choice
            case 'Yes'
                startFrame = handles.frameNum+1;
            case 'No, include all'
                startFrame = 1;
        end
        AreaAveraged = handles.Area(startFrame:end);
        for i = 2:length(AreaAveraged)-1
            AreaAveraged(i) = mean(handles.Area(i+startFrame-2:i+startFrame));
        end
        % use centered finite difference
        AreaSlopeVals = zeros(length(AreaAveraged),1);
        AreaSlopeVals(1) = (AreaAveraged(2)-AreaAveraged(1)) / (handles.time(startFrame+1)-handles.time(startFrame));
        AreaSlopeVals(end) = (AreaAveraged(end)-AreaAveraged(end-1)) / (handles.time(end) - handles.time(end-1));
        for i = 2:length(AreaAveraged)-1
            AreaSlopeVals(i) = (AreaAveraged(i+1)-AreaAveraged(i-1)) / (handles.time(startFrame+i)-handles.time(startFrame+i-2));
        end
        higherSlopeLogic = AreaSlopeVals > 0.4;
        for i = 1:length(higherSlopeLogic)-2
            if higherSlopeLogic(i) && all(higherSlopeLogic(i+1:i+2))
                startIdx = i;
                break
            end
        end
        if ~exist('startIdx')
            uiwait(msgbox('Could not find start of spreading'))
            startIdx = 1;
        end
        for i = startIdx:length(higherSlopeLogic)-2
            if ~higherSlopeLogic(i) && ~any(higherSlopeLogic(i+1:i+2))
                endIdx = i;
                break
            end
        end
        if ~exist('endIdx')
            uiwait(msgbox('Could not find end of spreading'))
            endIdx = length(higherSlopeLogic);
        end
        handles.startSpreadFrame = startIdx + startFrame - 1;
        handles.endSpreadFrame = endIdx + startFrame - 1;
        % set window and update readout
        handles.minWindowFrame = handles.startSpreadFrame;
        handles.maxWindowFrame = handles.endSpreadFrame;
        set(handles.minAreaSlider, 'Value', (handles.minWindowFrame-1)/handles.l)
        set(handles.maxAreaSlider, 'Value', (handles.maxWindowFrame-1)/handles.l)
        handles = updateReadout(handles);
    case 'Use current'
        handles.startSpreadFrame = handles.minWindowFrame;
        handles.endSpreadFrame = handles.maxWindowFrame;
end
handles.spreadingTimeFrame = handles.rawTime([handles.startSpreadFrame, handles.endSpreadFrame]);

% set zero displacement at start of spreading
handles.zeroDisplacementFrame = handles.startSpreadFrame;
if handles.startSpreadFrame == 1
    initX = mean(handles.correctedCentroids(1:3,1));
    initY = mean(handles.correctedCentroids(1:3,2));
else
    initX = mean(handles.correctedCentroids(handles.startSpreadFrame-1:handles.startSpreadFrame+1,1));
    initY = mean(handles.correctedCentroids(handles.startSpreadFrame-1:handles.startSpreadFrame+1,2));
end
handles.displacement = sqrt((handles.correctedCentroids(:,1) - initX).^2 + ...
    (handles.correctedCentroids(:,2) - initY).^2) * handles.PxScale;
axes(handles.DisplacementOverTime)
plot(handles.time,handles.displacement)
grid on
xlabel('Time [s]');
ylabel(sprintf('Displacement'));
handles.displacementPlotted=1;
%add displacement cursor
axes(handles.DisplacementOverTime)
hold on
handles.displacementCursor=plot(handles.time(handles.frameNum),handles.displacement(handles.frameNum),'x');
hold off

% find max roundness and max displacement on window
windowSize = 5;
displacementFilt = zeros(handles.l, 1);
roundnessFilt = zeros(handles.l, 1);
halfWindow = (windowSize-1)/2;
windowIdx = -halfWindow:halfWindow;
if handles.minWindowFrame > halfWindow
    startLoop = handles.minWindowFrame;
else
    startLoop = halfWindow + 1;
end
if handles.maxWindowFrame <= handles.l-halfWindow
    endLoop = handles.maxWindowFrame;
else
    endLoop = handles.l-halfWindow;
end
for i = startLoop:endLoop
    numIntRoundness = trapz(handles.time(windowIdx + i), handles.Roundness(windowIdx + i));
    numIntDisplacement = trapz(handles.time(windowIdx + i), handles.displacement(windowIdx + i));
    timeInt = handles.time(i+halfWindow) - handles.time(i-halfWindow);
    roundnessFilt(i) = numIntRoundness / timeInt;
    displacementFilt(i) = numIntDisplacement / timeInt;
end
[handles.maxRoundness, maxRoundIdx] = max(roundnessFilt);
[handles.maxDisplacement, maxDisplacementIdx] = max(displacementFilt);
handles.timeFrameMaxRoundness = handles.rawTime([maxRoundIdx-halfWindow,maxRoundIdx+halfWindow]);
handles.timeFrameMaxDisplacement = handles.rawTime([maxDisplacementIdx-halfWindow,maxDisplacementIdx+halfWindow]);
roundText = sprintf('Maximum Roundness = %.3f', handles.maxRoundness);
set(handles.MaxRoundnessText, 'String', roundText)
displacementText = sprintf('Maximum Displacement = %.3f', handles.maxDisplacement);
set(handles.MaxDisplacementText, 'String', displacementText)
% find overall max displacement
overallDisplFilt = zeros(handles.l,1);
for i = startLoop:handles.l-halfWindow
    numIntOverallDispl = trapz(handles.time(windowIdx + i), handles.displacement(windowIdx + i));
    timeInt = handles.time(i+halfWindow) - handles.time(i-halfWindow);
    overallDisplFilt(i) = numIntOverallDispl / timeInt;
end
[handles.maxOverallDispl, maxOverallDisplIdx] = max(overallDisplFilt);
handles.timeFrameMaxOverallDispl = handles.rawTime([maxOverallDisplIdx-halfWindow,maxOverallDisplIdx+halfWindow]);
handles = updateReadout(handles);
guidata(hObject,handles)

function handles = updateReadout(handles)
handles.minArea = handles.Area(handles.minWindowFrame);
handles.maxArea = handles.Area(handles.maxWindowFrame);
minTime = handles.time(handles.minWindowFrame);
maxTime = handles.time(handles.maxWindowFrame);
text = sprintf('Lower area limit is %.3fum^2 \nat time %.3fs\n\nUpper area limit is %.3fum^2 \nat time %.3fs',handles.minArea,minTime,handles.maxArea,maxTime);
set(handles.AreaText, 'String', text)

%get the roundness array sitting in between the min and max frames and find
%the average
windowTime = handles.time(handles.minWindowFrame:handles.maxWindowFrame);
windowRoundness = handles.Roundness(handles.minWindowFrame:handles.maxWindowFrame);
handles.roundnessAvg = trapz(windowTime', windowRoundness) / (windowTime(end)-windowTime(1));
%handles.roundnessAvg=mean(handles.Roundness(handles.minWindowFrame:handles.maxWindowFrame));
text = sprintf('Average Roundness = %.3f', handles.roundnessAvg);
set(handles.AverageRoundness, 'String', text)

%change the vertical lines for the AreaOverTime
axes(handles.AreaOverTime)
areaLim = get(gca, 'YLim');
hold on
if isfield(handles, 'lineAreaMin')
    delete(handles.lineAreaMin);
    delete(handles.lineAreaMax);
end
handles.lineAreaMin=line([minTime, minTime], areaLim, 'Color', 'green', 'LineStyle', '--');
handles.lineAreaMax=line([maxTime, maxTime], areaLim ,'Color', 'green', 'LineStyle', '--');
hold off

%plot the vertical lines in RoundnessOverTime
axes(handles.RoundnessOverTime)
roundnessLim = get(gca, 'YLim');
hold on
if isfield(handles, 'lineRoundnessMin')
    delete(handles.lineRoundnessMin);
    delete(handles.lineRoundnessMax);
end
handles.lineRoundnessMin=line([minTime, minTime], roundnessLim, 'Color', 'green', 'LineStyle', '--');
handles.lineRoundnessMax=line([maxTime, maxTime], roundnessLim, 'Color', 'green', 'LineStyle', '--');
hold off

%plot the vertical lines in DisplacementOverTime
axes(handles.DisplacementOverTime)
displacementLim = get(gca, 'YLim');
hold on
if isfield(handles, 'lineDisplacementMin')
    delete(handles.lineDisplacementMin);
    delete(handles.lineDisplacementMax);
end
handles.lineDisplacementMin=line([minTime, minTime], displacementLim, 'Color', 'green', 'LineStyle', '--');
handles.lineDisplacementMax=line([maxTime, maxTime], displacementLim, 'Color', 'green', 'LineStyle', '--');
hold off
