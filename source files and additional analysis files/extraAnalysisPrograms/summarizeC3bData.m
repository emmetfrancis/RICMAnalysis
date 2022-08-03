dataFolder = uigetdir('','Choose C3b data parent folder');
allFolders = dir(dataFolder);
allFolderNames = {allFolders.name};
% [timeFile, timePath] = uigetfile('*.dat','Choose start and end time data');
% startEnd = readcell(fullfile(timePath,timeFile));
% startEnd = cell2mat(startEnd(2:end,:));
folderNum = 1;
C3bAreaData = {};
figure
for i = 1:length(allFolderNames)
    curFileStruct = dir(strcat(dataFolder,'\',allFolderNames{i},'\*.txt'));
    if strcmp(allFolderNames{i},'.') || strcmp(allFolderNames{i},'..') || isempty(curFileStruct)
        continue
    end
%     if isempty(curFileStruct) || strcmp(allFolderNames{i},'10') || strcmp(allFolderNames{i},'17') || strcmp(allFolderNames{i},'18')
%         continue
%     end
    fileName = {curFileStruct.name};
    curData = load(fullfile(dataFolder,allFolderNames{i},fileName{1}));
%     if curData(1,2) > 50
%         timeCorrect = 20;
%         curData(:,1) = curData(:,1) - curData(1,1) + (curData(1,2)-20)/1.5;
% %         startIdx = find(curData(:,2)>50,1,'first');
% %         curData(:,1) = curData(:,1) - curData(startIdx,1);
%     else
        startIdx = find(curData(:,2)>100,1,'first');
        curData(:,1) = curData(:,1) - curData(startIdx,1);
%         curData(:,1) = curData(:,1) - startEnd(folderNum,1);
%     end
    C3bAreaData{folderNum} = curData;
    folderNum = folderNum+1;
    plot(curData(:,1),curData(:,2))
    hold on
end

%% find avg curve
timeSample = -150:12:250;
compileAreaCell = cell(1,length(timeSample));
figure
for i = 1:length(C3bAreaData)
    curTime = C3bAreaData{i}(:,1);
    curArea = C3bAreaData{i}(:,2);
    plot(curTime,curArea)
    hold on
    if curTime(1) > timeSample(1)
        startIdx = find(timeSample > curTime(1),1,'first');
    else
        startIdx = 1;
    end
    if curTime(end) < timeSample(end)
        endIdx = find(timeSample < curTime(end),1,'last');
    else
        endIdx = length(timeSample);
    end
    interpArea = interp1q(curTime,curArea,timeSample(startIdx:endIdx)');
    for j = startIdx:endIdx
        compileAreaCell{j} = [compileAreaCell{j},interpArea(j-startIdx+1)];
    end
end
avgAreaCell = zeros(1,length(timeSample));
stdAreaCell = zeros(1,length(timeSample));
for i = 1:length(avgAreaCell)
    avgAreaCell(i) = mean(compileAreaCell{i});
%     stdAreaCell(i) = std(compileAreaCell{i});
    stdAreaCell(i) = std(compileAreaCell{i})/sqrt(length(compileAreaCell{i}));
end
errorbar(timeSample+122,avgAreaCell,stdAreaCell,'LineWidth',2) 