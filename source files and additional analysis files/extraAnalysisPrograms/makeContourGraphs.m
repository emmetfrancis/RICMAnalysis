[file,path] = uigetfile('*.mat','Choose cell prints .mat file');
loadStruct = open(fullfile(path,file));
try cellPrints = loadStruct.cellPrints;
catch
    cellPrints = loadStruct.CellPrints;
end
[areaFile,areaPath] = uigetfile('*.txt',sprintf('Choose area file for %s',file));
areaMat = load(fullfile(areaPath,areaFile));
if ~any(mod(areaMat(:,1),1))
    frameRate = inputdlg('How many s per frame?');
    startFrame = inputdlg('Starting frame in figure (t = 0)?');
    time = (areaMat(:,1) - str2double(startFrame{1})) * str2double(frameRate{1});
%     time = time - time(1);
else
    time = areaMat(:,1) - areaMat(1,1);
end
% xCur = displMat(:,4);
% yCur = displMat(:,5);
[scaleFile,scalePath] = uigetfile('*.mat',sprintf('Load scale for %s',file));
getScale  = load(fullfile(scalePath,scaleFile));
scaleVal = getScale.scale;
% select every ~30 s
timeVals = [0,40,92,134,180];
idx = zeros(1,length(timeVals));
for i = 1:length(timeVals)
    [~,idx(i)] = min(abs(time-timeVals(i)));
end
% idx = 1;
% for i = 2:length(time)
%     if time(i) >= time(idx(end)) + 30
%         idx = [idx,i];
%     end
% end
cmap = parula(256);
figure
% subplot(1,2,1)
for i = 1:length(idx)
    curPrint = ~cellPrints(:,:,idx(i))';
    curPrint = flip(curPrint,2); 
    curColorIdx = round(time(idx(i))*255/time(idx(end))) + 1;
    curColor = cmap(curColorIdx,:);
    curTrace = scaleVal .* tracesFromBinary(curPrint,30,0);
    patch('XData',curTrace(:,1),'YData',curTrace(:,2),'FaceColor','none','EdgeColor',curColor,'LineWidth',3)
    hold on
    daspect([1 1 1])
end
colorbar('Ticks',[0 1],'TickLabels',{'t = 0 s', sprintf('t = %.0f s',time(idx(end)))})
xL = get(gca,'XLim');
yL = get(gca,'YLim');

frames = size(cellPrints,3);
centroids = zeros(frames,2);
for i = 1:frames
    s = regionprops(~cellPrints(:,:,i));
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
    centroids(i,:) = curCentroid;
end
xCur = centroids(:,1) * scaleVal;
yCur = centroids(:,2) * scaleVal;
yCur = size(cellPrints,1)*scaleVal-yCur;
% subplot(1,2,2)
for i = idx(1):idx(end)
    curColorIdx = round(time(i)*255/time(idx(end))) + 1;
    plot(xCur(i:i+1),yCur(i:i+1),'Color',cmap(curColorIdx,:),'LineWidth',3)
    daspect([1 1 1])
    hold on
end
% xLNew = [mean(xCur)-diff(xL)/8,mean(xCur)+diff(xL)/8];
% yLNew = [mean(yCur)-diff(yL)/8,mean(yCur)+diff(yL)/8];
% xlim(xLNew)
% ylim(yLNew)
xlim([0 scaleVal*size(cellPrints,2)])
ylim([0 scaleVal*size(cellPrints,1)])
set(gca,'FontSize',16)