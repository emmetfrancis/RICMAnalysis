%% view 2 cell contour outlines
imPath = uigetdir('','Choose image seq directory');
imSeq = dir(fullfile(imPath,'*.bmp'));
imFiles = {imSeq.name};
[printFile,printPath] = uigetfile('*.mat','Choose original print file');
oldPrints = load(fullfile(printPath,printFile));
oldPrints = oldPrints.cellPrints;
newPrintPath = uigetdir('','Choose new print directory');
newPrintSeq = dir(fullfile(newPrintPath,'*.tif'));
newPrintFiles = {newPrintSeq.name};

% load analysis file to find which frames were analyzed
[areaFile,areaPath] = uigetfile('*.txt','Choose area vs. frame number file');
areaData = load(fullfile(areaPath,areaFile));
frameNums = areaData(:,1);
oldAreaCompare = areaData(:,2);
oldAreas = zeros(size(frameNums));
newAreas = zeros(size(frameNums));

%%
% print each separately and save as colored bmps
savePath = uigetdir('','Choose folder to save new images');
figure
colormap gray
cmap = colormap;
oldColor = 'r';
newColor = 'g';
for i = 1:length(frameNums)
    imIdx  = frameNums(i);
    curIm = imread(fullfile(imPath,imFiles{imIdx}));
    imagesc(curIm,[90 190])
    daspect([1 1 1])
    hold on
    oldPrintCur = ~oldPrints(:,:,i); %inverse image so cell body is true
    newPrintCur = imread(fullfile(newPrintPath,newPrintFiles{i}));
    newPrintCur = logical(newPrintCur);
    oldAreas(i) = sum(oldPrintCur(:));
    oldTrace = tracesFromBinary(oldPrintCur,100,0);
    patch('XData',oldTrace(:,2),'YData',oldTrace(:,1),'EdgeColor',oldColor,'FaceColor','none','LineWidth',1)
    oldBoundary = boundaryFromVertices(oldTrace(:,1),oldTrace(:,2));
    [newTrace,newPrintRev] = tracesFromBinary(newPrintCur,100,0);
    newAreas(i) = sum(newPrintRev(:));
    patch('XData',newTrace(:,2),'YData',newTrace(:,1),'EdgeColor',newColor,'FaceColor','none','LineWidth',1)
    drawnow
    newBoundary = boundaryFromVertices(newTrace(:,1),newTrace(:,2));
    % create new image to save (same size)
%     curIm = imadjust(curIm,[90/255,190/255]);
%     rgbIm = ind2rgb(curIm,cmap);
%     for j = 1:length(oldBoundary)
%         rgbIm(oldBoundary(j,1),oldBoundary(j,2),1) = 1; % set red channel to 1 at old border
%     end
%     for j = 1:length(newBoundary)
%         rgbIm(newBoundary(j,1),newBoundary(j,2),2) = 1; % set green channel to 1 at new border
%     end
    rgbIm = print('-RGBImage');
    newFileName = sprintf('OverlayOL%d_%s',i,imFiles{imIdx});
%     imwrite(rgbIm,fullfile(savePath,newFileName))
    hold off
end