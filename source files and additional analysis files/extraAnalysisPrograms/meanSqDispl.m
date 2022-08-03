% mean square displacements of cell centroids
%% Data loading
% load displacement data into displStruct
% listFiles = dir('*.txt');
% fileNames = {listFiles.name};
% for i = 1:length(fileNames)
%     file = fileNames{i};
%     curSite = str2double(file(1));
%     num = length(displStruct(curSite).displMat)+1;
%     displStruct(curSite).displMat{num} = load(file);
%     displStruct(curSite).startTimes(num) = startTimes(i);
%     displStruct(curSite).endTimes(num) = endTimes(i);
% end

displStruct = struct;
for i = 1:4
    path = uigetdir(cd,sprintf('Choose folder with site %d data',i));
    allDays = dir(path);
    allDaysFolders = {allDays.name};
    curIdx = 1;
    [timeFile,timePath] = uigetfile('*.txt','Choose file with start and end times');
    cellIn = readcell(fullfile(timePath,timeFile));
    for j = 1:length(allDaysFolders)
        if strcmp(allDaysFolders{j}(1),'.')
            continue
        end
        curPath = strcat(path,'\',allDaysFolders{j});
        curTxtFiles  = dir(strcat(curPath,'\*.txt'));
        curTxtFiles = {curTxtFiles.name};
        for k = 1:length(curTxtFiles)
            if strcmp(curTxtFiles{k}(1),'.')
                continue
            end
            if any(isstrprop(curTxtFiles{k}(1:end-4),'alpha')) % Roundness files have RoundDispl tagged on
                curData = load(fullfile(curPath,curTxtFiles{k}));
                displStruct(i).folder{curIdx} = curPath;
                displStruct(i).file{curIdx} = curTxtFiles{k};
                if strcmp(cellIn{curIdx,1},'--')
                    startTime = nan;
                else
                    startTime = cellIn{curIdx,1};
                end
                if strcmp(cellIn{curIdx,2},'--')
                    endTime = nan;
                else
                    endTime = cellIn{curIdx,2};
                end
                displStruct(i).startTimes(curIdx) = startTime;
                displStruct(i).endTimes(curIdx) = endTime;
                displStruct(i).displMat{curIdx} = curData;
                curIdx = curIdx + 1;
            end
        end
    end
end

%% Plot curves
figure
for k = 1:4
    curSite = k;
    subplot(2,2,k)
    for i = 1:length(displStruct(curSite).startTimes)
        curMat = displStruct(curSite).displMat{i};
        tStart = displStruct(curSite).startTimes(i);
        tEnd = displStruct(curSite).endTimes(i);
        if isnan(tStart) || isnan(tEnd)
            continue
        end
        time = curMat(:,1);
        displ = curMat(:,3);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        plot(time(startIdx:endIdx)-time(startIdx),displ(startIdx:endIdx),'-o')
        hold on
    end
    xlim([0 100])
    ylim([0 5])
end

%% Calculate average curves
allSpeeds = cell(1,4);
for k = 1:4
    curSite = k;
    alignedDispl = cell(1,51);
    allSpeeds{k} = [];
    for i = 1:length(displStruct(curSite).startTimes)
        curMat = displStruct(curSite).displMat{i};
        tStart = displStruct(curSite).startTimes(i);
        tEnd = displStruct(curSite).endTimes(i);
        if isnan(tStart) || isnan(tEnd)
            continue
        end
        time = curMat(:,1);
        displ = curMat(:,3);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        relTime = time(startIdx:endIdx)-time(startIdx);
        displ = displ(startIdx:endIdx);
        maxTime = min([100,relTime(end)]);
        interpTimes = 0:2:maxTime;
        interpDispl = interp1q(relTime,displ,interpTimes');
        for j = 1:maxTime/2+1
            L = length(alignedDispl{j});
            alignedDispl{j}(L+1) = interpDispl(j)^2;
        end
        allSpeeds{k} = horzcat(allSpeeds{k},diff(interpDispl')/2);
    end
    displStruct(curSite).meanSqDispl = zeros(1,51);
    displStruct(curSite).stdSqDispl = zeros(1,51);
    for i = 1:51
        displStruct(curSite).meanSqDispl(i) = mean(alignedDispl{i});
        displStruct(curSite).stdSqDispl(i) = std(alignedDispl{i})/sqrt(length(alignedDispl{i}));
    end
    displStruct(curSite).alignedDispl = alignedDispl;
end
figure
time = 0:2:100;
titleCell = {'100%','10%','1%','0.1%'};
colorCell = {[0,0.5,0],[0.725,0.1412,0.4784],[0,0,1],[1,0,0]};
for i = 1:4
%     subplot(2,2,(4-i+1))
    curLine = errorbar(time(1:2:end),displStruct(i).meanSqDispl(1:2:end)...
        ,displStruct(i).stdSqDispl(1:2:end));
    hold on
%     title(sprintf('%s rabbit IgG',titleCell{i}))
    ylim([0 15])
    xlabel('rel time (s)')
    ylabel('mean square displacement (\mum^2)')
    set(curLine,'LineWidth',2)
    set(curLine,'Color',colorCell{i})
    curFit = polyfit(log(time(2:end)),log(displStruct(i).meanSqDispl(2:end)),1);
    slope(i) = curFit(1);
    plot(time(2:end),exp(polyval(curFit,log(time(2:end)))),'Color',colorCell{i},'LineWidth',1)
end
legend('100%','100% Fit','10%','10% Fit',...
    '1%','1% Fit','0.1%','0.1% Fit')
set(gca,'XScale','log','YScale','log')
xlim([10 100])
ylim([0.1 20])

%% Model fits
figure
slope = zeros(1,4);
persistFit = cell(1,4);
time = 0:2:100;
for k = 1:4
%     time = [];
%     totDispl = [];
%     for i = 1:51
%         curData = displStruct(k).alignedDispl{i};
%         time = horzcat(time,2*(i-1)*ones(1,length(curData)));
%         totDispl = horzcat(totDispl,curData);
%     end
    curLine = plot(time,displStruct(k).meanSqDispl,'*');
    hold on
    logFit = polyfit(log(time(2:25)),log(displStruct(k).meanSqDispl(2:25)),1);
    slope(k) = logFit(1);
    displVals = exp(polyval(logFit,log(time(2:end))));
    curColor = get(curLine,'Color');
    plot(time(2:end),displVals,'Color',curColor)
    powerType = fittype('a*t^b','independent','t','dependent','msd');
    fitTime = [];
    fitDispl = [];
    for i = 1:26
        fitTime = [fitTime,(i-1)*2*ones(1,length(displStruct(k).alignedDispl{i}))];
        fitDispl = [fitDispl,displStruct(k).alignedDispl{i}];
    end
    powerFit{k} = fit(fitTime',fitDispl',powerType,'StartPoint',[.005,1]);
%     persistType = fittype('2*S.^2.*P .* (t - P.*(1-exp(-t./P)))', 'independent', 't', 'dependent', 'msd'); 
%     persistFit{k} = fit(time', displStruct(k).meanSqDispl', persistType, 'StartPoint', [0.1, 10]);
%     plot(time,persistFit{k}(time),'Color',curColor)
end
legend('100%','100% Fit','10%','10% Fit','1%','1% Fit','0.1%','0.1% Fit')

%% Calculate direction autocorrelation
figure
for k = 1:4
    displStruct(k).lengthAC = cell(1,length(displStruct(k).displMat));
    displStruct(k).AC = cell(1,length(displStruct(k).displMat));
    displStruct(k).ACLogFit = cell(1,length(displStruct(k).displMat));
%     figure
    for i = 1:length(displStruct(k).displMat)
        curMat = displStruct(k).displMat{i};
        tStart = displStruct(k).startTimes(i);
        tEnd = displStruct(k).endTimes(i);
        time = curMat(:,1);
        displ = curMat(:,3);
        xCur = curMat(:,4);
        yCur = curMat(:,5);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        relTime = time(startIdx:endIdx)-time(startIdx);
        displ = displ(startIdx:endIdx);
        xCur = xCur(startIdx:endIdx);
        yCur = yCur(startIdx:endIdx);
        arcLength = zeros(1,length(xCur));
        for j = 2:length(xCur)
            arcLength(j) = arcLength(j-1) + sqrt((xCur(j)-xCur(j-1))^2 + (yCur(j)-yCur(j-1))^2);
        end
        maxLength = arcLength(end);
        %         lengthInterval = maxLength/1000;
        lengthInterval = .001;
        interpLength = (0:lengthInterval:maxLength)';
        interpX = interp1q(arcLength',xCur,interpLength);
        interpY = interp1q(arcLength',yCur,interpLength);
        angle = atan2(diff(interpY),diff(interpX));
        lengthAC = 0:lengthInterval:maxLength;
        AC = ones(1,length(lengthAC));
        AC(end-2:end-1) = 0;
        for j = 2:length(AC)-2
            shiftIdx = find(abs(interpLength-lengthAC(j))<maxLength/10000);
            shiftedAng = angle(shiftIdx:end);
            curAng = angle(1:length(shiftedAng));
            curLength = interpLength(1:length(shiftedAng));
            AC(j) = trapz(curLength,cos(shiftedAng-curAng))/(curLength(end));
        end
%         goodFit = false;
%         while ~goodFit
%             maxPlotIdx = find(lengthAC < maxLength/10,1,'last');
%             plotRange = 1:maxPlotIdx;
%             semilogy(lengthAC(plotRange),AC(plotRange))
%             hold on
%             [lengthStart,~] = selectPoint(lengthAC(plotRange),AC(plotRange));
%             [lengthEnd,~] = selectPoint(lengthAC(plotRange),AC(plotRange));
%             minFitIdx = find(abs(lengthStart-lengthAC)<maxLength/10000);
%             maxFitIdx = find(abs(lengthEnd-lengthAC)<maxLength/10000);
%             fitRange = minFitIdx:maxFitIdx;
%             curFit = polyfit(lengthAC(fitRange),log(AC(fitRange)),1);
%             semilogy(lengthAC(fitRange),exp(polyval(curFit,lengthAC(fitRange))))
%             hold off
%             choice = questdlg('Good fit or try again?','Good fit?','Good fit','Try again','Good fit');
%             if strcmp(choice,'Good fit')
%                 goodFit = true;
%             end
%         end
        displStruct(k).lengthAC{i} = lengthAC;
        displStruct(k).AC{i} = AC;
%         displStruct(k).ACLogFit{i} = curFit;
        % find persistence length parameter
        smoothAC = smooth(lengthAC,AC,11,'loess');
        decayLengthIdx = find(smoothAC < 0.8, 1, 'first');
        decayLength = lengthAC(decayLengthIdx);
        plot(lengthAC,AC)
        hold on
        plot(lengthAC,smoothAC,'LineStyle',':','LineWidth',1.5)
        plot([decayLength,decayLength],[0,1])
        xlim([0 max(lengthAC)/5])
        hold off
        displStruct(k).decayLength(i) = decayLength;
%         uiwait(msgbox(sprintf('Decay length for site %d cell %d',k,i)))
    end
end

%% mean direction AC
lengthVals = 0:.001:.12;
numAC = length(lengthVals);
figure
for k = 1:4
    ACMat = zeros(length(displStruct(k).AC),numAC);
    for i = 1:length(displStruct(k).AC)
        if length(displStruct(k).AC{i}) < numAC
            L = length(displStruct(k).AC{i});
        else
            L = numAC;
        end
        ACMat(i,1:L) = displStruct(k).AC{i}(1:L);
    end
    displStruct(k).meanLengthAC = mean(ACMat,1);
    displStruct(k).stdLengthAC = std(ACMat,1)./sqrt(size(ACMat,1));
%     plot(lengthVals,displStruct(k).meanLengthAC)
    curLine = errorbar(lengthVals(1:5:end),displStruct(k).meanLengthAC(1:5:end)...
        ,displStruct(k).stdLengthAC(1:5:end));
    title('Length autocorrelation')
    ylim([0 15])
    xlabel('length (\mum)')
    ylabel('autocorrelation')
    set(curLine,'LineWidth',2)
    set(curLine,'Color',colorCell{k})
    hold on
%     curFit = polyfit(timeVals,log(displStruct(k).meanTimeAC),1);
%     plot(timeVals,exp(polyval(curFit,timeVals)))
end
hold off
ylim([0 1])
legend('100%','10%','1%','0.1%')


%% Time AC
allSpeeds = cell(1,4);
allSpeedsNew = cell(1,4);
allSpeedsNewSeparate = cell(1,4);
pTime = cell(1,4);
AUCVals = cell(1,4);
totDist = cell(1,4);
surfStrings = {'100','10','1','0.1'};
figure
for k = 1:4
    displStruct(k).allTimeAC = cell(1,50);
    displStruct(k).allTimeAC{1} = 1;
%     subplot(2,2,k)
    title(sprintf('Direction autocorrelation on %s%% IgG',surfStrings{k}))
    xlabel('\Deltat (s)')
    ylabel('Direction autocorrelation')
    hold on
    totIt = length(displStruct(k).displMat);
    if k == 3
        totIt = 10;
    end
    for i = 1:totIt
        curMat = displStruct(k).displMat{i};
        tStart = displStruct(k).startTimes(i);
        tEnd = displStruct(k).endTimes(i);
        if isnan(tStart) || isnan(tEnd)
            continue
        end
        time = curMat(:,1);
        displ = curMat(:,3);
        xCur = curMat(:,4);
        yCur = curMat(:,5);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        altStart = find(smooth(diff(time(startIdx:end)),3) < 5, 1, 'first');
        if isempty(altStart)
            continue
        end
        startIdx = startIdx + altStart-1;
        altEnd = find(diff(time(startIdx:end)) > 10, 1, 'first');
        endIdx = min([altEnd+startIdx-1,endIdx]);
        relTime = time(startIdx:endIdx)-time(startIdx);
        displ = displ(startIdx:endIdx);
        xCur = xCur(startIdx:endIdx);
        yCur = yCur(startIdx:endIdx);
        stepSizes = sqrt(diff(xCur).^2 + diff(yCur).^2);
        stepVel = stepSizes ./ diff(relTime);
        if any(stepVel > 0.5)
            plot(relTime(2:end),stepSizes)
            fprintf('Did the stage move for site %d cell %d?\n',k,i)
        end
        speedInterp = 0:2:max(relTime);
        xSpeed = interp1q(relTime,xCur,speedInterp');
        ySpeed = interp1q(relTime,yCur,speedInterp');
        altSpeed = mean(sqrt(diff(xCur).^2+diff(yCur).^2)./diff(relTime)); %only including original data points
        origSpeed = mean(sqrt(diff(xSpeed).^2+diff(ySpeed).^2)/2); % with interpolation
        allSpeeds{k} = vertcat(allSpeeds{k},origSpeed);
        allSpeedsNew{k}(i) = altSpeed;
        for j = 2:length(xCur)
            curSpeed = sqrt((xCur(j)-xCur(j-1))^2+(yCur(j)-yCur(j-1))^2)/(relTime(j)-relTime(j-1));
            allSpeedsNewSeparate{k} = vertcat(allSpeedsNewSeparate{k},curSpeed);
        end
        areaData = spreadData(k).area{i};
        if length(areaData) ~= length(time)
            error('Area data doesn''t match')
        end
        midArea1 = 40;
        midArea2 = 90;
        if areaData(1) < midArea1 && any(areaData > midArea2)
            startMid = find(areaData > midArea1, 1, 'first');
            endMid = find(areaData > midArea2,1,'first');
%             midIdx = max([startMid,startIdx]):min([endMid,endIdx]);
            midIdx = startMid:endMid;
            xMid = curMat(midIdx,4);
            yMid = curMat(midIdx,5);
            totDistCur = sum(sqrt(diff(xMid).^2 + diff(yMid).^2));
            totDist{k} = [totDist{k},totDistCur];
        end
        
        % calculation of time AC with no interpolation
        maxIdx = round(max(relTime)/2);
        AC = cell(1,maxIdx);
        for j = 1:length(xCur)-2
            curAngle = atan2(yCur(j+1)-yCur(j),xCur(j+1)-xCur(j));
            for m = j+1:length(xCur)-1
                curGap = relTime(m)-relTime(j);
                curIdx = round(curGap/2)+1;
                angle = atan2(yCur(m+1)-yCur(m),xCur(m+1)-xCur(m));
                AC{curIdx} = horzcat(AC{curIdx},cos(angle-curAngle));
            end
        end
        
        meanAC = ones(1,length(AC));
        for j = 2:length(AC)
            meanAC(j) = mean(AC{j});  
        end
        timeAC = 0:2:2*(length(meanAC)-1);
        
%         angle = atan2(diff(xInterp),diff(yInterp));
%         timeAC = 0:.1:max(relTime);
%         AC = ones(1,length(timeAC));
%         AC(end-2:end-1) = 0;
%         for j = 2:length(AC)-2
%             shiftIdx = find(abs(timeInterp-timeAC(j))<max(relTime)/10000);
%             shiftedAng = angle(shiftIdx:end);
%             curAng = angle(1:length(shiftedAng));
%             curTime = timeInterp(1:length(shiftedAng));
%             AC(j) = trapz(curTime,cos(shiftedAng-curAng))/(curTime(end));
%         end
%     
%         goodFit = false;
%         while ~goodFit
%             maxPlotIdx = find(timeAC < max(relTime)/10,1,'last');
%             plotRange = 1:maxPlotIdx;
% %             semilogy(timeAC(plotRange),AC(plotRange))
% %             hold on
%             %         [lengthStart,~] = selectPoint(lengthAC(plotRange),AC(plotRange));
%             %         [lengthEnd,~] = selectPoint(lengthAC(plotRange),AC(plotRange));
%             %         minFitIdx = find(abs(lengthStart-lengthAC)<maxLength/10000);
%             %         maxFitIdx = find(abs(lengthEnd-lengthAC)<maxLength/10000);
%             %         fitRange = minFitIdx:maxFitIdx;
%             fitRange = 1:round(maxPlotIdx/5);
%             curFit = polyfit(timeAC(fitRange),log(AC(fitRange)),1);
% %             semilogy(timeAC(fitRange),exp(polyval(curFit,timeAC(fitRange))))
% %             hold off
% %             choice = questdlg('Good fit or try again?','Good fit?','Good fit','Try again','Good fit');
% %             if strcmp(choice,'Good fit')
% %                 goodFit = true;
% %             end
%               goodFit = true;
%         end
        displStruct(k).timeACRevised{i} = meanAC;
        for j = 2:50
            if length(meanAC) < j
                break
            end
            displStruct(k).allTimeAC{j} = horzcat(displStruct(k).allTimeAC{j},AC{j});
        end
%         displStruct(k).timeACLogFit{i} = curFit;
        subplot(1,2,1)
        plot(xCur,yCur)
        daspect([1 1 1])
        title(sprintf('Site %d cell %d',k,i))
        subplot(1,2,2)
        plot(timeAC,meanAC)
        hold on
        smoothedAC = smooth(meanAC,7);
        plot(timeAC,smoothedAC)
        plot([timeAC(1),timeAC(end)],[.4 .4])
%         intIdx = 1;
%         for j = 2:length(AC)
%             if meanAC(j) < 0
%                 break
%             elseif length(AC{j}) > 10
%                 intIdx = j;
%             end
%         end
        intIdx = 8;
        if length(AC) >= intIdx
            intTime = timeAC(intIdx);
            if length(AC{intIdx}) > 5
                AUC = trapz(timeAC(1:intIdx),meanAC(1:intIdx));
                AUCVals{k}(i) = AUC;
                pDist{k}(i) = AUC * allSpeedsNew{k}(i);
                myfun = @(tau,intTime,AUC) tau.*(1 - exp(-intTime./tau)) - AUC;
                try tauCur = fzero(@(tau) myfun(tau,intTime,AUC),[0,intTime]);
                    pTime{k} = [pTime{k},tauCur];
                catch
                    %nothing
                end
            end
        end
        xlim([0 30])
        hold off
%         title('Press enter to go to next')
%         f = gcf;
%         waitforbuttonpress
%         val=get(f,'CurrentKey');
%         while ~strcmp(val,'return')
%             waitforbuttonpress
%             val=get(f,'CurrentKey');
%         end
%         hold off
    end
end

%% mean time AC
timeVals = 0:2:50;
numAC = length(timeVals);
figure
for k = 1:4
    ACCell = cell(numAC,1);
    for i = 1:length(displStruct(k).timeACRevised)
        curAC = displStruct(k).timeACRevised{i};
        if isempty(curAC)
            continue
        end
        maxIdx = min(numAC,length(curAC));
        for j = 1:maxIdx
            ACCell{j} = horzcat(ACCell{j},curAC(j));
        end
    end
    for i = 1:numAC
        meanIdx = find(~isnan(ACCell{i}));
        displStruct(k).meanTimeAC(i) = mean(ACCell{i}(meanIdx));
        displStruct(k).stdTimeAC(i) = std(ACCell{i}(meanIdx))/sqrt(length(ACCell{i}(meanIdx)));
    end
    idx = 1:numAC;
%     curLine = errorbar(timeVals(idx),displStruct(k).meanTimeAC(idx),displStruct(k).stdTimeAC(idx));
    for i = 1:numAC
        meanTimeAC(i) = mean(displStruct(k).allTimeAC{i});
        stdTimeAC(i) = std(displStruct(k).allTimeAC{i})/sqrt(length(displStruct(k).allTimeAC{i}));
    end
    curLine = errorbar(timeVals,meanTimeAC,stdTimeAC);
%     set(curLine,'LineStyle','none')
    set(curLine,'LineWidth',2)
    set(curLine,'Marker','+')
    set(curLine,'MarkerSize',10)
%     set(curLine,'Color',colorCell{k})
    hold on
%     curFit = polyfit(timeVals,log(displStruct(k).meanTimeAC),1);
%     plot(timeVals,exp(polyval(curFit,timeVals)))
    displStruct(k).alignedTimeAC = ACCell;
end
hold off
xlim([2 50])
ylim([0 .8])
legend('100%','10%','1%','0.1%')
    
%% Autocorrelation fitting
fitSlope = zeros(1,4);
logFit = cell(1,4);
figure
for k = 1:4
    curAC = displStruct(k).avgAutoCorr;
    time = 0:2:100;
    fitSlope(k) = time'\log(curAC');
    logFit{k} = polyfit(time,log(curAC),1);
    curLine = plot(time,curAC,'*');
    color = get(curLine,'Color');
    hold on
%     plot(time,exp(time.*fitSlope(k)),'Color',color)
%     plot(time,exp(polyval(logFit{k},time)))
end
% legend('100%','100% Fit','10%','10% Fit','1%','1% Fit','0.1%','0.1% Fit')
legend('100%','10%','1%','0.1%')

%% ANOVA testing on MSD
timeVal = 50;
dataVec = [];
groupCell = {};
groupVals = [100,10,1,0.1];
groupVec = [];
timeVec = 0:2:100;
for k = 1:4
%     curDisplCell = displStruct(k).displMat;
%     for j = 2:length(displStruct(k).alignedDispl)
%         curDisplData = sqrt(displStruct(k).alignedDispl{j})./timeVec(j);
%         curDisplData = sqrt(displStruct(k).alignedDispl{timeVal/2 + 1});
        curDisplData = sqrt(displStruct(k).alignedTimeGapDispl{timeVal/2});
        for i = 1:length(curDisplData)
            %         curTime = curDisplCell{i}(:,1);
            %         curTime = curTime - curTime(1);
            %         curDispl = curDisplCell{i}(:,3);
            %         if 1%any(curTime > timeVal-2.1 & curTime < timeVal + 2.1)
            %             [~,curIdx] = min(abs(curTime-timeVal));
            dataVec = vertcat(dataVec,curDisplData(i));
            groupCell = vertcat(groupCell,{sprintf('Group %d',k)});
            groupVec = vertcat(groupVec,groupVals(k));
            %         end
        end
%     end
end
[p, tbl, stats] = anova1(dataVec, groupCell);
figure
[compareGroups, means] = multcompare(stats);
[rho, pval] = corr(groupVec, dataVec, 'Type', 'Spearman');

%% ANOVA on autocorrelation
timeVal = 8;
timeIdx = round(timeVal/2) + 1;
dataVec = [];
groupCell = {};
groupVals = [100,10,1,0.1];
groupVec = [];
for k = 1:4
%     curACCell = displStruct(k).alignedTimeAC;
%     curGroup = persistL{k};
    curGroup = AUCVals{k};
%     curGroup = curACCell{timeIdx};
%     curGroup = displStruct(k).decayLength;
%     curGroup = displStruct(k).allTimeAC{11};
%     curGroup = allSpeedsNewSeparate{k};
    for i = 1:length(curGroup)
        if curGroup(i) == 0
            continue
        end
        dataVec = vertcat(dataVec,curGroup(i));
%         dataVec = vertcat(dataVec, displStruct(k).AC{i}(121));
        groupCell = vertcat(groupCell,{sprintf('Group %d',k)});
        groupVec = vertcat(groupVec,groupVals(k));
    end
end
[p, tbl, stats] = anova1(dataVec, groupCell);
figure
[compareGroups, means] = multcompare(stats);
[rho, pval] = corr(groupVec, dataVec, 'Type', 'Spearman');

%% text files for Volkmar
path = uigetdir;
cd(path)
siteString = {'100R','10R','1R','point1R'};
for k = 1:4
    mkdir(siteString{k})
    for i = 1:length(displStruct(k).displMat)
        curMat = displStruct(k).displMat{i};
        tStart = displStruct(k).startTimes(i);
        tEnd = displStruct(k).endTimes(i);
        time = curMat(:,1);
        xCur = curMat(:,4);
        yCur = curMat(:,5);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        relTime = time(startIdx:endIdx)-time(startIdx);
        xCur = xCur(startIdx:endIdx);
        yCur = yCur(startIdx:endIdx);
        fileName = sprintf('%s_%d.txt',siteString{k},i);
        fid = fopen(fullfile(siteString{k},fileName),'w');
        fprintf(fid,'%.4f\t%.4f\t%.4f\r\n',[relTime,xCur,yCur]');
        fclose(fid);
    end
end

%% Test displacement as a function of time gap
figure
for k = 1:4
    displStruct(k).timeGaps = cell(1,length(displStruct(k).displMat));
    displStruct(k).gapDispl = cell(1,length(displStruct(k).displMat));
    subplot(2,2,k)
    for i = 1:length(displStruct(k).displMat)
        curMat = displStruct(k).displMat{i};
        tStart = displStruct(k).startTimes(i);
        tEnd = displStruct(k).endTimes(i);
        if isnan(tStart) || isnan(tEnd)
            continue
        end
        time = curMat(:,1);
        displ = curMat(:,3);
        xCur = curMat(:,4);
        yCur = curMat(:,5);
        [~,startIdx] = min(abs(time-tStart));
        [~,endIdx] = min(abs(time-tEnd));
        relTime = time(startIdx:endIdx)-time(startIdx);
        displ = displ(startIdx:endIdx);
        xCur = xCur(startIdx:endIdx);
        yCur = yCur(startIdx:endIdx);
        
        maxTime = relTime(end);
        interpTime = (0:0.1:maxTime)';
        interpX = interp1q(relTime,xCur,interpTime);
        interpY = interp1q(relTime,yCur,interpTime);
        arcLength = zeros(1,length(interpX));
        for j = 2:length(interpX)
            arcLength(j) = arcLength(j-1) + sqrt((interpX(j)-interpX(j-1))^2 + ...
                (interpY(j)-interpY(j-1))^2);
        end
        timeGaps = 0:2:maxTime;
        gapDispl = zeros(1,length(timeGaps));
        xGapDispl = zeros(1,length(timeGaps));
        yGapDispl = zeros(1,length(timeGaps));
        totGapDispl = zeros(1,length(timeGaps));
        allRawMSD = cell(1,length(timeGaps));
        allRawMSD{1} = [0];
        gapDispl(end) = arcLength(20*(length(timeGaps)-1)+1)-arcLength(1);
        xGapDispl(end) = (interpX(20*(length(timeGaps)-1)+1)-interpX(1))^2;
        yGapDispl(end) = (interpY(20*(length(timeGaps)-1)+1)-interpY(1))^2;
        totGapDispl(end) = xGapDispl(end) + yGapDispl(end);
        for j = 2:length(timeGaps)-1
            shiftIdx = find(abs(interpTime-timeGaps(j))<.01);
            shiftedArcLength = arcLength(shiftIdx:end);
            xShifted = interpX(shiftIdx:end);
            yShifted = interpY(shiftIdx:end);
            curArcLength = arcLength(1:length(shiftedArcLength));
            xVals = interpX(1:length(xShifted));
            yVals = interpY(1:length(yShifted));
            curTimeGaps = interpTime(1:length(shiftedArcLength));
            gapDispl(j) = trapz(curTimeGaps,shiftedArcLength-curArcLength)/(curTimeGaps(end));
            xGapDispl(j) = trapz(curTimeGaps,(xShifted-xVals).^2)/curTimeGaps(end);
            yGapDispl(j) = trapz(curTimeGaps,(yShifted-yVals).^2)/curTimeGaps(end);
            totGapDispl(j) = trapz(curTimeGaps,sqrt((xShifted-xVals).^2 + (yShifted-yVals).^2)./curTimeGaps(end));
        end
        for j=1:length(relTime)-1
            for m = j+1:length(relTime)
                timeDiff = relTime(m)-relTime(j);
                if mod(timeDiff,2) > 0.5 && mod(timeDiff,2) < 1.5
                    continue
                else
                    curIdx = round(timeDiff/2) + 1;
                    if curIdx > length(allRawMSD)
                        continue
                    else
                        curMSD = (xCur(m)-xCur(j))^2 + (yCur(m)-yCur(j))^2;
                        allRawMSD{curIdx} = [allRawMSD{curIdx},curMSD];
                    end
                end
            end
        end
        displStruct(k).timeGaps{i} = timeGaps;
        displStruct(k).gapDispl{i} = gapDispl;
        displStruct(k).xGapDispl{i} = xGapDispl;
        displStruct(k).yGapDispl{i} = yGapDispl;
        displStruct(k).gapMSD{i} = xGapDispl + yGapDispl;
        displStruct(k).totGapDispl{i} = totGapDispl;
        displStruct(k).allRawMSD{i} = allRawMSD;
        plot(timeGaps,gapDispl)
        hold on
        drawnow
    end
    uiwait(msgbox(sprintf('Surface %d, cell %d',k,i)))
    hold off
end

%% mean time gap displacement curve
timeGaps = 0:2:100;
numPoints = length(timeGaps);
figure
for k = 1:4
    GapDisplMat = zeros(length(displStruct(k).gapDispl),numPoints);
    curData = cell(1,numPoints);
    curXData = cell(1,numPoints);
    curYData = cell(1,numPoints);
    curMSD = cell(1,numPoints);
    curTotData = cell(1,numPoints);
    rawMSD = cell(1,numPoints);
    totIt = length(displStruct(k).timeGaps);
    if k == 3
        totIt = 10;
    end
    for i = 1:totIt
        if length(displStruct(k).gapDispl{i}) < numPoints
            L = length(displStruct(k).timeGaps{i});
        else
            L = numPoints;
        end
        for j = 1:L
            curData{j} = [curData{j},displStruct(k).gapDispl{i}(j)];
            curXData{j} = [curXData{j},displStruct(k).xGapDispl{i}(j)];
            curYData{j} = [curYData{j},displStruct(k).yGapDispl{i}(j)];
            curMSD{j} = [curMSD{j},displStruct(k).gapMSD{i}(j)];
            curTotData{j} = [curTotData{j},displStruct(k).totGapDispl{i}(j)];
            rawMSD{j} = [rawMSD{j},displStruct(k).allRawMSD{i}{j}];
        end
    end
    for i = 1:numPoints
        displStruct(k).meanGapDispl(i) = mean(curData{i});
        displStruct(k).stdGapDispl(i) = std(curData{i})/sqrt(length(curData{i}));
        displStruct(k).meanXGapDispl(i) = mean(curXData{i});
        displStruct(k).stdXGapDispl(i) = std(curXData{i})./sqrt(length(curXData{i}));
        displStruct(k).meanYGapDispl(i) = mean(curYData{i});
        displStruct(k).stdYGapDispl(i) = std(curYData{i})./sqrt(length(curYData{i}));
        displStruct(k).meanGapMSD(i) = mean(curMSD{i});
        displStruct(k).stdGapMSD(i) = std(curMSD{i})./sqrt(length(curMSD{i}));
        displStruct(k).meanTotGapDispl(i) = mean(curTotData{i});
        displStruct(k).stdTotGapDispl(i) = std(curTotData{i})./sqrt(length(curTotData{i}));
        displStruct(k).meanRawMSD(i) = mean(rawMSD{i});
        displStruct(k).stdRawMSD(i) = std(rawMSD{i})./sqrt(length(rawMSD{i}));
    end
%     plot(lengthVals,displStruct(k).meanLengthAC)
    subplot(1,2,1)
    curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanXGapDispl(1:2:end)...
        ,displStruct(k).stdXGapDispl(1:2:end));
    title('x MSD vs time gap')
    xlabel('time gap (s)')
    ylabel('x MSD (\mum^2)')
    set(curLine,'LineWidth',2)
    set(curLine,'Color',colorCell{k})
    hold on
    subplot(1,2,2)
    curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanYGapDispl(1:2:end)...
        ,displStruct(k).stdYGapDispl(1:2:end));
    title('y MSD vs time gap')
    xlabel('time gap (s)')
    ylabel('y MSD (\mum^2)')
    set(curLine,'LineWidth',2)
    set(curLine,'Color',colorCell{k})
    hold on
%     curFit = polyfit(timeVals,log(displStruct(k).meanTimeAC),1);
%     plot(timeVals,exp(polyval(curFit,timeVals)))
    displStruct(k).alignedTimeGapDispl = curTotData;
end
hold off
subplot(1,2,1)
legend('100%','10%','1%','0.1%')
subplot(1,2,2)
legend('100%','10%','1%','0.1%')

%% mean gap displ
figure
expType = fittype('a.*t.^b','independent','t');
fitMSD = zeros(4,2);
for k = 1:4
%     curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanGapMSD(1:2:end)...
%        ,displStruct(k).stdGapMSD(1:2:end));
%    curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanTotGapDispl(1:2:end)...
%        ,displStruct(k).stdTotGapDispl(1:2:end));
%     curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanGapDispl(1:2:end)...
%        ,displStruct(k).stdGapDispl(1:2:end));
    curLine = errorbar(timeGaps(1:2:end),displStruct(k).meanRawMSD(1:2:end)...
        ,displStruct(k).stdRawMSD(1:2:end),'LineStyle','None','LineWidth',1);
    hold on
%     title('Mean squared displacement vs time gap')
%     xlabel('time gap (s)')
%     ylabel('MSD (\mum^2)')
%     plot(timeGaps,displStruct(k).meanRawMSD,'LineWidth',2,'Color',colorCell{k})
    fitMSD(k,:) = polyfit(log(timeGaps(3:31)),log(displStruct(k).meanRawMSD(3:31)),1);
%     curFit = fit(timeGaps(3:31)',displStruct(k).meanRawMSD(3:31)',expType);
%     fitMSD(k,1) = curFit.b;
    plot(timeGaps,exp(fitMSD(k,2))*timeGaps.^(fitMSD(k,1)),'LineWidth',2,'Color',colorCell{k})
%     plot(timeGaps,curFit(timeGaps),'LineWidth',2,'Color',colorCell{k})
    set(curLine,'Color',colorCell{k})
    hold on
end
legend('100%','10%','1%','0.1%')
xlim([0 60])

%% print text files for raw msd
numPts = 0;
for i = 1:length(rawMSD)
    if length(rawMSD{i}) > numPts
        numPts = length(rawMSD{i});
    end
end
[MSDFile,MSDPath] = uiputfile('*.txt','Location to save raw MSD');
fid = fopen(fullfile(MSDPath,MSDFile),'w');
for i = 1:numPts
    for j = 1:length(rawMSD)
        if length(rawMSD{j}) < i
            fprintf(fid,'--');
        else
            fprintf(fid,'%.8f',rawMSD{j}(i));
        end
        if j == length(rawMSD)
            fprintf(fid,'\r\n');
        else
            fprintf(fid,'\t');
        end
    end
end
fclose(fid);