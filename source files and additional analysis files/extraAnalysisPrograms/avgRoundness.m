% average roundness over entire spreading window
[file,path] = uigetfile('*.txt','Find roundness and displacement data file');
dataMat = load(fullfile(path,file));
time = dataMat(:,1);
roundness = dataMat(:,2);
rangeInput = inputdlg({'Enter starting time for spreading','Enter ending time for spreading'});
startTime = str2double(rangeInput{1});
endTime = str2double(rangeInput{2});
[~,startIdx] = min(abs(time-startTime));
[~,endIdx] = min(abs(time-endTime));
avgRound = (1/(endTime-startTime)) * trapz(time(startIdx:endIdx),roundness(startIdx:endIdx));
plot(time,roundness)
hold on
plot(time(startIdx:endIdx),avgRound*ones(1,endIdx-startIdx+1))
title('Roundness vs time')
legend('Roundness','Average roundness')
xlabel('time (s)')
ylabel('Roundness')

%% already organized data...
window = 10;
for i = 1:4
    displStruct(i).avgRoundness = zeros(length(displStruct(i).displMat),1);
    displStruct(i).windowRoundness = zeros(length(displStruct(i).displMat),1);
    for j = 1:length(displStruct(i).displMat)
        time = displStruct(i).displMat{j}(:,1);
        roundness = displStruct(i).displMat{j}(:,2);
        startTime = displStruct(i).startTimes(j);
        endTime = displStruct(i).endTimes(j);
        if isnan(startTime) || isnan(endTime)
            continue
        else
            [~,startIdx] = min(abs(time-startTime));
            [~,endIdx] = min(abs(time-endTime));
            avgRound = (1/(endTime-startTime)) * trapz(time(startIdx:endIdx),roundness(startIdx:endIdx));
            displStruct(i).avgRoundness(j,1) = avgRound;
            midTime = (startTime+endTime)/2;
            startMid = midTime - window/2;
            endMid = midTime + window/2;
            midInterp = (startMid:.1:endMid)';
            roundInterp = interp1q(time,roundness,midInterp);
            windowRound = (1/(endMid-startMid)) * trapz(midInterp,roundInterp);
            displStruct(i).windowRoundness(j,1) = windowRound;
        end
    end
end