%% Staircase nums: L/R L/R

day = '19'
month = '12'
year = '16'
subj = 's025'


files = dir(['~/data/somGetDelta/' subj '/' year month day '*.mat'])
%%
load(['~/data/somGetDelta/' subj '/' files(end).name])

leftStairs = doStaircase('combine',stimulus.s(1), stimulus.s(3))
rightStairs = doStaircase('combine',stimulus.s(2), stimulus.s(4))

lThresh = doStaircase('getThreshold', leftStairs)
rThresh = doStaircase('getThreshold', rightStairs)



index1 = find(lThresh.weibullFitParams(:,1).y < 0.9)
leftVal = lThresh.weibullFitParams(:,1).x(max(index1));

index2 = find(rThresh.weibullFitParams(:,1).y < 0.9)
rightVal = rThresh.weibullFitParams(:,1).x(max(index2));



grandMean = mean([leftVal, rightVal])