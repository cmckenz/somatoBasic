% somdisPS.m
%
%        $Id:$ 
%      usage: somdisPS()
%         by: justin gardner
%       date: 06/03/14
%    purpose: Basic somato discrim task
%
function myscreen = somdisFreq()

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen
testing = 1;

myscreen = initScreen;

% get the last stimfile
stimfile = getLastStimfile(myscreen,'stimfileNum=-1', 'onlyToday=true')

% task parameters
task{1}.waitForBacktick = 1;
task{1}.segmin = [0.25 0.75 2 1]; %added a 250ms initial segment to set up stim sound
task{1}.segmax = [0.25 0.75 2 9];
task{1}.synchToVol = [0 0 0 1];
task{1}.getResponse = [1 1 1 0];
task{1}.parameter.pedestal = 0;
task{1}.parameter.frequency = [20 60 120 180];
task{1}.randVars.uniform.side = [-1 1];
task{1}.randVars.calculated.threshold = nan;
task{1}.randVars.calculated.correct = nan;
task{1}.randVars.calculated.responseSide = nan;
task{1}.random = 1;
task{1}.numTrials = 250;

if testing
  task{1}.waitForBacktick = 0;
  task{1}.segmin = [0.25 1 2 1];
  task{1}.segmax = [0.25 1 2 1];
  task{1}.parameter.pedestal = 0;
  task{1}.parameter.frequency = [20 60 120 180];
%  task{1}.parameter.pedestal = [-1 0.5];
  task{1}.synchToVol = [0 0 0 0];
end




% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,stimfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task)) && ~myscreen.userHitEsc
  % update the task
  [task myscreen phaseNum] = updateTask(task,myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  % pedestal stimulus strength
  ped = task.thistrial.pedestal;
  % stimulus frequency
  freq = task.thistrial.frequency;
  % which staircase? which delta?
  staircaseNum = find(stimulus.frequency==freq);
  [stimulus.delta stimulus.s(staircaseNum)] = doStaircase('testValue',stimulus.s(staircaseNum));
  if (stimulus.delta < 0) stimulus.delta = 0;end
  
  % if this is ped == -1, then it means to give max difference randomly betwen sides
  if ped == -1
    ped = 0;
    stimulus.delta = 1;
  end
  
  % add the delta to the correct side
  if task.thistrial.side == -1
    stimulus.stimLeft = ped + stimulus.delta;
    stimulus.stimRight = ped;
  else
    stimulus.stimLeft = ped;
    stimulus.stimRight = ped + stimulus.delta;
  end
  % make sure we do not go over the maximum amplitude
  if stimulus.stimLeft > 1,stimulus.stimLeft = 1;end
  if stimulus.stimRight > 1,stimulus.stimRight = 1;end
  
  % multiply to get the actual voltage amplitude
  stimulus.stimLeft = stimulus.stimLeft * stimulus.maxStim;
  stimulus.stimRight = stimulus.stimRight * stimulus.maxStim;
  
  freqStim = stimulus.stimBase{staircaseNum};
  
  stimWave = [freqStim(1,:)*stimulus.stimLeft; freqStim(2,:)*stimulus.stimRight];
  stimulus.stimInd = mglInstallSound(stimWave);
  mglSetSound(stimulus.stimInd, 'deviceID', stimulus.deviceID);
end

if task.thistrial.thisseg == 2
  % set up the two stimulation period
  mglPlaySound(stimulus.stimInd);
  % display what we are doing
  disp(sprintf('Trial %i: %0.3f vs %0.3f (pedestal: %0.2f frequency: %i delta: %0.2f side: %i)',task.trialnum,stimulus.stimLeft,stimulus.stimRight,task.thistrial.pedestal,task.thistrial.frequency, stimulus.delta,task.thistrial.side));
  % store delta
  task.thistrial.delta = stimulus.delta;
end

if task.thistrial.thisseg == 4
  % no response
  if isnan(task.thistrial.correct)
    task.thistrial.correct = false;
    % update appropriate staircase
    staircaseNum = find(stimulus.frequency==task.thistrial.frequency);
    stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
    %  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
    disp(sprintf('  No response. Delta: %0.2f',task.thistrial.delta));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% display fixation cross
mglClearScreen;
if any(task.thistrial.thisseg == [1 2 3])
%  mglFixationCross(1,1,[0 1 1]);
  if task.thistrial.side == -1
    leftColor = [0 1 0];
    rightColor = [1 0 0];
  else
    leftColor = [1 0 0];
    rightColor = [0 1 0];
  end
  % set colors only after the subject has responded
  if ~isnan(task.thistrial.correct)
    % display the chosen side in green or red contingent
    % on whether subject was correct or not
    if task.thistrial.responseSide == -1
      mglFillOval(-4,0,[0.5 0.5],leftColor);
    else
      mglFillOval(-4,0,[0.5 0.5],[1 1 1]);
    end
    if task.thistrial.responseSide == 1
      mglFillOval(4,0,[0.5 0.5],rightColor);
    else
      mglFillOval(4,0,[0.5 0.5],[1 1 1]);
    end
  else
    mglFillOval(-4,0,[0.5 0.5],[1 1 1]);
    mglFillOval(4,0,[0.5 0.5],[1 1 1]);
  end

end  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if task.thistrial.gotResponse < 1
  % caluclate response
  task.thistrial.responseSide = (task.thistrial.whichButton-1.5)*2;
  if task.thistrial.responseSide == task.thistrial.side
    task.thistrial.correct = true;
  else
    task.thistrial.correct = false;
  end
  % update appropriate staircase
  staircaseNum = find(stimulus.frequency==task.thistrial.frequency);
  stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
%  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
  disp(sprintf('  Response: %i Correct: %i Delta: %0.2f',task.thistrial.responseSide,task.thistrial.correct,task.thistrial.delta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,stimfile)

% set the maximum amplitude
stimulus.maxStim = 1;

% determine frequencies
params = getTaskParameters(myscreen,task);
frequency = params.originalTaskParameter.frequency;
numFreqs = length(frequency);
stimulus.frequency = frequency;
stimulus.numFreqs = numFreqs;
stimulus.deviceID = 1;
stimfile = [];


%Basic stimulus waveform - Sine Wave, Amplitude = 1 (or = 2 peak to peak),
% Offset from beginning of file by 0.1s


startTime = 0.1;
eventLength = 0.4;
sampleRate = 8192;

stimulus.stimBase = {};

for ii = 1:numFreqs
    
    freq = stimulus.frequency(ii);

    somEvents{1}(1) = startTime;
    somEvents{1}(2) = eventLength;
    somEvents{1}(3) = freq;
    somEvents{1}(4) = 0;

    somEvents{2}(1) = startTime;
    somEvents{2}(2) = eventLength;
    somEvents{2}(3) = freq;
    somEvents{2}(4) = 1;

stimulus.stimBase{ii} = somTrial(somEvents, sampleRate);

end


% first time, initialize the staircases,
if isempty(stimfile)
  for iStaircase = 1:numFreqs
    stimulus.s(iStaircase) = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=0.3','initialStepsize=0.1','nTrials=40','stepRule=levitt','minThreshold=0');
  end
% subsequent times, continue staircases from where we left off
else
  oldStimulus = stimfile.stimulus;
  % make sure all is the same
  if ~isfield(oldStimulus,'frequency') || ~isfield(oldStimulus,'numFreqs') || ~isfield(oldStimulus,'s') || ~isequal(oldStimulus.frequency,frequency) || ~isequal(oldStimulus.numFreqs,numFreqs) || ~(length(oldStimulus.s) == numFreqs)
    disp(sprintf('(somdisFreq) !!! Previous stimfile has different conditions than this run, so restarting staircases !!!'));
    stimulus = myInitStimulus(stimulus,myscreen,task,[]);
    return
  end
  % if all is the same, then just reinit the staircases with the previous ones
  for iStaircase = 1:numFreqs
    stimulus.s(iStaircase) = doStaircase('init',oldStimulus.s(iStaircase));
  end
end

  
  


