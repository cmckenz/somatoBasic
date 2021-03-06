% somdis.m
%
%        $Id:$ 
%      usage: somdisNS()
%         by: justin gardner/cameron mckenzie
%       date: 06/03/14
%    purpose: Basic somato discrim task - removed staircase, changed timing for scanner
%
function myscreen = somdisNS()

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen
testing = false;


myscreen = initScreen;

% get the last stimfile
stimfile = getLastStimfile(myscreen,'stimfileNum=-1');

% task parameters
task{1}.waitForBacktick = 1;
task{1}.segmin = [0.25 0.75 2 2 2]; %added a 250ms initial segment to set up stim sound
task{1}.segmax = [0.25 0.75 6 2 6]; %set stim/play stim/response delay/response interval/inter-trial interval
task{1}.synchToVol = [1 0 1 0 0]; %synch before stimulus, responses
task{1}.getResponse = [0 0 0 1 0];
task{1}.parameter.pedestal = [0.5 0.25 0.125 0];
task{1}.parameter.stimulus = [1 2 3];
task{1}.randVars.uniform.side = [-1 1];
task{1}.randVars.calculated.threshold = nan;
task{1}.randVars.calculated.correct = nan;
task{1}.randVars.calculated.responseSide = nan;
task{1}.random = 1;
task{1}.numtrials = 250;

if testing
  task{1}.waitForBacktick = 0;
  task{1}.segmin = [0.25 1 2 1 0.25];
  task{1}.segmax = [0.25 1 2 1 0.25];
  task{1}.parameter.pedestal = [-1 0.5 0.25 0.125 0];
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
myscreen.feedback = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

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
  % time, so that we can preciesly set the two stimulation intervals
  timeNow = mglGetSecs;
  % pedestal stimulus strength
  ped = task.thistrial.pedestal;
  stimNum = task.thistrial.stimulus;
  % delta 
  pedNum = find(stimulus.pedestal==ped);
  
  stimulus.delta = stimulus.stimVal(pedNum, stimNum);
  
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
  
  stimWave = [stimulus.stimBase(1,:)*stimulus.stimLeft; stimulus.stimBase(2,:)*stimulus.stimRight];
  stimulus.stimInd = mglInstallSound(stimWave);
  mglSetSound(stimulus.stimInd, 'deviceID', stimulus.deviceID);
end

if task.thistrial.thisseg == 2
  % set up the two stimulation period
  mglPlaySound(stimulus.stimInd);
  % display what we are doing
  disp(sprintf('Trial %i: %0.3f vs %0.3f (pedestal: %0.2f delta: %0.2f side: %i)',task.trialnum,stimulus.stimLeft,stimulus.stimRight,task.thistrial.pedestal,stimulus.delta,task.thistrial.side));
  % store delta
  task.thistrial.delta = stimulus.delta;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% display fixation cross
mglClearScreen;

if myscreen.feedback; %set feedback flag in myscreen to true for feedback


if any(task.thistrial.thisseg == [1 2 4])
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

end

if ~ myscreen.feedback; %feedback flag set to false?
    if any(task.thistrial.thisseg == [1 2])
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
%Display what is happening with the task
  disp(sprintf('  Response: %i Correct: %i Delta: %0.2f',task.thistrial.responseSide,task.thistrial.correct,task.thistrial.delta));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,stimfile)

% set the maximum amplitude
stimulus.maxStim = 1;

%Basic stimulus waveform - Sine Wave, Amplitude = 1 (or = 2 peak to peak),
% Offset from beginning of file by 0.1s
startTime = 0.01;
eventLength = 0.4;
freq = 80;
sampleRate = 8192;

somEvents{1}(1) = startTime;
somEvents{1}(2) = eventLength;
somEvents{1}(3) = freq;
somEvents{1}(4) = 0;

somEvents{2}(1) = startTime;
somEvents{2}(2) = eventLength;
somEvents{2}(3) = freq;
somEvents{2}(4) = 1;

stimBase = somTrial(somEvents, sampleRate);


% determine pedestals
params = getTaskParameters(myscreen,task);
pedestal = params.originalTaskParameter.pedestal;
numPedestal = length(pedestal);
stimulus.pedestal = pedestal;
stimulus.numPedestal = numPedestal;
stimulus.stimBase = stimBase;
stimulus.deviceID = 2;
stimfile = [];

threshold1 = 0.05;
threshold2 = 0.05;
threshold3 = 0.05;
threshold4 = 0.05;

stimulus.stimVal(1,1:3) = threshold1*[1 2 4];
stimulus.stimVal(2,1:3) = threshold2*[1 2 4];
stimulus.stimVal(3,1:3) = threshold3*[1 2 4];
stimulus.stimVal(4,1:3) = threshold4*[1 2 4];
  
  


