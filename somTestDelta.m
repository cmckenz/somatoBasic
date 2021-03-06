% somdis.m
%
%        $Id:$ 
%      usage: somTestDelta(delta)
%         by: cam mckenzie
%       date: dec 16
%    purpose: somFixedDelta with only one block for experimental induction
%
%               Observer discriminates which of two temporal intervals
%               contains an increment over the base stimulus intensity.
%               Two piezo buzzers are used, only one of which has increment
%               applied. If attended condition, arrow indicates left (-1)/
%               right (1) buzzer has increment, else both equally likely.
% 
%           This takes a fixed delta and uses for all conditions
%
%               Start by running somGetDelta to get delta for middle
%               pedestal
%
function myscreen = somTestDelta(delta)

%this is the correct setting for dubonnet... (now redundant)
% system('osascript -e "set volume 6"') USE VOLUME SETTING IN SCREEN
% SETTINGS: (MGLEDITSCREENPARAMS)

% initalize the screen (uncomment whichever you want)
myscreen = initScreen('somato');
%myscreen = initScreen('test');


% task parameters
task{1}.waitForBacktick = 0;
task{1}.segmin = [1 0.5 0.2 0.5 0.3 2 1]; %attention cue, stimulus interval 1, stim off, stimulus interval 2, pause, response cue, feedback
task{1}.segmax = [1 0.5 0.2 0.5 0.3 2 1.5];
task{1}.synchToVol = [0 0 0 0 0 0 0];
task{1}.getResponse = [0 0 0 1 1 1 0];
%task{1}.parameter.pedestal = [0 0.125 0.25 0.5]; old parameters
task{1}.parameter.pedestal = [0.2 0.45 0.7];
task{1}.randVars.block.distractPed = task{1}.parameter.pedestal;
task{1}.parameter.side = [-1 1];
task{1}.parameter.attention = [0 1];
task{1}.randVars.uniform.interval = [0 1];
task{1}.randVars.calculated.correct = nan;
task{1}.randVars.calculated.responseInterval = nan;
task{1}.freq = 80;
task{1}.onTime = 0.5;
task{1}.waitTime = 0.2;
task{1}.random = 1;
task{1}.numBlocks = 1;





% initialize the task
for phaseNum = 1:length(task)
  [task{phaseNum} myscreen] = initTask(task{phaseNum},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);

stimfile = [];
stimulus = myInitStimulus(stimulus,myscreen,task,stimfile, delta);


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

    %pedestal
    ped = task.thistrial.pedestal;

    % stimulation level on distractor side?
    distractPed = task.thistrial.distractPed;

  % add the delta to the correct location
  if task.thistrial.side == -1
    stimulus.stimLeft = ped + stimulus.delta;
    stimulus.stimRight = distractPed;
    stimulus.foilLeft = ped;
    stimulus.foilRight = distractPed;
  else
    stimulus.stimLeft = distractPed;
    stimulus.stimRight = ped + stimulus.delta;
    stimulus.foilLeft = distractPed;
    stimulus.foilRight = ped;
  end
  
  % make sure we do not go over the maximum amplitude
  if stimulus.stimLeft > 1
      disp('WARNING: LEFT STIMULUS GREATER THAN 1 --> REDUCING STIMULUS INTENSITY');
      disp(sprintf('L Stimulus of %.2f reduced to 1', stimulus.stimLeft));
      stimulus.stimLeft = 1;
  end
  if stimulus.stimRight > 1
      disp('WARNING: RIGHT STIMULUS GREATER THAN 1 --> REDUCING STIMULUS INTENSITY');
      disp(sprintf('R Stimulus of %.2f reduced to 1', stimulus.stimRight));
      stimulus.stimRight = 1;
  end
  % multiply to get the actual voltage amplitude
  stimulus.stimLeft = stimulus.stimLeft * stimulus.maxStim;
  stimulus.stimRight = stimulus.stimRight * stimulus.maxStim;
  
  pulseWave = [stimulus.stimBase(1,:)*stimulus.stimLeft; stimulus.stimBase(2,:)*stimulus.stimRight];
  foilWave = [stimulus.stimBase(1,:)*stimulus.foilLeft; stimulus.stimBase(2,:)*stimulus.foilRight];
    
  
  stimulus.stimInd = mglInstallSound(pulseWave);
  stimulus.foilInd = mglInstallSound(foilWave);
  
  mglSetSound(stimulus.stimInd, 'deviceID', stimulus.deviceID);
  mglSetSound(stimulus.foilInd, 'deviceID', stimulus.deviceID);
end

if task.thistrial.thisseg == 2
    if task.thistrial.interval == 0 %play the stimulus (stim + ped on one side) if first interval is stimulus interval
        mglPlaySound(stimulus.stimInd);
    else
        mglPlaySound(stimulus.foilInd); %else play the foil (just pedestals)
    end
  % set up the two stimulation period
  
  % display what we are doing
  disp(sprintf('Trial %i: %0.3f vs %0.3f (pedestal: %0.3f distractor: %0.3f delta: %0.3f side: %i interval: %i)',task.trialnum,stimulus.stimLeft,stimulus.stimRight,task.thistrial.pedestal,task.thistrial.distractPed, stimulus.delta,task.thistrial.side, task.thistrial.interval));
  % store delta
  task.thistrial.delta = stimulus.delta;
end

if task.thistrial.thisseg == 4
    if task.thistrial.interval == 1 %play the stimulus if second interval (seg 4) is stimulus interval
        mglPlaySound(stimulus.stimInd);
    else
        mglPlaySound(stimulus.foilInd); %else play the foil (just pedestals)
    end
  % set up the two stimulation period
  
  % display what we are doing
  disp('Delivering stimulus for interval two.');
  % store delta
  task.thistrial.delta = stimulus.delta;
end

if task.thistrial.thisseg == 7
  % no response
  if isnan(task.thistrial.correct)
    task.thistrial.correct = false;

    disp(sprintf('  No response. Delta: %0.2f',task.thistrial.delta));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback(task, myscreen)

global stimulus

% display attention cue
mglClearScreen;

if any(task.thistrial.thisseg == [1 5 6])
    mglFixationCross(0.75, 0.75, [1 1 1]); % White cross for pre-stim and response intervals
end

%Stimulation interval: green fixation cross

if any(task.thistrial.thisseg == [2 4]);
    mglFixationCross(0.75, 0.75, [0 1 0]);
end


if any(task.thistrial.thisseg == [1 2 3 4])
    
    if task.thistrial.attention == 0;
        leftCue = [1 1 1];
        rightCue = [1 1 1];
    end

if task.thistrial.attention == 1;
    if task.thistrial.side == -1;
        leftCue = [1 1 1];
        rightCue = [0 0 0];
    end
    if task.thistrial.side == 1;
        leftCue = [0 0 0];
        rightCue = [1 1 1];
    end
end

%where to draw the left arrow
  x = [3     4.5    3];
  y = [0.75  0   -0.75];
  
  mglPolygon(-x, y, leftCue);
  mglPolygon(x, y, rightCue);
end

%3rd segment: blank period

if task.thistrial.thisseg == 6 %response cue... where was the stimulus?

%  mglFixationCross(1,1,[0 1 1]);
  if task.thistrial.side == -1
    leftCue = [1 1 1];
    rightCue = [0 0 0];
  end
  if task.thistrial.side == 1
    leftCue = [0 0 0];
    rightCue = [1 1 1];
  end
  
  %where to draw the left arrow
  x = [3     4.5    3];
  y = [0.75  0   -0.75];
  
  
    % set colors only after the subject has responded
  if ~isnan(task.thistrial.correct)
      if task.thistrial.correct
          leftCue = leftCue.*[0 1 0];
          rightCue = rightCue.*[0 1 0];
      else
          leftCue = leftCue.*[1 0 0];
          rightCue = rightCue.*[1 0 0];
      end
  end
  
  mglPolygon(-x, y, leftCue);
  mglPolygon(x, y, rightCue);

end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback(task,myscreen)

global stimulus

if task.thistrial.gotResponse < 1
  % caluclate response
  task.thistrial.responseInterval = (task.thistrial.whichButton-1);
  if task.thistrial.responseInterval == task.thistrial.interval
    task.thistrial.correct = true;
  else
    task.thistrial.correct = false;
  end
  
  
  disp(sprintf('  Response: %i Correct: %i Delta: %0.2f',task.thistrial.responseInterval,task.thistrial.correct,task.thistrial.delta));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,stimfile, delta)

% set the maximum amplitude
stimulus.maxStim = 1;

%Basic stimulus waveform - Sine Wave, Amplitude = 1 (or = 2 peak to peak),
% Offset from beginning of file by 0.01s
onTime = task{1}.onTime
waitTime = task{1}.waitTime
freq = task{1}.freq;
sampleRate = 8192;

stimBase = sin(freq * (2*pi*(1:round(onTime*sampleRate))/sampleRate));

stimBase = [stimBase; stimBase];

% waitTime = 0.2; %Wait time now set as length of 3rd segment - no
% associated audio waveform

% waitBase = zeros(2, round(waitTime*sampleRate));

stimulus.stimBase = stimBase;
% stimulus.waitBase = waitBase;

stimulus.deviceID = 1; %deviceID for sound

stimulus.delta = delta;
