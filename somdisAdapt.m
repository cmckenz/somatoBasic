% somdis.m
%
%        $Id:$ 
%      usage: somdisAdapt()
%         by: cam mckenzie
%       date: jul/aug 16
%    purpose: Somato discrimination task with attentional manipulation
%
%               Tactile discrimination on one hand after adaptation.
%   Standard params: adapt at 0.2, 0.35 or 0.5
%                    pedestals at same



function myscreen = somdisAdapt()

%this is the correct setting for dubonnet... (now redundant)
% system('osascript -e "set volume 6"') USE VOLUME SETTING IN SCREEN
% SETTINGS: (MGLEDITSCREENPARAMS)

testing = 0;
adaptTime = 20;
adaptLevel = 0.2;
% initalize the screen
if testing == 1
    myscreen = initScreen('test');
else
    myscreen = initScreen('somato');
end



% if oldStairs
%     % get the last stimfile
%     stimfile = getLastStimfile(myscreen,'stimfileNum=-1');
% else
    stimfile = [];
% end
%task parameters - phase 1, adaptation

task{1}.waitForBacktick = 0;
task{1}.seglen = adaptTime + 5;
task{1}.getResponse = 0;
task{1}.synchToVol = 0;
task{1}.parameter.adaptLevel = adaptLevel;
task{1}.freq = 80;
task{1}.numTrials = 1;


% task parameters - phase 2, testing...
task{2}.waitForBacktick = 0;
task{2}.segmin = [1 0.5 0.2 0.5 0.3 2 1 5 0.5]; %fix cross, stimulus interval 1, stim off, stimulus interval 2, pause, response cue, feedback, adapt topup, pause
task{2}.segmax = [1 0.5 0.2 0.5 0.3 2 1.5 5 0.5]; %Adaptation topup: segment 8
task{2}.synchToVol = [0 0 0 0 0 0 0 0 0];
task{2}.getResponse = [0 0 0 1 1 1 0 0 0];
task{2}.parameter.pedestal = [0.2 0.35 0.5];
task{2}.parameter.adaptLevel = task{1}.parameter.adaptLevel;
task{2}.parameter.side = -1; %just a reminder
task{2}.parameter.interval = [0 1];
task{2}.randVars.calculated.threshold = nan;
task{2}.randVars.calculated.correct = nan;
task{2}.randVars.calculated.responseInterval = nan;
task{2}.freq = 80;
task{2}.onTime = 0.5;
task{2}.waitTime = 0.2;
task{2}.random = 1;
task{2}.numTrials = 90;





% initialize the task
[task{1} myscreen] = initTask(task{1},myscreen,@startSegmentAdapt,@screenUpdateCallback);
[task{2} myscreen] = initTask(task{2},myscreen,@startSegmentCallback,@screenUpdateCallback,@responseCallback);


% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task,stimfile);


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
% function that gets called for the adaptation block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task myscreen] = startSegmentAdapt(task, myscreen)

global stimulus

%Basic stimulus waveform - Sine Wave, Amplitude = 1 (or = 2 peak to peak),
% Offset from beginning of file by 0.01s
onTime = task.seglen - 5;
freq = task.freq;
sampleRate = 8192;

stimulus.adaptLevel = task.parameter.adaptLevel

adaptBase = stimulus.adaptLevel * sin(freq * (2*pi*(1:round(onTime*sampleRate))/sampleRate));
soundNum = mglInstallSound(adaptBase);
mglSetSound(soundNum, 'deviceID', stimulus.deviceID);
mglPlaySound(soundNum);
disp('Presenting adaptor stimulus');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback(task, myscreen)

global stimulus;

if task.thistrial.thisseg == 1
  
  % pedestal stimulus strength
  ped = task.thistrial.pedestal;
  % delta 
  
  % delta 
  staircaseNum = find(stimulus.pedestal==ped);
 

  [stimulus.delta stimulus.s(staircaseNum)] = doStaircase('testValue',stimulus.s(staircaseNum));
  if (stimulus.delta < 0) stimulus.delta = 0;end

  % make sure we do not go over the maximum amplitude
  if (ped + stimulus.delta) > 1
      disp('WARNING: STIMULUS GREATER THAN 1 --> REDUCING STIMULUS INTENSITY');
      disp(sprintf('L Stimulus of %.2f reduced to 1', (ped + stimulus.delta)));
      stimulus.delta = 1 - ped;
  end
  
  pulseWave = [stimulus.stimBase*(ped + stimulus.delta)];
  foilWave = [stimulus.stimBase*ped];
  
  % multiply to get the actual voltage amplitude
  pulseWave = pulseWave * stimulus.maxStim;
  foilWave = foilWave * stimulus.maxStim;
  
    
  
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
  % set up the two stimulation periods
  
  % display what we are doing
  disp(sprintf('Trial %i: (pedestal: %0.3f delta: %0.3f interval: %i)',task.trialnum, task.thistrial.pedestal, stimulus.delta, task.thistrial.interval));
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
    % update appropriate staircase
    staircaseNum = find(stimulus.pedestal==task.thistrial.pedestal);
  
    stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
    %  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
    disp(sprintf('  No response. Delta: %0.2f',task.thistrial.delta));
  end
end


if task.thistrial.thisseg == 8 %Adaptation topup
    mglPlaySound(stimulus.adaptNum);
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



%3rd segment: blank period

if task.thistrial.thisseg == 6 %response cue... where was the stimulus?

  feedback = [1 1 1]; %all other things equal... fix cross is white
  
    % set colors only after the subject has responded
  if ~isnan(task.thistrial.correct)
      if task.thistrial.correct
          feedback = [0 1 0];
      else
          feedback = [1 0 0];
      end
  end
  
  mglFixationCross(0.75, 0.75, feedback)

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
  staircaseNum = find(stimulus.pedestal==task.thistrial.pedestal);
  stimulus.s(staircaseNum) = doStaircase('update',stimulus.s(staircaseNum),task.thistrial.correct);
%  threshold = doStaircase('threshold',stimulus.s(staircaseNum));
  disp(sprintf('  Response: %i Correct: %i Delta: %0.2f',task.thistrial.responseInterval,task.thistrial.correct,task.thistrial.delta));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task,stimfile)

% set the maximum amplitude
stimulus.maxStim = 1;

%Basic stimulus waveform - Sine Wave, Amplitude = 1 (or = 2 peak to peak),
% Offset from beginning of file by 0.01s
onTime = task{2}.onTime;
freq = task{2}.freq;
sampleRate = 8192;

stimBase = sin(freq * (2*pi*(1:round(onTime*sampleRate))/sampleRate));

% waitTime = 0.2; %Wait time now set as length of 3rd segment - no
% associated audio waveform

stimulus.stimBase = stimBase;
% stimulus.waitBase = waitBase;

stimulus.deviceID = 1; %deviceID for sound


%Adaptor topup
onTime = 5;
freq = task{2}.freq;
sampleRate = 8192;

adaptLevel = task{2}.parameter.adaptLevel;
adaptBase = adaptLevel*sin(freq * (2*pi*(1:round(onTime*sampleRate))/sampleRate));
stimulus.adaptNum = mglInstallSound(adaptBase);
mglSetSound(stimulus.adaptNum, 'deviceID', stimulus.deviceID);

% initialize the staircases


% determine pedestals
params = getTaskParameters(myscreen,task{2});
pedestal = params.originalTaskParameter.pedestal;
numPedestal = length(pedestal);
stimulus.pedestal = pedestal;
stimulus.numPedestal = numPedestal;

stimfile = [];
% first time, initialize the staircases,
for iStaircase = 1:numPedestal
    stimulus.s(iStaircase) = doStaircase('init','upDown','nup=1','ndown=2','initialThreshold=0.3','initialStepsize=0.1','nTrials=30','stepRule=levitt','minThreshold=0');
end
