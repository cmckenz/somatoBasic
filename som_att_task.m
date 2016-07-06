% somdis.m
%
%        $Id:$ 
%      usage: som_att_task()
%         by: justin gardner/cameron mckenzie
%       date: 05/07/16
%    purpose: Simple somato strength localizer - no responses, attentional
%              visual task added, no pedestals
%
function myscreen = som_att_task()

% check arguments
if ~any(nargin == [0])
  help taskTemplate
  return
end

% initalize the screen

myscreen = initScreen;

% task parameters for somato stimulation
task{1}{1}.waitForBacktick = 1;
task{1}{1}.segmin = [0.25 0.75 1]; %added a 250ms initial segment to set up stim sound
task{1}{1}.segmax = [0.25 0.75 9]; %set stim/play stim/response delay/response interval/inter-trial interval
task{1}{1}.synchToVol = [1 0 0]; %synch before stimulus, responses
task{1}{1}.getResponse = [0 0 0];
task{1}{1}.parameter.pedestal = [0];
task{1}{1}.parameter.stimulus = [1 2 3];
task{1}{1}.randVars.uniform.side = [-1 1];
task{1}{1}.random = 1;

testing = true;

if testing
  task{1}{1}.waitForBacktick = 0;
  task{1}{1}.segmin = [0.25 0.75 1]; %added a 250ms initial segment to set up stim sound
  task{1}{1}.segmax = [0.25 0.75 3]; %set stim/play stim/response delay/response interval/inter-trial interval
  task{1}{1}.parameter.pedestal = [0];
  task{1}{1}.synchToVol = [0 0 0];
end

%task parameters for attentional task
task{2}{1}.segmin = 0.8
task{2}{1}.segmax = 1.07%display letters for a not exactly half a second (which is TR - avoid simple divisibility by TR)
task{2}{1}.parameter.letterNum = [1:26];
task{2}{1}.random = 1;
task{2}{1}.synchToVol = 0;
task{2}{1}.randVars.uniform.lure = [zeros(1,9) 1];



% initialize the task
for phaseNum = 1:length(task{1})
    [task{1}{phaseNum} myscreen] = initTask(task{1}{phaseNum}, myscreen, @startSegmentCallback1, @screenUpdateCallback1, @responseCallback1);
    [task{2}{phaseNum} myscreen] = initTask(task{2}{phaseNum}, myscreen, @startSegmentCallback2, @screenUpdateCallback2, @responseCallback2);
end

% init the stimulus
global stimulus;
myscreen = initStimulus('stimulus',myscreen);
stimulus = myInitStimulus(stimulus,myscreen,task);
%myscreen.feedback = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run the eye calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%myscreen = eyeCalibDisp(myscreen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main display loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phaseNum = 1;
while (phaseNum <= length(task{1})) && ~myscreen.userHitEsc
  % update the task
  [task{1} myscreen phaseNum] = updateTask(task{1},myscreen,phaseNum);
  [task{2} myscreen phaseNum] = updateTask(task{2},myscreen,phaseNum);
  % flip screen
  myscreen = tickScreen(myscreen,task);
end

% if we got here, we are at the end of the experiment
myscreen = endTask(myscreen,task);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called at the start of each segment - somato
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = startSegmentCallback1(task, myscreen)

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
% function that gets called at the start of each segment - somato
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [task myscreen] = startSegmentCallback2(task, myscreen)

%not much to do here really...
if task.thistrial.lure == 1
    display('Target trial')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback1(task, myscreen)

%nothing to do here really...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that gets called to draw the stimulus each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = screenUpdateCallback2(task, myscreen)

mglClearScreen();

if task.thistrial.lure == 0
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    lNum = task.thistrial.letterNum;
    mglTextDraw(letters(lNum), [0 0], 0, 0)
end

if task.thistrial.lure == 1
    mglTextDraw( '8', [0 0], 0, 0)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback1(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    responseCallback    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [task myscreen] = responseCallback2(task,myscreen)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to init the stimulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimulus = myInitStimulus(stimulus,myscreen,task)

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
pedestal = params{1}.originalTaskParameter.pedestal;
numPedestal = length(pedestal);
stimulus.pedestal = pedestal;
stimulus.numPedestal = numPedestal;
stimulus.stimBase = stimBase;
stimulus.deviceID = 2;

threshold1 = 0.05;
threshold2 = 0.05;
threshold3 = 0.05;
threshold4 = 0.05;

stimulus.stimVal(1,1:3) = threshold1*[1 2 4];
stimulus.stimVal(2,1:3) = threshold2*[1 2 4];
stimulus.stimVal(3,1:3) = threshold3*[1 2 4];
stimulus.stimVal(4,1:3) = threshold4*[1 2 4];
  
  


