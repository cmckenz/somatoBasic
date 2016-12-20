function [] = playUnityStimulus(stimTime, stimFreq, numChannels, deviceID);
% playUnityStimulus(stimTime, stimFreq, numChannels, deviceID);
%
% This function delivers a somatosensory (or auditory) stimulus of length 
% "stimulus" with the peak value of 1, or 2 peak to peak. It may
% be helpful for checking device calibration with the oscillo.
%
% author: cam mckenzie, august 2016

myscreen = initScreen('test');

sampleRate = 8192;

baseArray = 2*pi*(1:(round(stimTime)*sampleRate))/sampleRate;
toneArray = 0.5*sin(stimFreq*baseArray);

[a b] = size(toneArray);

while a < numChannels
    toneArray = [toneArray; toneArray(1,:)];
    [a b] = size(toneArray);
end

toneNum = mglInstallSound(toneArray);
mglSetSound(toneNum, 'deviceID', deviceID)
mglPlaySound(toneNum);

mglClose;