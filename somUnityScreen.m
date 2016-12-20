function [] = somUnityScreen(onTime, deviceID, screenName)

%Just to test the sound... set a time for test wave and a deviceID. If
%deviceID is set to -1 should show available devices.

myscreen = initScreen(screenName);


freq = 80
sampleRate = 8192;

stimBase = sin(freq * (2*pi*(1:round(onTime*sampleRate))/sampleRate));

stimBase = [stimBase; stimBase];

soundNum = mglInstallSound(stimBase);

if deviceID == -1
    mglSetSound(soundNum, 'deviceID')
    else
    mglSetSound(soundNum, 'deviceID', deviceID)
    mglPlaySound(deviceID)
    pause(onTime)
end

mglClose