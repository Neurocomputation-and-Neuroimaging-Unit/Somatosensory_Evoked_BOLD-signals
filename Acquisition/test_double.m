%% try double pulse
ioObj = io64;
status = io64(ioObj);
address = hex2dec('3FF8');
data_out = 128;
data_last = 0;

catch_ISI = 0.05; % in s, that would be 50 ms

io64(ioObj,address,data_out); % send
pause(0.001)
io64(ioObj,address,data_last); % send
pause(catch_ISI)
io64(ioObj,address,data_out); % send again if catch trial
pause(0.001)
io64(ioObj,address,data_last); %3 send