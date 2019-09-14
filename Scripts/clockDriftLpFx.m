function [aa bb cc, UnAidedPerf] = clockDriftLpFx(simStructNew,simStructTDOA,simStructBaseline,...
    TimeThresh, SimThresh)



[aa, ~]=runClassPerfLoop(simStructBaseline,TimeThresh);
[cc, ~]= runClassPerfLoop(simStructNew, TimeThresh, SimThresh);
[bb, UnAidedPerf]= runClassPerfLoop(simStructTDOA, TimeThresh, SimThresh);


end