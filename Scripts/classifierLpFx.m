function [aa bb cc, UnAidedPerf] = classifierLpFx(simStructNew,simStructTDOA,simStructBaseline,...
    TimeThresh, SimThresh, betaParm1, beta2)


simStructNew.betaParm1= betaParm1;
simStructNew.betaParm2 = beta2;

truthTable=createSpeciesPreds(simStructNew);


simStructNew.betaParm1 = betaParm1;
simStructNew.betaParm2 = beta2;
simStructNew.truthTable =truthTable;

simStructTDOA.betaParm1 = betaParm1;
simStructTDOA.betaParm2 = beta2;
simStructTDOA.truthTable =truthTable;

simStructBaseline.betaParm1 = betaParm1;
simStructBaseline.betaParm2 = beta2;
simStructBaseline.truthTable =truthTable;

[aa, UnAidedPerf2]=runClassPerfLoop(simStructBaseline,...
    TimeThresh);
[cc, UnAidedPerf1]= runClassPerfLoop(simStructNew, TimeThresh, SimThresh);
[bb, UnAidedPerf]= runClassPerfLoop(simStructTDOA, TimeThresh, SimThresh);


end