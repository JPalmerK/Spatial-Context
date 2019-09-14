function [aa bb cc] = classifierLpFx(simStructNew,simStructTDOA,simStructBaseline,...
    TimeThresh, SimThresh, betaParm1, beta2)




simStructNew.betaParm1 = betaParm1;
simStructNew.betaParm2 = beta2;

simStructTDOA.betaParm1 = betaParm1;
simStructTDOA.betaParm2 = beta2;

simStructBaseline.betaParm1 = betaParm1;
simStructBaseline.betaParm2 = beta2;

[ExpScoresMeth, UnAidedPerf]= runClassPerfLoop(simStructTDOA, TimeThresh, SimThresh);
[ExpScoresMeth, UnAidedPerf]= runClassPerfLoop(simStructNew, TimeThresh, SimThresh);
[ExpScoresMeth, UnAidedPerf]=runClassPerfLoop(simStructBaseline,...
    TimeThresh);


end