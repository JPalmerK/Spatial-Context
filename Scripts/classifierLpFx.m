function [aa bb cc] = classifierLpFx(simStructNew,simStructTDOA,simStructBaseline,...
    TimeThresh, SimThresh, betaParm1, beta2)

            simStructNew.betaParm1 = betaParm1;
            simStructNew.betaParm2 = beta2;
            
            simStructTDOA.betaParm1 = betaParm1;
            simStructTDOA.betaParm2 = beta2;
            
            simStructBaseline.betaParm1 = betaParm1;
            simStructBaseline.betaParm2 = beta2;
            
            bb= runClassPerfLoop(simStructTDOA, TimeThresh, SimThresh);
            cc= runClassPerfLoop(simStructNew, TimeThresh, SimThresh);
            aa =runClassPerfLoop(simStructBaseline,...
                TimeThresh);
end