%% Figure 1 - optimal leaving times vs participant behaviour

RR_Leave = [21.8678 18.5632]; % this is the average background RR at which subjects should leave, for Rich and Poor env... - 
... respectively.
    
clear optLT
A=[32.5 45 57.5];a=0.075;
for e=1:2 %each env
    for p=1:3 % each patch type
        optLT(e,p)=(log(RR_Leave(e)/A(p)))/-a;
    end
end