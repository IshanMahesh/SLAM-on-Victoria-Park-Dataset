function initialize_new_landmark(z, R)

global Param;
global State;

Landmark = zeros(3,1);
State.Ekf.Observed_landmarks = [State.Ekf.Observed_landmarks,z(3,1)];
Landmark(1,1)= State.Ekf.predMu(1) + z(1,1)*cos(z(2,1)+State.Ekf.predMu(3));  %%APPLY OBSERVATION MODEL
Landmark(2,1)= State.Ekf.predMu(2) + z(1,1)*sin(z(2,1)+State.Ekf.predMu(3));
Landmark(3,1)= z(3,1);
State.Ekf.predMu = [State.Ekf.predMu; Landmark]; %%AUGMENT MU

if Param.sim==1
    A=[12345678923456789, 0; 0,12345678923456789];
else
    A=[10,0 0;0,10 0;0 0 1];
end
State.Ekf.predSigma = [State.Ekf.predSigma, zeros(size(State.Ekf.predSigma,1),3); zeros(3,size(State.Ekf.predSigma,2)), A]; %%AUGMENT SIGMA