function doGraphics_Sarvesh(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;
plot(State.Ekf.aug_mu(1,:), State.Ekf.aug_mu(2,:), '--b');
hold on
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on;

robot_graph =plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3);

BB=100;
axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);
textLabels = false; 
N=length(State.Ekf.Observed_landmarks);
if N > 0
    for i = 1:N
        plotcov2d(State.Ekf.mu(3+2*i-1,1), State.Ekf.mu(3+2*i,1), State.Ekf.Sigma(3+2*i-1:3+2*i,3+2*i-1:3+2*i), 'g', 0, 'y', 0.5, 3);
    end
end
title('Victoria Park')
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos(b+tr-pi/2);
    yl = yr + r*sin(b+tr-pi/2);
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end

hold off;