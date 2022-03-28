%% Hw4-2
clear;clc;close;warning off;

% Read file
filename = "Hw4-2.xls"; sheet = "Sheet1"; 
sensor_data = xlsread(filename, sheet);
n = length(sensor_data);

% Initial setting
g = 9.80665; u = 1;
q = zeros(4,n+1); q(:,1) = [1;0;0;0];
ss = [sensor_data']; se = zeros(3,n);
atd = zeros(3,n+1); atd(:,1) = [0;0;0];

ang = deg2rad(0);
q_g = [cos(ang/2);0;0;sin(ang/2)];

% Cost fun & Jacobian & gradient
e = @(q1,q2,q3,q4,ax,ay,az) [-2*g*(q2*q4-q1*q3) - ax;
                             -2*g*(q1*q2+q3*q4) - ay;
                             -2*g*(0.5-(q2^2)-(q3^2)) - az];
                         
Jc = @(q1,q2,q3,q4) (-2*g).*[-q3 q4 -q1 q2;
                              q2 q1 q4 q3;
                              0 -2*q2 -2*q3 0];
                          
C = @(q1,q2,q3,q4,ax,ay,az,g) (ax-2*g*(q1*q3 - q2*q4))*(ax - 2*g*(q1*q3 - q2*q4)) ...
                            + (ay + 2*g*(q1*q2 + q3*q4))*(ay + 2*conj(g)*(q1*q2 + q3*q4)) ...
                            + (az - 2*g*(q2^2 + q3^2 - 1/2))*(az - 2*g*(q2^2 + q3^2 - 1/2));
 
grd_c = @(q1,q2,q3,q4,ax,ay,az,g) [2*g*q2*(ay + 2*g*(q1*q2 + q3*q4)) - 2*g*q3*(ax - 2*g*(q1*q3 - q2*q4));
                                   2*g*q4*(ax - 2*g*(q1*q3 - q2*q4)) + 2*g*q1*(ay + 2*g*(q1*q2 + q3*q4)) - 4*g*q2*(az - 2*g*(q2^2 + q3^2 - 1/2));
                                   2*g*q4*(ay + 2*g*(q1*q2 + q3*q4)) - 2*g*q1*(ax - 2*g*(q1*q3 - q2*q4)) - 4*g*q3*(az - 2*g*(q2^2 + q3^2 - 1/2));
                                   2*g*q2*(ax - 2*g*(q1*q3 - q2*q4)) + 2*g*q3*(ay + 2*g*(q1*q2 + q3*q4))];
% Estimate attitude
for i = 1:n
    
    q(:,i+1) = gredient_decsent(q_g, ss(:,i), g, u, C, grd_c);
    Q = quaternion(q(1,i+1), q(2,i+1), q(3,i+1), q(4,i+1));
    atd(:,i+1) = quat2eul(Q,"XYZ")'.*(180/pi);
    G = [0,0,-g];
    ptrot = rotateframe(Q, G);    
    se(:,i) = ptrot';
end

% Least square error
E = 0;

for o = 1:n
    
    error = ss(:,o) - se(:,o);
    E = E + norm(error);
end
fprintf("Least square error: %.4f \n",E);

% Plot attitude trajectory
tl = zeros(n+1,3); ro = q';

for j = 1:n+1
    
    ax = plotTransforms(tl(j,:),ro(j,:));
    pause(1/30);
end

% Gradient decsent
function [q] = gredient_decsent(q, s, g, u, C, grd_c)
    
    Grd_c = grd_c(q(1,1), q(2,1), q(3,1), q(4,1), s(1,1), s(2,1), s(3,1), g);
    while norm(Grd_c) > 0.005
        
        Grd_c = grd_c(q(1,1), q(2,1), q(3,1), q(4,1), s(1,1), s(2,1), s(3,1), g);
        u = BTLS(q, s, g, u, C, Grd_c);
        q = q - u*(Grd_c/norm(Grd_c));
    end
end

function [u] = BTLS(q, s, g, u, C, grd_c)

    alpha = 0.1;
    beta =  0.707;

    delta_x = -(grd_c/norm(grd_c));
    dx_g = (grd_c')*delta_x;
    
    Co = C(q(1,1), q(2,1), q(3,1), q(4,1), s(1,1), s(2,1), s(3,1), g);
    Cn = C(q(1,1)+u*delta_x(1,1), q(2,1)+u*delta_x(2,1), q(3,1)+u*delta_x(3,1), ...
           q(4,1)+u*delta_x(4,1), s(1,1), s(2,1), s(3,1), g);
    
    while Cn > (Co + alpha*u*dx_g)
        
        u = u*beta;
        Cn = C(q(1,1)+u*delta_x(1,1), q(2,1)+u*delta_x(2,1), q(3,1)+u*delta_x(3,1), ...
               q(4,1)+u*delta_x(4,1), s(1,1), s(2,1), s(3,1), g);
    end
end