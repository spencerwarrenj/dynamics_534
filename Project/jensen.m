function [wx,wy,wz,psi,theta,phi]=jensen(wx0,wy0,wz0,psi0,theta0,phi0,t,Mx,My,Mz);

% OG Author:   Mark Colton
% Edited: Spencer Jensen
% Date:     3/29/21

% Initial conditions
y0 = [psi0,theta0,phi0,wx0,wy0,wz0];
%change degrees to rad
y0 = y0*pi/180;

t_vec = t;

% Call 'ode45' to solve ODE.  'ode45' calls the function
% 'csm' repeatedly, which returns the two 
% derivatives at each time step.  'ode45' uses the returned
% derivatives to calculate the solution at each time step.
% The solution is returned in the variable y.
[t,y] = ode45(@csm_function,t, y0);

% Extract solution from y
psi = y(:,1)*180/pi;
theta = y(:,2)*180/pi;
phi = y(:,3)*180/pi;
wx = y(:,4)*180/pi;
wy = y(:,5)*180/pi;
wz = y(:,6)*180/pi;

    function ydot = csm_function(t, y)
        % This m-file contains the state derivatives function
        % for the solution of a system of equations of motion
        % using the matrix approach. 
        
        % OG_author: Mark Colton
        % Edited: Spencer Jensen
        % 3/29/20

        %general case
        I = [40823.1019533526,1537.80758482821,-3179.29746413299;1537.80758482821,90593.5331952687,128.577121406529;-3179.29746413299,128.577121406529,98742.9021241113];
        %simplified case
%         I = [40482.0123511000,0,0;0,90358.3606572137,0;0,0,98636.9850599137]
        
        % Calculate the M matrix for the current time step (note
        % that M is a function of the states, specifically x2,
        % which is phi)
        M = [eye(3),zeros(3);zeros(3),I];

        % Calculate the F vector for the current time step
        t_small = find(t>=t_vec,1)
        t_big = find(t<t_vec,1)
        mx = (Mx(t_small)*(t_vec(t_big)-t)+Mx(t_big)*(t-t_vec(t_small)))/(t_vec(t_big)-t_vec(t_small))
        my = (My(t_small)*(t_vec(t_big)-t)+My(t_big)*(t-t_vec(t_small)))/(t_vec(t_big)-t_vec(t_small))
        mz = (Mz(t_small)*(t_vec(t_big)-t)+Mz(t_big)*(t-t_vec(t_small)))/(t_vec(t_big)-t_vec(t_small))
        
        mom = [mx;my;mz];
        jangle = [y(1);y(2);y(3)];
        omega = [y(4);y(5);y(6)];
        
        F = [1/cos(jangle(2))*(omega(2)*sin(jangle(3))+omega(3)*cos(jangle(3))); ...
            omega(2)*cos(jangle(3))-omega(3)*sin(jangle(3)); ...
            1/cos(jangle(2))*(omega(2)*sin(jangle(2))*sin(jangle(3))+omega(3)*sin(jangle(2))*cos(jangle(3)))+omega(1); ...
            mom-cross(omega,I*omega)];

        % Calculate the state derivatives using the matrix approach.
        % Note that 'M\F' is another way of calculating 'inv(M)*F',
        % but allows MATLAB to select the most numerically efficient
        % and accurate method for calculating the inverse.
        ydot = M\F;
    end
end
