function dx=cart_equationsofmotion(t,x,par)
%function dx=cart_equationsofmotion(t,x,par)
%This function contains the 2D equations of motion of a three-wheeled cart
%driving around on the ground, as used in homework
%assignments of the course "Advanced Dynamics" at TUD.
%
%Inputs: 
%t (time) 
%state vector x: Cartesian x/y positions and velocities of cart center 
%               of mass & heading angle theta
%par: parameter struct containing length and mass properties of the cart
%
%Output: dx: derivatives of the states
%
%Author: H. Vallery, October 2014

%----------------------------
%extract parameters:
%----------------------------

%mass parameters of the cart:
m=par.m;% mass of cart
Is=par.Is;% moment of inertia about center of mass

%geometry of the cart:
l=par.length_cart;%[m], length of the cart
%d=par.width_cart;%[m] width of the cart, not needed


%----------------------------
%read current states:
%----------------------------

sX=x(1);
sY=x(2);
theta=x(3);

dsX=x(4);
dsY=x(5);
omega=x(6);

Ts = 0.01;
n = min(round(t/Ts) + 1, par.N);
Fx = par.F(n, 1);
Fy = par.F(n, 2);

%----------------------------
%Lagrange:
%----------------------------

A=[m 0 -2*Is/l*sin(theta);
    0 m 2*Is/l*cos(theta);
    sin(theta) -cos(theta) l/2];

b=[ Fx*cos(theta) - 2*Fy*sin(theta);
    Fx*sin(theta)+2*Fy*cos(theta);
    -dsX*omega*cos(theta)-dsY*omega*sin(theta)];

ddq=A\b;
ddsX=ddq(1);
ddsY=ddq(2);
domega=ddq(3);

%----------------------------
%calculate state derivatives:
%----------------------------

dx=[dsX;dsY;omega;ddsX;ddsY;domega];

