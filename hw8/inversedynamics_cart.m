function [Fx,Fy]=inversedynamics_cart(theta,ddsX,ddsY,domega,par)
%This function calculates forces Fx and Fy as functions of given motion
%(not checking for violation of constraints!)
%par: parameter struct containing length and mass properties of the cart
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

%----------------------------
%Lagrange:
%----------------------------
A=[m 0 -2*Is/l*sin(theta);
    0 m 2*Is/l*cos(theta);
    sin(theta) -cos(theta) l/2];

%b=[ Fx*cos(theta) - 2*Fy*sin(theta);
%    Fx*sin(theta)+2*Fy*cos(theta);
%    -dsX*omega*cos(theta)-dsY*omega*sin(theta)];

%re-sort the equations to solve for Fx and Fy:
ddq=[ddsX;ddsY;domega];

AF=[cos(theta) -2*sin(theta);
    sin(theta) 2*cos(theta)];

bF=A(1:2,:)*ddq;

F=AF\bF;
Fx=F(1);
Fy=F(2);



