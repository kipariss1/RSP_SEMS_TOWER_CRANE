clearvars

%%%parametry
L = 0.2; %m
g = 9.810; %ms^-2

%%%stavovy popis - symbolicky
%    beta       alpha    x_w        theta     beta_dot  
syms x1         x2       x3         x4        x5          real

%    alpha_dot  x_w_dot  theta_dot  x_w_ddot  theta_ddot
syms x6         x7       x8         u1        u2          real


%vyjadreni popisu ve stavove formulaci
dx1 = x5;
dx2 = x6;
dx3 = x7;
dx4 = x8;
dx5 = -1/(2*L)*(2*g*cos(x2)*sin(x1) + 4*x7*x8*cos(x1) - L*x8^2*sin(2*x1)*(cos(x2))^2 - 2*x3*x8^2*sin(x1)*sin(x2) ...
- 4*L*x6*x8*cos(x2)*(cos(x1))^2 + L*x6^2*sin(2*x1) + 2*sin(x1)*sin(x2)*u1 + 2*x3*cos(x1)*u2 - 2*L*sin(x2)*u2);
dx6 = -1/(L*cos(x1))*(g*sin(x2) + x3*x8^2*cos(x2) - L*x8^2*sin(x2)*cos(x1)*cos(x2) + 2*L*x5*x8*cos(x1)*cos(x2) ...
- 2*L*x5*x6*sin(x1) - cos(x2)*u1 + L*sin(x1)*cos(x2)*u2);
dx7 = u1;
dx8 = u2;

%vektor stavovych promennych
x = [x1 x2 x3 x4 x5 x6 x7 x8]';
dx = [dx1 dx2 dx3 dx5 dx6 dx7 dx8]';

u = [u1 u2]';

%definice pracovniho bodu
x1s = 0;
x2s = 0;
x3s = 0.3; %max 0.6
x4s = 0;
x5s = 0;
x6s = 0;
x7s = 0; %0.02;
x8s = 0; %12*3.1416/180;
xs = [x1s x2s x3s x4s x5s x6s x7s x8s]';

u1s = 0;
u2s = 0;
us = [u1s,u2s]';

%linearizace
f1 = dx1;
f2 = dx2;
f3 = dx3;
f4 = dx4;
f5 = dx5;
f6 = dx6;
f7 = dx7;
f8 = dx8;

A = jacobian([f1;f2;f3;f4;f5;f6;f7;f8],x);
A_num = subs(A,x,xs);
A_num = subs(A_num,u,us);
A_num = double(A_num);

B = jacobian([f1;f2;f3;f4;f5;f6;f7;f8],u);
B_num = subs(B,x,xs);
B_num = double(B_num);

%matice C a D
C = [1 0 0 0 0 0 0 0
     0 1 0 0 0 0 0 0
     0 0 1 0 0 0 0 0
     0 0 0 1 0 0 0 0];

D = zeros(4, 2);

%%%Stavovy popis - souhrn matic A, B, C, D
crane_ss = ss(A_num, B_num, C, D);

%TF = tf(A_num, B_num,C,D)
TF = tf(crane_ss);
figure(1);
pzmap(crane_ss);

figure(2);
step(crane_ss)
% odezva na vstup 1 (x_w) u vystupu 2 (alpha) je kmitava (coz je logicky)
% odezva na vstup 2 (theta) u vystupu 1 (beta) je taky kmitava
% (coz je taky logicky) (podelny uhel reaguje na posuv voziku, ...
% pricny uhel reaguje na otaceni jerabu)
% !soustava je kmitava, na mezi stability coz odpovida mape nul a polu!

figure(3);
nyquist(crane_ss);

EigCrane = eig(A_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State space --> Transfer function model%
%           Decoupling                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

crane_tf = tf(crane_ss);

%transfer function G from 1. input to 2. output 
TFI1toQ2 = tf(crane_tf.Numerator(2,1), crane_tf.Denominator(2,1));
%transfer function G from 2. input to 1. output
TFI2toQ1 = tf(crane_tf.Numerator(1,2), crane_tf.Denominator(1,2));


%Vysetrovani stability

[Gm12, Pm12] = margin(TFI1toQ2);
[Gm21, Pm21] = margin(TFI2toQ1);
figure(4);
margin(TFI1toQ2);
figure(5);
margin(TFI2toQ1);

%%PID%%
rp = 10;
ri = 49.05;
rd = 0;

reg_tf = tf([rp ri rd],[1 0]);

%%Open-loop systems%%

crane_tf_ol_12 = TFI1toQ2*reg_tf;
crane_tf_ol_21 = TFI2toQ1*reg_tf;

%%Closed-loop systems%%

crane_tf_cl_12 = feedback(crane_tf_ol_12, 1);
crane_tf_cl_21 = feedback(crane_tf_ol_21, 1);
[Gm12c, Pm12c] = margin(crane_tf_cl_12);
[Gm21c, Pm21c] = margin(crane_tf_cl_21);
figure(6);
margin(crane_tf_cl_12);
figure(7);
margin(crane_tf_cl_21);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%END%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



