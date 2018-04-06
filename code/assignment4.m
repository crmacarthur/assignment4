%% Question 1
% $V_i = V_{in}$
%
% $G_1(V_2-V_1)+C(\frac{d(V_2-V_1)}{dt} )+G_2V_2-I_L=0$
%
% $$V_2-V_3-L \frac {dI_l}{dt}=0$$
%
% $$-I_L+G_3V_3=0$$
%
% $$V_4-\alpha I_3=0$$
%
% $$G_3V_3-I_3=0$$
%
% $$G_4(V_O-V_4)+G_O V_O =0$$
%
% in frequency domain
%
% $V_i = V_{in}$
%
% $G_1(V_2-V_1)+C(j \omega (V_2-V_1) )+G_2V_2-I_L=0$
%
% $$V_2-V_3-L (j \omega)I_l=0$$
%
% $$-I_L+G_3V_3=0$$
%
% $$V_4-\alpha I_3=0$$
%
% $$G_3V_3-I_3=0$$
%
% $$G_4(V_O-V_4)+G_O V_O =0$$





%variables [v1 v2 v3 v4 v5 Iin IL]
close all
clear

g1=1;
g2=0.5;
g3=.1;
g4=10;
g5=1e-3;
a=100;
c=0.25;
l=0.2;






% G matrix
G=[1,0,0,0,0,0,0;
   -g2,g1+g2,-1,0,0,0,0;
   0,1,0,-1,0,0,0;
   0,0,-1,g3,0,0,0;
   0,0,0,0,-a,1,0;
   0,0,0,g3,-1,0,0;
   0,0,0,0,0,-g4,g4+g5]

% C matrix
C=[0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0;
   0,0,-l,0,0,0,0;
   0,0,0,0,0,0,0;
   0,0,0,0,0,0,0;
   0,0,0,0,0,0,0;
   0,0,0,0,0,0,0]

q=[];
w=[];

% DC Sweep from -10 to 10
for vin=-10:10
    %F Vector
    F=[vin;0;0;0;0;0;0];
    V=G\F;
    q=[q,V(4)];
    w=[w,V(7)];
end
figure(1)
plot (-10:10,q)
title('Voltage at V3 vs Vin for Vin =-10 to 10 Vdc')
xlabel('Vin (V)')
ylabel('Voltage at V3 (V)')


figure(2)
plot (-10:10,w)
title('Vout vs. Vin for Vin =-10 to 10 Vdc')
xlabel('Vin (V)')
ylabel('Vout (V)')

q=[];
w=[];
vin=1;
for z=1:1000
    F=[vin;0;0;0;0;0;0];
    V=(G+z*1j*C)\F;
    w=[w,V(7)];
end

figure(3)
semilogx(1:1000,abs(w))
title('V0 vs. frequency with 1V amplitude')
xlabel('Frequency (rad/s)')
ylabel('Vout (V)')


figure(4)
semilogx(1:1000,20*log10(abs(w)))
title('V0 vs. frequency with 1V amplitude')
xlabel('Frequency (rad/s)')
ylabel('Vout (dB)')


w=[];
for z=1:1000
    c=0.25+0.05*randn();
    F=[vin;0;0;0;0;0;0];
    C= [0,0,0,0,0,0,0;
       -c,c,0,0,0,0,0;
       0,0,-l,0,0,0,0;
       0,0,0,0,0,0,0;
       0,0,0,0,0,0,0;
       0,0,0,0,0,0,0;
       0,0,0,0,0,0,0];
   V=(G+pi*1j*C)\F;
   w=[w,V(7)];
end


figure(5)
histogram(20*log10(sqrt(real(w).^2+imag(w).^2)))
title('Histogram of gain')
xlabel('Values of gain (dB)')
ylabel('count')

%% question 2

% This is an amplifier circuit.  it has a bandpass response due to the
% capcitor and inductor.  The resistor parallel to the capacitor limits the
% low frequency filtering, so there will be less than a 1st order drop off
% for low frequency components.  The high frequency components will
% experience first order drop off due to the inductor


% Create the voltage step (F1)
t = linspace(0,1,1000);
vin=1;
F1=0;
for i=1:31
    F1(i,1:7)=[0,0,0,0,0,0,0];
end
for i=32:1000
    F1(i,1:7)=[vin;0;0;0;0;0;0];
end


v1=0;
dt=0.001;
% voltage step output
v1(1:7,1)=(C/dt+G)^-1 *(F1(1,:)');
for i=2:1000
    v1(:,i)=(C/dt+G)^-1 *(C*v1(:,i-1)/dt+F1(i,:)');   
end

figure(6)

plot(t,v1(7,:))
hold on 
plot(t,F1(:,1))
legend('Vout','input step')
title('DC step responce')
xlabel('Time (s)')
ylabel('Voltage (V)')

% create the sin input finction
v2=0;
t = linspace(0,1,1000);
f=1/.03;
F2=0;
for i=1:1000
    F2(i,1:7)=[sin(2*pi*f*t(i)),0,0,0,0,0,0];
end

% sin output
v2(1:7,1)=(C/dt+G)^-1 *(F2(1,:)');
for i=2:1000
    v2(:,i)=(C/dt+G)^-1 *(C*v2(:,i-1)/dt+F2(i,:)');   
end

figure(7)
plot(t,v2(7,:))
hold on 
plot(t,F2(:,1))
legend('Vout','input signal')
title('DC sin responce')
xlabel('Time (s)')
ylabel('Voltage (V)')


%guassian function

v3=0;
t = linspace(0,1,1000);
f=1/.03;
F3=0;
for i=1:1000
    F3(i,1:7)=[exp(-1/2*((t(i)-0.06)/0.03)^2),0,0,0,0,0,0];
end

%guassian output
v3(1:7,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:1000
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end



figure(8)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
legend('Vout','input signal')
title('DC Guassian responce')
xlabel('Time (s)')
ylabel('Voltage (V)')

%guassian output larger timestep
dt=0.001;
t = linspace(0,1,100);
v4(1:7,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:100
    v4(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(9)
plot(t,v4(7,:))
hold on 
t = linspace(0,1,1000);
plot(t,F3(:,1))
legend('Vout','input signal')
title('DC Guassian responce')
xlabel('Time (s)')
ylabel('Voltage (V)')

figure(10)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(v1(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F1(:,1)))))
legend('Vout','input signal')
title('DC Step Fourier Transform')
xlabel('Frequency(rad/s)')
ylabel('Voltage (V)')

figure(11)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(v2(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F2(:,1)))))
legend('Vout','input signal')
title('Sin Fourier Transform')
xlabel('Frequency (rad/s)')
ylabel('Voltage (V)')

figure(12)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(v3(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F3(:,1)))))
title('Guassian Fourier Transform')
legend('Vout','input signal')
xlabel('Frequency (rad/s)')
ylabel('Voltage (V)')

% As we can see in figure 9, when the time step is enlarged the simulation
% behaves unexpectedly.  The guassian peak becaomes very delayed in time on
% the output, and much wider than expected.

%% Question 3

iin=0.001*randn();
c2=.00001;

% G matrix
G=[1,0,0,0,0,0,0,0;
   -g2,g1+g2,-1,0,0,0,0,0;
   0,1,0,-1,0,0,0,0;
   0,0,-1,g3,0,0,0,-1;
   0,0,0,0,-a,1,0,0;
   0,0,0,g3,-1,0,0,0;
   0,0,0,0,0,-g4,g4+g5,0;
   0,0,0,0,0,0,0,1];

% C matrix
C=[0,0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0,0;
   0,0,-l,0,0,0,0,0;
   0,0,0,c2,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0];

F=[vin;0;0;0;0;0;0;iin];

v3=0;
t = linspace(0,1,1000);
dt=0.001;
f=1/.03;
F3=0;

for i=1:1000
    iin=0.001*randn();
    F3(i,1:8)=[exp(-1/2*((t(i)-0.06)/0.03)^2),0,0,0,0,0,0,iin];
end

v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:1000
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(13)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with noise')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')

figure(14)
semilogy(linspace(-500,500,1000),fftshift(abs(fft(v3(7,:)))))
hold on 
semilogy(linspace(-500,500,1000),fftshift(abs(fft(F3(:,1)))))
title('Guassian Fourier Transform with Noise')
legend('Vout','input signal')
xlabel('Frequency (rad/s)')
ylabel('Voltage (V)')

% change the c value
c2=.00002;
C=[0,0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0,0;
   0,0,-l,0,0,0,0,0;
   0,0,0,c2,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0];

v3=0;
v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:1000
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(15)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with cn=0.00002')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')

c2=.0002;
C=[0,0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0,0;
   0,0,-l,0,0,0,0,0;
   0,0,0,c2,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0];

v3=0;
v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:1000
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(16)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with cn=0.0002')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')

c2=.002;
C=[0,0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0,0;
   0,0,-l,0,0,0,0,0;
   0,0,0,c2,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0];

v3=0;
v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:1000
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(17)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with cn=0.002')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')

% Decreasing the capacitor increases the noise in the output.  As the
% capcacitor value is increased to 0.002, the noise decreases until it can
% not be seen.

c2=.00001;
C=[0,0,0,0,0,0,0,0;
   -c,c,0,0,0,0,0,0;
   0,0,-l,0,0,0,0,0;
   0,0,0,c2,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0;
   0,0,0,0,0,0,0,0];

% different time step

v3=0;
t = linspace(0,1,750);
dt=0.0013;
f=1/.03;
F3=0;
for i=1:750
    iin=0.001*randn();
    F3(i,1:8)=[exp(-1/2*((t(i)-0.06)/0.03)^2),0,0,0,0,0,0,iin];
end

v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:750
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(18)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with noise with dt =.0013 ')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')



v3=0;
t = linspace(0,1,500);
dt=0.2;
f=1/.03;
F3=0;
for i=1:500
    iin=0.001*randn();
    F3(i,1:8)=[exp(-1/2*((t(i)-0.06)/0.03)^2),0,0,0,0,0,0,iin];
end

v3(1:8,1)=(C/dt+G)^-1 *(F3(1,:)');
for i=2:500
    v3(:,i)=(C/dt+G)^-1 *(C*v3(:,i-1)/dt+F3(i,:)');   
end

figure(19)
plot(t,v3(7,:))
hold on 
plot(t,F3(:,1))
title('Plot of guassian source with noise and dt=0.2')
legend('Vout','input signal')
xlabel('Time (s)')
ylabel('Voltage (V)')

% When the time step is increased, we see the shape of the output change.
% At higher time steps, the area of the guassian which goes below 0V
% disappears.  The plots are also less noisy than the plots with lower time
% steps

%% Question 4 
% In order to implement the non linearity, introduce a new vector for the
% non linear components.  With the new vector, the equation to solve for V
% must be include the extra b vector added to GX.
% to create the B vector, have a matrix where the V4 component is 
% $\beta*(V_3*G_3)^2 + \gamma (V_3*G_3)^3$
% the new equation is: 
% $$C\frac{dV}{dt}+GV+B=F$
% in order to do this we need to do a newton iteration at every step to
% make sure Vj converges at every time step.
