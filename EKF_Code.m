clear variables; 
close all;
clc;
%% Extended Kalman Filter Matlab

load('system_model','-mat','Times','HorizontalPositionm','h','fn','time','data','h10','h_20')

%% Set up time step and parameters

del_t = 0.005;            %time step
runT = 5;              %run time
nT = runT/del_t;    %number of samples
t = 0:del_t:runT;

n = 2;                  %number of states
nz = 3;                 %number of measurements

%% Constructing State Vectors
x = zeros(n,nT);        %x_0 = [0;0]
%x(1,1) = 20;
%x(2,1) = 2;
P = zeros(n,n,nT);      %covariance matrix
z = zeros(nz,nT);    %measurement vector

for i = 1:nT
    P(:,:,i) = eye(n);
end 

%% System Model
syms d

Q = diag([0.050 0.5]);
R = diag([0.0880 0.1018 0.1018]);
G = [1 del_t;
     0 1];

%% Sensor Model
H = zeros(nz,n,nT);
h = zeros(nz,nT);

% Sensor Equations
h_fun20= h_20(1,1)*d^3 + h_20(1,2)*d^2 + h_20(1,3)*d + h_20(1,4);
h_fun10 = h10(1,1)*d^3 + h10(1,2)*d^2 + h10(1,3)*d + h10(1,4);
    
H_fun20 = diff(h_fun20);
H_fun10 = diff(h_fun10);

%% Setting up Matrices
xpred = x;
Ppred = P;
z(1,:) = data_1(:,2)';
z(2,:) = data_1(:,4)';
z(3,:) = data_1(:,6)';


%% Simulated Sensor
 
 xsim(1:n,1:nT)=0;
 xsim(1,1) = 20;
 xsim(2,1)=2;
 
 for k = 2:nT
    
    xsim(:,k) = [20+2*t(k-1);
                  2;];
               
     xsim(:,k)= xsim(:,k) + mvnrnd(zeros(n,1),Q,1)';
    
     z(:,k) = [double(subs(h_fun20,d,xsim(1,k)));
              double(subs(h_fun10,d,xsim(1,k)));
               double(subs(h_fun10,d,xsim(1,k)));];
     z(:,k)= z(:,k)+ mvnrnd(zeros(nz,1), R, 1)'; 
     
     
 end

%% Begin simulation
for k = 2:nT
    
    %Predict
    %xpred(:,k) = double(subs(g,t,timevec(k-1))); 
    %Ppred(:,:,k) = double(subs(G,t,timevec(k-1)))*P(:,:,k-1)*double(subs(G,t,timevec(k-1))) + Q;
    xpred(:,k) = [x(1,k-1) + x(2,k-1)*del_t;
                    x(2,k-1)];
    Ppred(:,:,k) = G*P(:,:,k-1)*G' + Q;
    
    %Calculate H and h   
    h(1,k) = double(subs(h_fun20,d,xpred(1,k)));
    h(2,k) = double(subs(h_fun10,d,xpred(1,k)));
    h(3,k) = double(subs(h_fun10,d,xpred(1,k)));
    
    H(1,1,k) = double(subs(H_fun20,d,xpred(1,k)));
    H(2,1,k) = double(subs(H_fun10,d,xpred(1,k)));
    H(3,1,k) = double(subs(H_fun10,d,xpred(1,k)));
    
    %Update - assume z has been set
    K = Ppred(:,:,k)*H(:,:,k)'*inv((H(:,:,k)*Ppred(:,:,k)*H(:,:,k)' + R));
    x(:,k) = xpred(:,k) + K*(z(:,k) - h(:,k)); %z(:,k) is either actual measurements or simulated measurements
    P(:,:,k) = (eye(n) - K*H(:,:,k)*Ppred(:,:,k));
    
end

% %% Phone Motion
% vec10 = zeros(1,115);
% time10 = zeros(1,115);
% for i = 1:115
%     vec10(i,1) = 10;
%     time10(i,1) = 5.78;
% end
% 
% HorizontalPositionm = HorizontalPositionm + vec10;
% Times = Times - time10;

%% Plots
hold on
figure(1)
title('Distance vs Time with Object Motion Data - Poor Sensing Conditions')
xlabel('Time(s)')
ylabel('Distance(cm)')
plot(t_fil,x(1,:),'DisplayName','Extended Kalman Filter')
%plot(Times(61:115),HorizontalPositionm(61:115),'DisplayName','Phone Motion App')
legend('show')
%plot(t,xsim(1,:),'DisplayName','Simulated Model')

figure(2)
title('Velocity vs Time with Object Motion Data - Poor Sensing Conditions')
xlabel('Time(s)')
ylabel('Velocity(cm/s)')
hold on
%plot(t,xsim(2,:),'DisplayName','Simulated Model')
plot(t_fil,x(2,:),'DisplayName','Extended Kalman Filter')
legend('show')





