clear all
close all

%% parameters
N = 50; %number of presynaptic inputs
del_x = 1.87; %length of track (meters)
del_t = 16.1; %traversal time of track (seconds)
dt = 1*10^(-3); %timestep 
v = del_x/del_t; %velocity
t = 0:dt:del_t - dt*del_t; %time vector
D_vect_x = (-N/2:N/2 - 1)*(del_t/N);
D_vect_t = (t-del_t/2);
x = v*t; %place vector
sigma = .9/(3*sqrt(2)); %sigma for presynaptic gaussian
alpha = 7.5;

t2 = 0:dt/2:del_t - dt*del_t/2; %copy for presynaptic creation
x2 = v*t2; %copy for presynaptic creation
sigma2 = sigma/2; %copy for presynaptic creation


[R_1] = deal(zeros(N,length(x))); %presynaptic input array
[R_i] = deal(zeros(2*N,length(x2))); %copy for presynaptic creation
[LTP,LTD,LTP_W,LTD_W] = deal(zeros(N,length(t))); %trace arrays. weight modified trace arrays
[IP,ID] = deal(zeros(1,N)); %overlap vectors

tau_p = .2; %time constant for LTP trace
tau_d = 1.5; %time constant for LTD trace
eta_p = .2; %activation constant for LTP trace
eta_d = 200; %activation constant for LTD trace
eta_W = .0006; %learning rate
max_p = 2.2; %maximum for LTP trace
max_d = 2; %maximum for LTD trace
tau_P = .4; %time constant for plateau potential
p_mag = 1; %magnitude of plateau potential

num_laps = 10; %number of laps
num_tests = 1; %number of simulations (different parameters can be chosen before each test
circular = 2; %set to 1 for circular condition, 2 for linear

tr_loc = round((N)/2); %location of traces for plotting


%% creation of presynaptic inputs
for i = 1:2*N
    R_i(i,:) = circshift(exp(-((x2-(del_x/2))).^2/(2*(sigma2^2))),i*round((length(x2)/(2*N))));
end

for i = 1:N
    R_1(i,:) = circshift(exp(-((x-(del_x/2))).^2/(2*(sigma^2))),i*round((length(x)/(N))));
end

R_i2 = [R_i(N + 1:2*N,:);R_i(1:N,:)];
R_i2 = R_i2((2*N)/4:3*(2*N)/4,round(length(x2)/4):round(3*length(x2)/4));
R_i2 = R_i2(1:50,1:length(R_1));
R_2 = [R_1(N/2 + 1:N,:);R_1(1:N/2,:)];

if circular == 1
    R_i= R_2; %use for circular
else
    R_i = R_i2; %use for linear
end


%% main program

for test = 1:num_tests
    W_center = 4; %center of existing postsynaptic pf
    W = 0 + .0*exp(-(((1:N) - W_center).^2)/65);% initial weights

    p_loc = (N/2)*(del_x/N); %location of plateau potential
    p_time = p_loc/v; %time of pp
    p_t_ind = round(p_time/dt); %index of pp
    
    V_t = zeros(num_laps,length(t)); %ramp amplitude array
    V_t(1,:) = alpha*(W*R_i); %initial ramp amplitude
    
    P = zeros(1,length(t)); %instructive signal vector
    
    for lap = 2:num_laps
        dW = 0;
        for i = circular:length(t) 
            if i == 1
                t_minus = length(t);
            else
                t_minus = i-1;
            end
            t_val = t(t_minus);
            x_val = v*t_val;
            x_minus = find(x == x_val);
            dLTP = (-LTP(:,t_minus) + eta_p*R_i(:,x_minus).*(max_p-LTP(:,t_minus)))*(dt/tau_p);
            dLTD = (-LTD(:,t_minus) + eta_d*R_i(:,x_minus).*(max_d-LTD(:,t_minus)))*(dt/tau_d);
            LTP(:,i) = LTP(:,t_minus) + dLTP;
            LTD(:,i) = LTD(:,t_minus) + dLTD;
            LTP_W(:,i) = LTP(:,i).*(1-W');
            LTD_W(:,i) = LTD(:,i).*(W');

            if i == round(p_time/dt) %time of plateau potential
                P(i) = p_mag;
            else
                P(i) = P(t_minus) - P(t_minus)*(dt/tau_P);
            end

            dW = dW + eta_W*P(i)*(LTP_W(:,i)-LTD_W(:,i)); %weight update
        end
        if lap == 3 && test == num_tests
            IP = LTP*P';
            ID = LTD*P';
            W_fixed = (IP./(IP + ID))'; %numerical fixed point
            V_new = alpha*W_fixed*R_i; %fixed point ramp amplitude
        end
        W = W + dW';
        W(W<0) = 0;
        
        V_t(lap,:) = alpha*(W*R_i); %ramp amplitude at lap "lap"
    end
end

%% plotting
figure
subplot(1,3,1)
imagesc(R_i)
title('Presynaptic inputs')
xlabel('time (ms)')
ylabel('Neuron number')

subplot(1,3,2)
plot(D_vect_x,W_fixed)
title('Fixed point W')
xlabel('D(s)')
ylabel('synaptic strength')
xlim([-8 8])

subplot(1,3,3)
plot(D_vect_t,V_new)
hold on
plot(D_vect_t,V_t(1,:), 'r--')
hold on
plot(D_vect_t,V_t(2,:), 'b--')
hold on
plot(D_vect_t,V_t(4,:), 'g--')
hold on
plot(D_vect_t,V_t(6,:), 'y--')
hold off
title('ramp amplitude')
xlabel('D(s)')
ylabel('Ramp amplitude (mV)')
xlim([-8 8])