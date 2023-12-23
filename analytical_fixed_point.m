clear all
close all

%% parameters
del_x = 1.87; %length of track
del_t = 6.1; %time of track traversal
t_trial = del_t;
gamma = 3;
tau_p = .5; %time constant for LTP trace
tau_d = 1.5; %time constant for LTD trace
eta_p = .25; %activation constant for LTP trace
eta_d = 200; %activation constant for LTD trace
alpha = 1; %activation constant for LTD trace
T_max_p = 2.2; %maximum value for LTP trace
T_max_d = 2; %maximum value for LTD trace
T_0_p = 0; %basal level of LTP trace
T_0_d = 1.5; %basal level of LTD trace
T_max_py = T_max_p - T_0_p;
T_max_dy = T_max_d - T_0_d;
tau_P = .4; %time constant of instructive signal
dur_rf = 1; %duration of rectangular place field
t_1 = del_t/2 - dur_rf/2; %start of rectangular place field
t_2 = t_1 + dur_rf; %end of rectangular place field
t_0_list = 0:del_t/100:del_t;
[W,I_p_list, I_d_list] = deal(zeros(1,length(t_0_list)));
W_ind = 1;

%% analytical solution for 3 conditions (see appendix)
for t_0 = t_0_list
    if t_0 < t_1
        factor_p = gamma*eta_p*alpha*T_max_py/(1 + eta_p*alpha);
        factor_d = gamma*eta_d*alpha*T_max_dy/(1 + eta_d*alpha);
        term_1_p = tau_P*exp(t_0/tau_P)*(exp(-t_1/tau_P) - exp(-t_2/tau_P));
        term_1_d = tau_P*exp(t_0/tau_P)*(exp(-t_1/tau_P) - exp(-t_2/tau_P));
        factor_2_p = exp((1+eta_p*alpha)*t_1/tau_p + t_0/tau_P)/((1 + eta_p*alpha)/tau_p + 1/tau_P);
        factor_2_d = exp((1+eta_d*alpha)*t_1/tau_d + t_0/tau_P)/((1 + eta_d*alpha)/tau_d + 1/tau_P);
        term_2_p = exp(-((1 + eta_p*alpha)/tau_p + 1/tau_P)*t_2) - exp(-((1 + eta_p*alpha)/tau_p + 1/tau_P)*t_1);
        term_2_d = exp(-((1 + eta_d*alpha)/tau_d + 1/tau_P)*t_2) - exp(-((1 + eta_d*alpha)/tau_d + 1/tau_P)*t_1);
        factor_3_p = (exp(t_0/tau_P + ((1 + eta_p*alpha)*t_1/tau_p) - eta_p*alpha*t_2/tau_p) - exp(t_2/tau_p + t_0/tau_P))/(1/tau_p + 1/tau_P);
        factor_3_d = (exp(t_0/tau_P + ((1 + eta_d*alpha)*t_1/tau_d) - eta_d*alpha*t_2/tau_d) - exp(t_2/tau_d + t_0/tau_P))/(1/tau_d + 1/tau_P);
        term_3_p = exp(-(1/tau_p + 1/tau_P)*t_trial) - exp(-(1/tau_p + 1/tau_P)*t_2);
        term_3_d = exp(-(1/tau_d + 1/tau_P)*t_trial) - exp(-(1/tau_d + 1/tau_P)*t_2);
        term_4_p = gamma*T_0_p*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        term_4_d = gamma*T_0_d*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        
        I_p = factor_p*(term_1_p + factor_2_p*term_2_p + factor_3_p*term_3_p) + term_4_p;
        I_d = factor_d*(term_1_d + factor_2_d*term_2_d + factor_3_d*term_3_d) + term_4_d;
        
        W(W_ind) = I_p/(I_p + I_d);
        I_p_list(W_ind) = I_p;
        I_d_list(W_ind) = I_d;
        
    elseif (t_1 <= t_0) && (t_0 < t_2)
        factor_p = gamma*eta_p*alpha*T_max_py/(1 + eta_p*alpha);
        factor_d = gamma*eta_d*alpha*T_max_dy/(1 + eta_d*alpha);
        term_1_p = tau_P*(1 - exp(-(t_2-t_0)/tau_P));
        term_1_d = tau_P*(1 - exp(-(t_2-t_0)/tau_P));
        factor_2_p = exp((1+eta_p*alpha)*t_1/tau_p + t_0/tau_P)/((1 + eta_p*alpha)/tau_p + 1/tau_P);
        factor_2_d = exp((1+eta_d*alpha)*t_1/tau_d + t_0/tau_P)/((1 + eta_d*alpha)/tau_d + 1/tau_P);
        term_2_p = exp(-((1 + eta_p*alpha)/tau_p + 1/tau_P)*t_2) - exp(-((1 + eta_p*alpha)/tau_p + 1/tau_P)*t_0);
        term_2_d = exp(-((1 + eta_d*alpha)/tau_d + 1/tau_P)*t_2) - exp(-((1 + eta_d*alpha)/tau_d + 1/tau_P)*t_0);
        factor_3_p = (exp(t_0/tau_P + ((1 + eta_p*alpha)*t_1/tau_p) - eta_p*alpha*t_2/tau_p) - exp(t_2/tau_p + t_0/tau_P))/(1/tau_p + 1/tau_P);
        factor_3_d = (exp(t_0/tau_P + ((1 + eta_d*alpha)*t_1/tau_d) - eta_d*alpha*t_2/tau_d) - exp(t_2/tau_d + t_0/tau_P))/(1/tau_d + 1/tau_P);
        term_3_p = exp(-(1/tau_p + 1/tau_P)*t_trial) - exp(-(1/tau_p + 1/tau_P)*t_2);
        term_3_d = exp(-(1/tau_d + 1/tau_P)*t_trial) - exp(-(1/tau_d + 1/tau_P)*t_2);
        term_4_p = gamma*T_0_p*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        term_4_d = gamma*T_0_d*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        
        I_p = factor_p*(term_1_p + factor_2_p*term_2_p + factor_3_p*term_3_p) + term_4_p;
        I_d = factor_d*(term_1_d + factor_2_d*term_2_d + factor_3_d*term_3_d) + term_4_d;
        
        W(W_ind) = I_p/(I_p + I_d);
        I_p_list(W_ind) = I_p;
        I_d_list(W_ind) = I_d;
        
    elseif (t_2 <= t_0) && (t_0 < t_trial)
        factor_p = gamma*eta_p*alpha*T_max_py/(1 + eta_p*alpha);
        factor_d = gamma*eta_d*alpha*T_max_dy/(1 + eta_d*alpha);
        factor_3_p = (exp(t_0/tau_P + ((1 + eta_p*alpha)*t_1/tau_p) - eta_p*alpha*t_2/tau_p) - exp(t_2/tau_p + t_0/tau_P))/(1/tau_p + 1/tau_P);
        factor_3_d = (exp(t_0/tau_P + ((1 + eta_d*alpha)*t_1/tau_d) - eta_d*alpha*t_2/tau_d) - exp(t_2/tau_d + t_0/tau_P))/(1/tau_d + 1/tau_P);
        term_3_p = exp(-(1/tau_p + 1/tau_P)*t_trial) - exp(-(1/tau_p + 1/tau_P)*t_0);
        term_3_d = exp(-(1/tau_d + 1/tau_P)*t_trial) - exp(-(1/tau_d + 1/tau_P)*t_0);
        term_4_p = gamma*T_0_p*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        term_4_d = gamma*T_0_d*tau_P*exp(t_0/tau_P)*(exp(-t_0/tau_P) - exp(-t_trial/tau_P));
        
        I_p = factor_p*(factor_3_p*term_3_p) + term_4_p;
        I_d = factor_d*(factor_3_d*term_3_d) + term_4_d;
        
        W(W_ind) = I_p/(I_p + I_d);
        I_p_list(W_ind) = I_p;
        I_d_list(W_ind) = I_d;
    end
    W_ind = W_ind + 1;
end

%% plotting
figure
subplot(1,2,1)
plot(t_0_list-del_t/2,I_p_list*1000, 'r--')
hold on
plot(t_0_list-del_t/2,I_d_list*1000, 'b--')
xlabel('D')
ylabel('I_{p} and I_{d}')
legend('I_{p}','I_{d}')
xlim([-3.5 3.5])
xticks([-3:1:3])
% ylim([0 .65])

subplot(1,2,2)
plot(t_0_list-del_t/2,W, 'k-')
xlabel('D(s)')
ylabel('Fixed point value of W')
xlim([-3.5 3.5])
xticks([-3:1:3])
        
        