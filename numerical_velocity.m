clear all
close all

linear = 0; %set to one for linear, zero for circular track

dt = .0001; %timestep
R = .3;
L = 2*pi*R;

tau_p = .5; %time constant for LTP trace
T_max_p = 2.5; %maximum for LTP trace
eta_p = 1; %activation constant for LTP trace

tau_d = 1.5; %time constant for LTD trace
T_max_d = 2; %maximum for LTD trace
eta_d = 200; %activation constant for LTD trace

Sig_RF = .15;%.9/(3*sqrt(2)); %sigma for presynaptic gaussian
RF_center_x = R*pi;

A_p = 3; %magnitude of plateau potential
tau_I = .4; %time constant for plateau potential


V_arr = 0.075:0.1:1.05;

vv = 1;
vv_max = length(V_arr);
x_sample = 0:L/50:L;

W_mat = zeros(vv_max,length(x_sample));

width_arr = zeros(vv_max,3);

while vv <= vv_max
    V = V_arr(vv);
    max_time = L/V;
    st_I_x = .3;
    st_I_t = st_I_x/V;
    RF_center_t = RF_center_x/V;
    
    dx = dt*V;
    t_line = 0:dt:max_time;
    x_line = V*t_line;
    ii_size = size(t_line);
    ii_size = ii_size(2);
    U_arg = zeros(1,ii_size);
    U_arg_p = zeros(1,ii_size);
    U_arg_d = zeros(1,ii_size);
    x_C = Circ_dist(RF_center_x,x_line,R);
    
    ii = 2;
    
    while ii <= ii_size
        str = x_C(ii);
        U_arg(ii) = U_arg(ii-1) + dt*RF_C0(str);
        ii = ii+1;
    end
    
    U_arg_p = eta_p*U_arg+t_line;
    U_arg_d = eta_d*U_arg+t_line;
    
    U_p = exp(U_arg_p/tau_p);
    U_d = exp(U_arg_d/tau_d);
    
    H_p = zeros(1,ii_size);
    H_d = zeros(1,ii_size);
    
    jj = 2;
    
    while jj <= (ii_size)
        X_curr = x_C(jj);
        H_p(jj) = H_p(jj-1) + dt*eta_p*RF_C0(X_curr)*U_p(jj);
        H_d(jj) = H_d(jj-1) + dt*eta_d*RF_C0(X_curr)*U_d(jj);
        jj = jj+1;
    end
    
    T_p1 = H_p./U_p;
    T_p1 = T_p1*T_max_p/tau_p;
    
    T_d1 = H_d./U_d;
    T_d1 = T_d1*T_max_d/tau_d;
    
    yp_t = P(t_line,A_p,st_I_t,tau_I,max_time);
    
    if linear == 1
        C_p = 0;%T_p1(end)/(1 - U_p(end)^(-1));
        C_d = 0;%T_d1(end)/(1 - U_d(end)^(-1));
    else
        C_p = T_p1(end)/(1 - U_p(end)^(-1));
        C_d = T_d1(end)/(1 - U_d(end)^(-1));
    end
    
    T_p = T_p1 + C_p./U_p;
    T_d = T_d1 + C_d./U_d;
    
    
    down_sample_factor = 1000;
    D_sample = x_sample - RF_center_x;
    D_sample_t = D_sample/V;
    
    D_size = size(D_sample);
    D_size = D_size(2);
    
    I_p_arr = zeros(1,D_size);
    I_d_arr = zeros(1,D_size);
    W_fp = zeros(1,D_size);
    
    dd = 1;
    
    while dd <= D_size
        D = D_sample_t(dd);
        st_I_t = D + RF_center_t;
        I_p = 0;
        I_d = 0;
        
        ii = 2;
        while ii <= ii_size
            I_p = I_p + dt*T_p(ii)*P(t_line(ii),A_p,st_I_t,tau_I,max_time);
            I_d = I_d + dt*T_d(ii)*P(t_line(ii),A_p,st_I_t,tau_I,max_time);
            ii = ii + 1;
        end
        I_p_arr(dd) = I_p;
        I_d_arr(dd) = I_d;
        W_fp(dd) = I_p/(I_p+I_d);
        dd = dd + 1;
    end
    
if vv == 1
    figure
    plot(D_sample,W_fp)
    hold on
else
    hold on
    plot(D_sample,W_fp)
    hold on
end

vv = vv + 1;
    
end

legend



function y = Circ_dist(RF_center,x_line,R)
    theta1=RF_center ./ R;
    theta2=x_line ./ R;
    theta=(((theta2 - theta1) < pi) .* (theta2 - theta1) +  ((theta2 - theta1) > pi) .* (theta2 - theta1 -2*pi ) );
    y = (R .* theta);
end

function y = RF_C0(x)
    sig_R = .15;
    y = exp((-(x.^2 )./(2*(sig_R ^2))));
end


function y = P_circ(x,A_p,st,tau,L)
y = (x>st)*A_p.*exp(-(x-st)/tau)...
    + (x <= st)*A_p.*exp(-(L-st+x)/tau);
end

function y = P(x,A_p,st,tau,L)
y = (x>st)*A_p.*exp(-(x-st)/tau);
end