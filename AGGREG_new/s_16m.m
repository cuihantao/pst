% m.file to simulate the 16 machine, 68 bus system in 
%   LOTDYS report using the Matlab Power System Toolbox
jay = sqrt(-1);

pst_var % set up global variable 
data16m % load input data from m.file
basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = 0 ;     % machine 1 is reference

while(0)
% solve for loadflow - loadflow parameter
tol = 1e-9;   % tolerance for convergence
iter_max = 30; % maximum number of iterations
vmin = 0.5;   % voltage minimum
vmax = 1.5;   % voltage maximum
acc = 1.0;   % acceleration factor
[bus_sol,line_flw] = ...
       loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2);
bus = bus_sol;  % solved loadflow solution needed for
                 % initialization
save data16m bus line mac_con
end %while(0)
load data16m

% simulation
t_switch(1) = 0;     % all time in second+s, start time
t_switch(2) = 0.04;  % time to apply fault
t_switch(3) = 0.04; % time to clear fault, about 4.8 cycles
t_switch(4) = 2.00;  % end time
h = 0.01; % integration stepsize
stepsize = h;
k_switch(1) = round((t_switch(2)-t_switch(1))/stepsize)+1;
k_switch(2) = round((t_switch(3)-t_switch(1))/stepsize)+1;
k_switch(3) = round((t_switch(4)-t_switch(1))/stepsize)+1;

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
[Y_red,V_rec] = red_ybus(bus,line);             % pre-fault 
                                         % admittance matrix
% create bus matrix with short circuit on bus 29
bus_f = bus;
bus_f(29,6) = 100000.;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);   % fault-on
						    % admittance matrix
line_pf = line;
% remove line 28-29
line_pf(45,:)=[28 29 0 1000. 0 1. 0 ]; % high line impedance
[Y_red_pf,V_rec_pf] = red_ybus(bus,line_pf);  % post-fault
                                    % admittance matrix

% step 2: initialization
flag = 0;
f = mac_em(0,1,bus,flag);    % all machine electro-
                             % mechanical model
keyboard
% step 3: perform a predictor-corrector integration 
for k = 1:k_switch(3)+1
  % step 3a: network solution
  % mach_ref(k) = mac_ang(syn_ref,k);
  mach_ref(k) = 0;
  flag = 1;
  f = mac_em(0,k,bus,flag); % network-machine interface
  psi = psi_re(:,k) + jay*psi_im(:,k); 
  if k >= k_switch(2)
     cur = Y_red_pf*psi; % network solution 
     v = V_rec_pf*psi;   % bus voltage reconstruction
  elseif k >= k_switch(1)
     cur = Y_red_f*psi; % network solution 
     v = V_rec_f*psi;   % bus voltage reconstruction
  else
     cur = Y_red*psi; % network solution 
     v = V_rec*psi;   % bus voltage reconstruction
  end
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
  % step 3b: compute dynamics and integrate
  flag = 2;
  pmech(:,k) = pmech(:,1); % constant mechanical input power
  f = mac_em(0,k,bus,flag); % dynamics calculation
  if k ~=k_switch(3)+1
    j = k+1;
    % following two statements are predictor steps
    mac_ang(:,j) = mac_ang(:,k) + h*dmac_ang(:,k); 
    mac_spd(:,j) = mac_spd(:,k) + h*dmac_spd(:,k);
    edprime(:,j) = edprime(:,k) + h*dedprime(:,k);
    eqprime(:,j) = eqprime(:,k) + h*deqprime(:,k);
    flag = 1;
    % mach_ref(j) = mac_ang(syn_ref,j);
    mach_ref(j) = 0;
    f = mac_em(0,j,bus,flag); 
    psi = psi_re(:,j) + jay*psi_im(:,j);
    if k >= k_switch(2)
      cur = Y_red_pf*psi;
    elseif k >= k_switch(1)
      cur = Y_red_f*psi;
    else
      cur = Y_red*psi;
    end
    cur_re(:,j) = real(cur); cur_im(:,j) = imag(cur);
    pmech(:,j) = pmech(:,k);
    flag = 2;
    f = mac_em(0,j,bus,flag);
    % following two statements are corrector steps
    mac_ang(:,j) = mac_ang(:,k) +..
                   h*(dmac_ang(:,k)+dmac_ang(:,j))/2.;
    mac_spd(:,j) = mac_spd(:,k) +..
                   h*(dmac_spd(:,k)+dmac_spd(:,j))/2.;
    edprime(:,j) = edprime(:,k) +..
                   h*(dedprime(:,k)+dedprime(:,j))/2.;
    eqprime(:,j) = eqprime(:,k) +..
                   h*(deqprime(:,k)+deqprime(:,j))/2.;
  end
end

t = [0:stepsize:t_switch(4)]'; % time
