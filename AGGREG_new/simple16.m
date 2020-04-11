% simple16.m
% hybrid3.m simulates the 50 machine 145 bus system from Iowa State University
% This routine simulates the system using the Matlab Power System Toolbox
% and swing.m before and after the fault is clear, respectively.
clear
jay = sqrt(-1);

% set up global variables
mac_var

basrad = 2*pi*60; % system frequency is 60 Hz
basmva = 100;     % 100 MVA base
syn_ref = -1;     % absolute machine angle (zero reference)

input('time to clear fault: (sec)     ');
t_switch(3) = ans;
input('simulation time: (sec)     ');
t_switch(4) = ans;

lfw = input('Do you need to recompute a loadflow solution? [y/n]   ','s');
% prjt = input('Do you want initial condition projection? [y/n]   ','s');
%if prjt == 'y',
%  ptype = input('projection:linear(1), partial nonlin.(2), nonlin.(3) or heavy integ(4)   ');
%end

if lfw == 'y',
  data16m     % load input data from m.file
  %% solve for loadflow
  % loadflow parameter
  tol = 1e-9;   % tolerance for convergence
  iter_max = 30; % maximum number of iterations
  vmin = 0.5;   % voltage minimum
  vmax = 1.5;   % voltage maximum
  acc = 1.0;   % acceleration factor %
  [bus_sol,line_flw] = loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2); %%
  bus = bus_sol;  % loadflow solution needed for initialization
  save data16m bus line mac_con asyn_damp
else
  load data16m
end
[nbus dummy] = size(bus);
[nmach dummy] = size(mac_con);

n = 5;     % number of groups
%U = form_u(bus,line,n)     % grouping matrix

%load gpmat
%UU = U;    % prepare the global grouping matrix

bus_display = [1 8 9];     % buses to be displayed
% construct bus voltage reduction matrix
bus_red = zeros(length(bus_display),nbus);
for i = 1:length(bus_display)
  bus_red(i,bus_display(i)) = 1;
end; %

% simulation
t_switch(1) = 0;      % all time in seconds, start time
t_switch(2) = 0.0;    % time to apply fault
%t_switch(3) = 0.23;   % time to clear fault, about 5 cycles
%t_switch(4) = 5.;    % end time
tol = 0.0005;

% step 1: construct reduced Y matrix - all loads are assumed to be
%         of constant impedance load
while(0)
[Y_red,V_rec] = red_ybus(bus,line);                 % pre-fault
                                                    % admittance matrix
end % while(0)
% create a short circuit at bus 29
%while(0)
bus_f = bus;
%bus_f(29,6) = 100000.;
bus_f(28,6) = 100000.;
[Y_red_f,V_rec_f] = red_ybus(bus_f,line);           % fault-on
						    % admittance matrix
% remove the line connecting buses 15 (28) and 16 (29)
bus_pf = bus;
line_pf = line;
line_pf(length(line)+1,:) = [28 29 -0.0014 -0.0151 -0.2490 0. 0.];
[Y_red_pf,V_rec_pf] = red_ybus(bus_pf,line_pf);        % post-fault
						    % admittance matrix
save yred Y_red_f V_rec_f Y_red_pf V_rec_pf
%end %while(0)
%load yred
% step 2: initialization
flag = 0;
f = mac_em(0,1,bus,flag);   % all machine electro-mechanical model
% keyboard

% step 3: perform integration
k = 1;
pow = 1/3;
t = t_switch(1);
hmax = (t_switch(4) - t)/5.;
hmin = (t_switch(4) - t)/20000.;
h = (t_switch(4) - t)/100.;
y = [mac_ang(:,1);mac_spd(:,1)];
tout = t;
tau = tol*max(norm(y, 'inf'), 1);

while (t < t_switch(3)) & (h >=hmin)
  if t + h > t_switch(4), h = t_switch(4) - t; end
  if (t < t_switch(3) & t + h > t_switch(3)), h = t_switch(3) - t; end

  % step 3a: network solution
    mach_ref(k) = 0;
  flag = 1;
  f = mac_em(0,k,bus,flag);
  psi = psi_re(:,k) + jay*psi_im(:,k);
  cur = Y_red_f*psi;
  v = V_rec_f*psi;
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
  bus_v(:,k) = abs(bus_red*v); 
  bus_ang(:,k) = angle(bus_red*v)*180/pi;
  % step 3b: compute dynamics and integrate
  pmech(:,k) = pmech(:,1);
  f = mac_em(0,k,bus,2);
% keyboard
  s1 = [dmac_ang(:,k);dmac_spd(:,k)];
  mac_ang(:,k+1) = mac_ang(:,k)+h*dmac_ang(:,k);
  mac_spd(:,k+1) = mac_spd(:,k)+h*dmac_spd(:,k);
  if syn_ref == 0
     % weighted average of machine angles
     mach_ref(k+1) = sum(mac_ang(:,k+1).*mac_con(:,16)./mac_pot(:,1))/H_sum;
  elseif syn_ref == -1
     mach_ref(k+1) = 0;
  else
     mach_ref(k+1) = mac_ang(syn_ref,k+1);
  end
  edprime(:,k+1)=edprime(:,k);
  eqprime(:,k+1)=eqprime(:,k);

 % fprintf('angle = %g, speed = %g. \n',  mac_ang(2,k+1), mac_spd(2,k+1))
  f = mac_em(0,k+1,bus,1);
  psi = psi_re(:,k+1) + jay*psi_im(:,k+1);
  cur = Y_red_f*psi;
  v = V_rec_f*psi;
  cur_re(:,k+1) = real(cur); cur_im(:,k+1) = imag(cur);
  bus_v(:,k+1) = abs(bus_red*v); 
  bus_ang(:,k+1) = angle(bus_red*v)*180/pi;

  pmech(:,k+1) = pmech(:,1);
  f = mac_em(0,k+1,bus,2);
  %fprintf('t = %g,  dangle = %g, dspeed = %g. \n', t, dmac_ang(2,k+1), dmac_spd(2,k+1))
  %pause
  s2 = [dmac_ang(:,k+1);dmac_spd(:,k+1)];
  mac_ang(:,k+2) = mac_ang(:,k)+h*(dmac_ang(:,k)+dmac_ang(:,k+1))/4;
  mac_spd(:,k+2) = mac_spd(:,k)+h*(dmac_spd(:,k)+dmac_spd(:,k+1))/4;
  if syn_ref == 0
    % weighted average of machine angles
    mach_ref(k+2) = sum(mac_ang(:,k+2).*mac_con(:,16)./mac_pot(:,1))/H_sum;
  elseif syn_ref == -1
      mach_ref(k+2) = 0;
  else
      mach_ref(k+2) = mac_ang(syn_ref,k+2);
  end
  edprime(:,k+2) = edprime(:,k);
  eqprime(:,k+2) = eqprime(:,k);

  f = mac_em(0,k+2,bus,1);
  psi = psi_re(:,k+2) + jay*psi_im(:,k+2);
  cur = Y_red_f*psi;
  v = V_rec_f*psi;
  cur_re(:,k+2) = real(cur); cur_im(:,k+2) = imag(cur);
  bus_v(:,k+2) = abs(bus_red*v); bus_ang(:,k+2) = angle(bus_red*v)*180/pi;
  pmech(:,k+2) = pmech(:,1);
  f = mac_em(0,k+2,bus,2);
  s3 = [dmac_ang(:,k+2);dmac_spd(:,k+2)];

  if t < t_switch(4)
      % Estimate the error and the acceptable error
      delta = norm(h*(s1 - 2*s3 + s2)/3,'inf');
      tau = tol*max(norm(y,'inf'),1.0);

     % Update the solution only if the error is acceptable
      if delta <= tau
	 t = t+h;
	 tout = [tout;t];
	 mac_ang(:,k+1) = mac_ang(:,k) + h*(dmac_ang(:,k) + 4*dmac_ang(:,k+2)...
			   + dmac_ang(:,k+1))/6;
	 mac_spd(:,k+1) = mac_spd(:,k) + h*(dmac_spd(:,k) + 4*dmac_spd(:,k+2)...
			   + dmac_spd(:,k+1))/6;
	 k = k + 1;
	% fprintf('t = %g,  angle = %g, speed = %g. \n', t, mac_ang(2,k), mac_spd(2,k))
	 fprintf('t = %g \n', t)
      end
    
      % Update the step size
      if delta ~= 0.0
         h = min(hmax, 0.9*h*(tau/delta)^pow);
      end
  end
end

%keyboard

t0 = t; tf = t_switch(4);
while (t < t_switch(4)) & (h >=hmin)
  if t + h > t_switch(4), h = t_switch(4) - t; end

  % step 3a: network solution
  if syn_ref == 0
    % weighted average of machine angles
    mach_ref(k) = sum(mac_ang(:,k).*mac_con(:,16)./mac_pot(:,1))/H_sum;
  elseif syn_ref == -1
    mach_ref(k) = 0;
  else
    mach_ref(k) = mac_ang(syn_ref,k);
  end
  flag = 1;
  f = mac_em(0,k,bus,flag);
  psi = psi_re(:,k) + jay*psi_im(:,k);
  cur = Y_red_pf*psi;
  v = V_rec_pf*psi;
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
  bus_v(:,k) = abs(bus_red*v); 
  bus_ang(:,k) = angle(bus_red*v)*180/pi;
  % step 3b: compute dynamics and integrate
  pmech(:,k) = pmech(:,1);
  f = mac_em(0,k,bus,2);
  s1 = [dmac_ang(:,k);dmac_spd(:,k)];
  mac_ang(:,k+1) = mac_ang(:,k)+h*dmac_ang(:,k);
  mac_spd(:,k+1) = mac_spd(:,k)+h*dmac_spd(:,k);
  if syn_ref == 0
     % weighted average of machine angles
     mach_ref(k+1) = sum(mac_ang(:,k+1).*mac_con(:,16)./mac_pot(:,1))/H_sum;
  elseif syn_ref == -1
     mach_ref(k+1) = 0;
  else
     mach_ref(k+1) = mac_ang(syn_ref,k+1);
  end
  edprime(:,k+1)=edprime(:,k);
  eqprime(:,k+1)=eqprime(:,k);

 % fprintf('angle = %g, speed = %g. \n',  mac_ang(2,k+1), mac_spd(2,k+1))
  f = mac_em(0,k+1,bus,1);
  psi = psi_re(:,k+1) + jay*psi_im(:,k+1);
  cur = Y_red_pf*psi;
  v = V_rec_pf*psi;
  cur_re(:,k+1) = real(cur); cur_im(:,k+1) = imag(cur);
  bus_v(:,k+1) = abs(bus_red*v); 
  bus_ang(:,k+1) = angle(bus_red*v)*180/pi;

  pmech(:,k+1) = pmech(:,1);
  f = mac_em(0,k+1,bus,2);
  %fprintf('t = %g,  dangle = %g, dspeed = %g. \n', t, dmac_ang(2,k+1), dmac_spd(2,k+1))
  %pause
  s2 = [dmac_ang(:,k+1);dmac_spd(:,k+1)];
  mac_ang(:,k+2) = mac_ang(:,k)+h*(dmac_ang(:,k)+dmac_ang(:,k+1))/4;
  mac_spd(:,k+2) = mac_spd(:,k)+h*(dmac_spd(:,k)+dmac_spd(:,k+1))/4;
  if syn_ref == 0
    % weighted average of machine angles
    mach_ref(k+2) = sum(mac_ang(:,k+2).*mac_con(:,16)./mac_pot(:,1))/H_sum;
  elseif syn_ref == -1
      mach_ref(k+2) = 0;
  else
      mach_ref(k+2) = mac_ang(syn_ref,k+2);
  end
  edprime(:,k+2) = edprime(:,k);
  eqprime(:,k+2) = eqprime(:,k);

  f = mac_em(0,k+2,bus,1);
  psi = psi_re(:,k+2) + jay*psi_im(:,k+2);
  cur = Y_red_pf*psi;
  v = V_rec_pf*psi;
  cur_re(:,k+2) = real(cur); cur_im(:,k+2) = imag(cur);
  bus_v(:,k+2) = abs(bus_red*v); bus_ang(:,k+2) = angle(bus_red*v)*180/pi;
  pmech(:,k+2) = pmech(:,1);
  f = mac_em(0,k+2,bus,2);
  s3 = [dmac_ang(:,k+2);dmac_spd(:,k+2)];

  if t < t_switch(4)
      % Estimate the error and the acceptable error
      delta = norm(h*(s1 - 2*s3 + s2)/3,'inf');
      tau = tol*max(norm(y,'inf'),1.0);

     % Update the solution only if the error is acceptable
      if delta <= tau
	 t = t+h;
	 tout = [tout;t];
	 mac_ang(:,k+1) = mac_ang(:,k) + h*(dmac_ang(:,k) + 4*dmac_ang(:,k+2)...
			   + dmac_ang(:,k+1))/6;
	 mac_spd(:,k+1) = mac_spd(:,k) + h*(dmac_spd(:,k) + 4*dmac_spd(:,k+2)...
			   + dmac_spd(:,k+1))/6;
	 k = k + 1;
	% fprintf('t = %g,  angle = %g, speed = %g. \n', t, mac_ang(2,k), mac_spd(2,k))
	 fprintf('t = %g \n', t)
      end
    
      % Update the step size
      if delta ~= 0.0
         h = min(hmax, 0.9*h*(tau/delta)^pow);
      end
  end
end
% save results
%save result2 tout mac_ang mac_spd aggre_var diff_var

mac_angd = mac_ang(:,1:k)*180/pi;
mac_a = mac_angd - ones(16,1)*mac_angd(16,:);
mac_b = [mac_a(1:3,:); mac_a(5:6,:)];
spd_b = [mac_spd(1:3,1:k); mac_spd(1:3,1:k)];
save allvar
