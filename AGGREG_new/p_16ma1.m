%  p_16ma1.m
% m.file to for the computing perturbations
%   for sv_16ma1.m

  % step 3a: network solution
  flag = 1;
  f = mac_tra(0,k,bus,flag);                         % machine 1
  f = exc_dc12(0,k,bus,flag,stepsize);               % machine 1 exciter
  psi = psi_re(:,k) + jay*psi_im(:,k);               % generator voltage
  V_nc = nc_load(bus,flag,Y22,Y21,psi,V_nc,1e-14);
  cur = Y11*psi + Y12*V_nc;                          % generator current
  V_b = rec_V1*psi + rec_V2*V_nc;                    % bus voltage
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
%  v(bus_ord,1) = [V_nc; V_b];
% step 3b: compute dynamics 
  pmech(:,k) = pmech(:,1);
  exc_sig(:,k) = zeros(9,1);
%  vex(1,k) = vex(1,1);
  flag = 2;
  f = mac_tra(0,k,bus,flag);                         % machine 1
  f = exc_dc12(0,k,bus,flag,stepsize);               % machine 1 exciter
if j <= 63 
  a(1:9,j)   = dmac_ang(:,1) / pert;                 % building A matrix
  a(10:18,j) = dmac_spd(:,1) / pert;
  a(19:27,j) = deqprime(:,1) / pert;
  a(28:36,j) = dedprime(:,1) / pert;
  a(37:45,j) = dEfd(:,1) / pert;
  a(46:54,j) = dV_R(:,1) / pert;
  a(55:63,j) = dR_f(:,1) / pert;
  c(1,j) = ( abs  ( V_nc(1,1) ) - y1 ) / pert;       % building C matrix
  c(2,j) = ( angle( V_nc(1,1) ) - y2 ) / pert;
  c(3,j) = ( abs  ( V_nc(2,1) ) - y3 ) / pert;
  c(4,j) = ( angle( V_nc(2,1) ) - y4 ) / pert;
 else
     b(1:9,ii)   = dmac_ang(:,1) / pert;             % building B matrix
     b(10:18,ii) = dmac_spd(:,1) / pert;
     b(19:27,ii) = deqprime(:,1) / pert;
     b(28:36,ii) = dedprime(:,1) / pert;
     b(37:45,ii) = dEfd(:,1) / pert;
     b(46:54,ii) = dV_R(:,1) / pert;
     b(55:63,ii) = dR_f(:,1) / pert;
     d(1,ii) = ( abs  ( V_nc(1,1) ) - y1 ) / pert;   % building D matrix
     d(2,ii) = ( angle( V_nc(1,1) ) - y2 ) / pert;
     d(3,ii) = ( abs  ( V_nc(2,1) ) - y3 ) / pert;
     d(4,ii) = ( angle( V_nc(2,1) ) - y4 ) / pert;
end
