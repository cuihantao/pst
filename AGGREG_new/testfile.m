% load data16m
load iamed17
%[bus_sol,line_flow] = loadflow(bus,line,1e-10,5,.5,1.5,1,..
%                               'n',2);
%while(0)
bus(1,1)= 300;
line(1,1) = 300;
line(2,2) = 300;
line(6,2) = 300;
%end
bus_list=[2:1:154]';
bus_list=[2:1:29, 300]';
flag = 0; % eliminate
flag=1;  % retain
% tol=10e-8;
[bus_red,line_red]=reduce(bus,line,bus_list,flag);
 keyboard
% save temp bus_red line_red
 [bus_sol,line_flow] = loadflow(bus_red,line_red, ..
                    1e-10,5,.5,1.5,1,'n',2);
