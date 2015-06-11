% To run, save as binaryEx_synch.m or copy and paste in a MatLab m-file.
% Change the variables InitCond or t_stop if desired.  Save and run code 
% (under Debug menu).

function binaryEx_synch

% This function runs the binary model synchronously
%
% Species Identifiers:
% (1)egf (2)hrg (3)e    gfr (4)raf (5)pi3k (6)erk (7)akt.
% To run, press F5 or choose Debug-> save and run.
%
% The following variables can be altered:
% InitCond: initial condition of each species given in species identifiers
% above
% t_stop: number of time steps to run simulation.


InitCond = [1 1 0 0 0 0 0];
t_stop = 10;

% Do not alter below

% set the initial condition
y(1,:) = InitCond;
step = 1;
% this loop calculates the values of each species during
while step < t_stop +1
    % the following line calls the rules specified below.
    y(step+1,:) = Rules(y(step,:));
    fprintf('%d',step +1) 
    y(step+1,:) 
    step = step +1;
end

% plot timecourse of Erk and Akt
figure('Name','Binary_synch')

plot((0:1:t_stop),y(:,6),'k-o','LineWidth',2,'MarkerSize',8)
hold on
plot((0:1:t_stop),y(:,7),'r:d','LineWidth',2,'MarkerSize',5)
set(gca,'fontsize',18)
xlabel('TimeStep','fontsize',30)
ylabel('Value of Species','fontsize',30)
legend({'Erk','Akt'},'fontsize',30)



%%Rules
function y = Rules(x)
y = x;

%egfr
if x(1) + x(2) > 0
    y(3) = 1;
else
    y(3) = 0;
end

%raf
if x(3) + x(7) > 0
    y(4) = 1;
else
    y(4) = 0;
end

%pi3k
if x(3) - x(6) > 0
    y(5) = 1;
else
    y(5) = 0;
end

%akt
if x(5) > 0
    y(7) = 1;
else
    y(7) = 0;
end

%erk
if x(4) > 0
    y(6) = 1;
else
    y(6) = 0;
end
