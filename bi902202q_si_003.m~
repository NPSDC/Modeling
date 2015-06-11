% To run, save as booleanEx_synch.m or copy and paste in a MatLab m-file.
% Change the variables InitCond or t_stop if desired.  Save and run code 
% (under Debug menu).

function booleanEx_synch

% This function runs the binary rules in boolean form synchronously.  
% The only difference between this function and the "binary_sync" function
% is that the rules are specified using Boolean "truth" statements rather
% than equations.  Note that many more efficient ways of running boolean
% simulations are possible, but we have demonstrated the use using
% "if-then" statements because they are much more intuitive.
%
% Species Identifiers:
% (1)egf (2)hrg (3)egfr (4)raf (5)pi3k (6)erk (7)akt.
% To run, press F5 or choose Debug-> save and run.
% The following variables can be altered:
%
% InitCond: initial condition of each species given in species identifiers
% above
%
% t_stop: number of time steps to run simulation.


InitCond = [0 1 0 0 0 0 0];
t_stop = 10;

% Do not alter below

% set the initial condition
y(1,:) = InitCond;
step = 1;
% this loop calculates the values of each species during
while step < t_stop +1
    y(step+1,:) = Rules(y(step,:));
    step = step +1;
end

% plot timecourse of erk and akt
figure('Name','Boolean_Synch')
hold on
plot((0:1:t_stop),y(:,7),'c-.d','LineWidth',2,'MarkerSize',13)
plot((0:1:t_stop),y(:,6),'k.-','LineWidth',2,'MarkerSize',20)
set(gca,'fontsize',18)
xlabel('TimeStep','fontsize',30)
ylabel('Value of Species','fontsize',30)
legend({'Akt','Erk'},'fontsize',30)



%%Rules
function y = Rules(x)
y = x;

%egfr
if x(1) == 1 || x(2) == 1
    y(3) = 1;
else
    y(3) = 0;
end

%raf
if x(3) ==1 || x(7) == 1
    y(4) = 1;
else
    y(4) = 0;
end

%pi3k
if x(3) == 1 && ~x(6) == 1
    y(5) = 1;
else
    y(5) = 0;
end

%akt
if x(5) == 1
    y(7) = 1;
else
    y(7) = 0;
end

%erk
if x(4) == 1
    y(6) = 1;
else
    y(6) = 0;
end
