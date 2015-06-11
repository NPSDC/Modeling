% To run, save as multiEx_asynch.m or copy and paste in a MatLab m-file.
% Change the variables InitCond or t_stop if desired.  Save and run code 
% (under Debug menu).

function multiEx_asynch

% This function runs the multi-level model asynchronously
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


InitCond = [1 1 0 0 0 0 0];
t_stop = 20;

% Do not alter below

% Enumerate the possible orders, to be chosen randomly if updatingScheme is
% asynch
possOrders = perms(1:1:5);

% set the initial condition
y(1,:) = InitCond;
step = 1;
% this loop calculates the values of each species during 
while step < t_stop +1    
    % a random order is chosen
    order  = possOrders(randi(size(possOrders,  1)),:);
    % This function calls the rules specified below
    y(step+1,:) = Rules(y(step,:),order);
    fprintf('%d',step);
    y(step + 1,:)
    step = step +1;
end

% plot timecourse of erk and akt
figure('Name','Multi_asynch')
hold on
plot((0:1:t_stop),y(:,6),'k-o','LineWidth',2,'MarkerSize',8)
plot((0:1:t_stop),y(:,7),'r:d','LineWidth',2,'MarkerSize',5)
set(gca,'fontsize',18)
xlabel('TimeStep','fontsize',30)
ylabel('Value of Species','fontsize',30)
legend({'Erk','Akt'},'fontsize',30)



%%Rules
function y = Rules(x,order)
% Because of the random ordering, each rule must now be evaluated one at a time. 
y = x;

for i = 1:length(order) %length(order) is 5
    y = EachRule(y,order(i));
end

function y = EachRule(x,Choose)

y = x;

if Choose == 1
    %egfr
    if x(1) + 0.5* x(2) >= 1
        y(3) = 2;
    elseif x(1) +0.5*x(2) > 0
        y(3) = 1;
    else
        y(3) = 0;
    end
    
elseif Choose == 2
    %raf
    if x(3) + 2* x(7) > 1
        y(4) = 1;
    else
        y(4) = 0;
    end
    
elseif Choose == 3
    %pi3k
    if x(3) - x(6) > 0
        y(5) = 1;
    else
        y(5) = 0;
    end
    
elseif Choose == 5
    %akt
    if x(5) > 0
        y(7) = 1;
    else
        y(7) = 0;
    end
    
elseif Choose == 4
    %erk
    if x(4) > 0
        y(6) = 1;
    else
        y(6) = 0;
    end
end