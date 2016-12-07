% This example extends the basic MPC example by additional parametric
% inequalities. The inequalities are defined by a time-varying 2x2 matrix
% that is defined by 2 parameters (instead of 4 as one would expect).
%
%   R(k)*x <= R(k)*xmax
%
% where k is the simulation step and the rotation matrix is defined by 
%
%    R(k) = [ cos(k*w)  -sin(k*w)    = [ r1  -r2
%             sin(k*w)   cos(k*w) ]      r2   r1 ] 
%
% where k is the simulation step, and w a fixed number. Hence we have
%
%   r1 := cos(k*w)
%   r2 := sin(k*w)
%
% Overall, the following problem is solved at each time step:
% 
%  min   xN'*P*xN + sum_{i=1}^{N-1} xi'*Q*xi + ui'*R*ui
% xi,ui
%       s.t. x1 = x
%            x_i+1 = A*xi + B*ui  for i = 1...N-1
%            xmin <= xi <= xmax   for i = 1...N
%            umin <= ui <= umax   for i = 1...N
%            R(k)*xi <= R(k)*xmax for i = 1...N
%
% and P is solution of Ricatti eqn. from LQR problem

clear; clc;

%% MPC problem data

% system matrices
A = [1.1 1; 0 1];
B = [1; 0.5];
[nx,nu] = size(B);

% horizon
N = 10;

% cost matrics
Q = eye(2);
R = eye(1);
[~,P] = dlqr(A,B,Q,R);

% constraints
umin = -0.5;     umax = 0.5;
xmin = [-5; -5]; xmax = [5; 5];

%% Build MPC problem in Yalmip

% Initialize objective and constraints of the problem
cost = 0;
const = [];

% Cell arrays for x_0, x_1, ..., x_N and u_0, ..., u_N-1
x = num2cell(sdpvar(nx,N+1),1);
u = num2cell(sdpvar(nu,N),1);
sdpvar r1 r2

for i = 1:N        
    % cost
    if( i == N )
        cost = cost + 0.5*x{i+1}'*P*x{i+1} + 0.5*u{i}'*R*u{i};
    else
        cost = cost + 0.5*x{i+1}'*Q*x{i+1} + 0.5*u{i}'*R*u{i};
    end

    % model
    const = [const, x{i+1} == A*x{i} + B*u{i}];

    % bounds
    const = [const, umin <= u{i} <= umax];
    const = [const, xmin <= x{i+1} <= xmax]; 
    
    % rotation constraint R(k)*x <= R(k)*xmax
    const = [const, [r1, -r2; r2, r1]*x{i+1} <= [r1, -r2; r2, r1]*xmax];
     
end


%% Create controller object (generates code)
codeoptions = getOptions('FORCESsolver');
%codeoptions.printlevel = 0; % switch off FORCES printing (default is 2)

parameters = { x{1}, r1, r2 };
controller = optimizerFORCES(const, cost, codeoptions, parameters, u{1});


%% Simulate
x1 = [-4; 2];
kmax = 30;
X = zeros(nx,kmax+1); X(:,1) = x1;
U = zeros(nu,kmax);
problem.z1 = zeros(2*nx,1);
w = 1E-2;
for k = 1:kmax
    
    % Evaluate controller function for parameters
    r1 = cos(w*k); r2 = sin(w*k);
    [U(:,k),exitflag,info] = controller{ {X(:,k), r1, r2} };
    
    % Always check the exitflag in case something went wrong in the solver
    if( exitflag == 1 )
        fprintf('Time step %2d: FORCES took %2d iterations and %5.3f ', k,  info.it, info.solvetime*1000);
        fprintf('milliseconds to solve the problem.\n');
    else
        info
        error('Some problem in solver');
    end
    
    % State update
    X(:,k+1) = A*X(:,k) + B*U(:,k);
end


%% plot
figure(1); clf;
subplot(2,1,1); grid on; title('states'); hold on;
plot([1 kmax], [xmax xmax]', 'r--'); plot([1 kmax], [xmin xmin]', 'r--');
ylim(1.1*[min(xmin),max(xmax)]); stairs(1:kmax,X(:,1:kmax)');
subplot(2,1,2);  grid on; title('input'); hold on;
plot([1 kmax], [umax umax]', 'r--'); plot([1 kmax], [umin umin]', 'r--');
ylim(1.1*[min(umin),max(umax)]); stairs(1:kmax,U(:,1:kmax)');
