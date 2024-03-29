% Basic MPC example demonstrating the use of Yalmip to formulate MPC 
% problems and FORCESPRO to solve them very quickly.
%
% In this example, we will have the matrices A and B as parameters. This
% often occurs when a system is (re-)linearized around an operating point.
%
% Simple MPC - double integrator example for use with FORCESPRO
% 
%  min   xN'*P*xN + sum_{i=0}^{N-1} xi'*Q*xi + ui'*R*ui
% xi,ui
%       s.t. x0 = x(t)
%            x_i+1 = A*xi + B*ui  for i = 0...N-1
%            xmin <= xi <= xmax   for i = 1...N
%            umin <= ui <= umax   for i = 0...N-1
%
% Q, R, and P are suitable cost matrices and A and B are
% parameters to the problem.
%
% Note: due to 1-based indexing in Matlab, we use 1...N+1 instead of 0...N
%       as indices for state and input trajectory
%
% This file is part of the y2f project: http://github.com/embotech/y2f, 
% a project maintained by embotech under the MIT open-source license.
%
% (c) Gian Ulli and embotech AG, Zurich, Switzerland, 2013-2023.

clear; clc;

%% MPC problem data

% system dimensions
nx = 2; nu = 1;

% horizon
N = 10;

% cost matrices
Q = eye(2);
R = eye(1);
P = 10*Q; % we don't know A and B yet, so we cannot compute the LQR cost and use this instead

% constraints
umin = -0.5;     umax = 0.5;
xmin = [-5; -5]; xmax = [5; 5];

%% Build MPC problem in Yalmip

% Define variables
x0 = sdpvar(nx,1); % initial state
X = sdpvar(nx,N+1,'full'); % state trajectory: x1,...,xN (columns of X)
U = sdpvar(nu,N,'full'); % input trajectory: u0,...,u_{N-1} (columns of U)
A = sdpvar(nx,nx,'full'); % system matrix - parameter
B = sdpvar(nx,nu,'full'); % input matrix - parameter

% Initialize objective and constraints of the problem
cost = 0;
const = x0 == X(:,1);

% Assemble MPC formulation
for i = 1:N        
    % cost
    if( i < N )
        cost = cost + 0.5*X(:,i+1)'*Q*X(:,i+1) + 0.5*U(:,i)'*R*U(:,i);
    else
        cost = cost + 0.5*X(:,N+1)'*P*X(:,N+1) + 0.5*U(:,N)'*R*U(:,N);
    end
    
    % model
    const = [const, X(:,i+1) == A*X(:,i) + B*U(:,i)];

    % bounds
    const = [const, umin <= U(:,i) <= umax];
    const = [const, xmin <= X(:,i+1) <= xmax];
end

%% Create controller object (generates code)
% for a complete list of codeoptions, see 
% https://forces.embotech.com/Documentation/solver_options/index.html
codeoptions = getOptions('parametricDynamics_solver'); % give solver a name
parameters = { x0, A, B };
parameterNames = { 'xinit', 'Amatrix', 'Bmatrix' };
controller = optimizerFORCES(const, cost, codeoptions, parameters, U(:,1), parameterNames, {'u0'} );


%% Simulate
x1 = [-4; 2];
kmax = 30;
X = zeros(nx,kmax+1); X(:,1) = x1;
U = zeros(nu,kmax);
problem.z1 = zeros(2*nx,1);

for k = 1:kmax
    
    % Set system matrices
    A = [1.1 1; 0 1];
    B = [1; 0.5];
    
    % Evaluate controller function for parameters
    [U(:,k),exitflag,info] = controller{ X(:,k), A, B };
    
    % Always check the exitflag in case something went wrong in the solver
    if( exitflag == 1 )
        fprintf('Time step %2d: FORCESPRO took %2d iterations and %5.3f ', k,  info.it, info.solvetime*1000);
        fprintf('milliseconds to solve the problem.\n');
    else
        disp(info);
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
