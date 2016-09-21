# Y2F - YALMIP to FORCES Pro Interface

This project provides a simple MATLAB interface that connects [YALMIP](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.WhatIsYALMIP)
with [FORCES Pro](https://www.embotech.com/FORCES-Pro). It combines YALMIP's intuitiveness with the high efficiency of FORCES Pro for rapid development.

## Installation

Simply download the code to the desired location and add the `Y2F` folder to your [MATLAB search path](http://ch.mathworks.com/help/matlab/ref/addpath.html). 

The Y2F interface requires a working YALMIP installation. See [http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Installation](http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Installation) for instructions on how to install YALMIP.

The code has been tested with YALMIP release 20150919. Older versions might work but have not been tested.

## Example Usage

Consider the following linear MPC problem with lower and upper bounds on state and inputs, and a terminal cost term:

![\begin{aligned}\text{minimize} \quad & x_N^T P x_N + \sum_{i=0}^{N-1} x_i^T Q x_i + u_i^T R u_i \\ \text{s.t.} \quad & x_0 = x(t) \\& x_{i+1} = Ax_i + Bu_i \\& \underline{x} \leq x_i \leq \bar{x} \\& \underline{u} \leq u_i \leq \bar{u}\end{aligned}](http://www.sciweavers.org/tex2img.php?eq=%5Cbegin%7Baligned%7D%0A%5Ctext%7Bminimize%7D%20%5Cquad%20%26%20x_N%5ET%20P%20x_N%20%2B%20%5Csum_%7Bi%3D0%7D%5E%7BN-1%7D%20x_i%5ET%20Q%20x_i%20%2B%20u_i%5ET%20R%20u_i%20%5C%5C%20%0A%5Ctext%7Bs.t.%7D%20%5Cquad%20%26%20x_0%20%3D%20x%28t%29%20%5C%5C%0A%26%20x_%7Bi%2B1%7D%20%3D%20Ax_i%20%2B%20Bu_i%20%5C%5C%0A%26%20%5Cunderline%7Bx%7D%20%5Cleq%20x_i%20%5Cleq%20%5Cbar%7Bx%7D%20%5C%5C%0A%26%20%5Cunderline%7Bu%7D%20%5Cleq%20u_i%20%5Cleq%20%5Cbar%7Bu%7D%0A%5Cend%7Baligned%7D&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0)

This problem is parametric in the initial state x(t), and the first input u0 is typically applied to the system after a solution has been obtained. The following code generates a solver that returns u0, which can then be applied to the system:

```
% Objective and constraints of the problem
cost = 0;
const = [];

% Cell arrays for x_0, x_1, ..., x_N and u_0, ..., u_N-1
x = num2cell(sdpvar(nx,N+1),1);
u = num2cell(sdpvar(nu,N),1);

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
end

controller = optimizerFORCES(const, cost, codeoptions, x{1}, u{1});
```

The generated solver can then be called using curly braces:

```
u0 = controller{x0};
```

## Limitations

The Y2F interface only supports convex quadratically constrained programs (QCQPs). FORCES' NLP solver is currently not supported.

## License

The code is licensed under the MIT License. For more information see LICENSE file.
