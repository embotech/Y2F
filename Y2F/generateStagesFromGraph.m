function [stages, params, standardParamValues, forcesParamMap] = generateStagesFromGraph( G, H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,qcqpParams,yalmipParamMap )
%GENERATESTAGESFROMGRAPH Generates stages from path graph containing all variables of the problem
% qcqpParams will be modified to contain a link to the params they are linked to
%   Input:
%       yalmipParamMap:         Map that maps YALMIP parameters to value matrices
%                                   1st row: index of matrix with values,
%                                   2nd row: index of element inside matrix
%       
%   Output:
%       stages:                 FORCES stages
%       params:                 FORCES parameters
%       standardParamValues:    "standard" values for FORCES parameters. 
%                               Not every element of the matrix has be 
%                               covered by a QCQP parameter. Param values
%                               get added to this value.
%       forcesParamMap:         Map that maps FORCES parameters to parameter value matrices
%       	.(param_name)(1,.):   index of element that's affected
%           .(param_name)(2,.):   factor by which value gets multiplied
%           .(param_name)(3,.):   matrix that contains value
%           .(param_name)(4,.):   points to element inside value matrix

if ~checkIfGraphIsPathGraph(G)
    error('Path is not a path graph. Other structures are not supported.')
end

% Find all variable names
last_idx = G.vertices{1};
all_idx = sort(cell2mat(G.vertices));

% Prepare stages
stages = MultistageProblem(G.n);
params = struct('name',{},'maps2stage',{},'maps2data',{},'type',{},...
    'maps2mat',{},'structure',{},'variables',{}); % empty param struct
standardParamValues = struct;
p = 1;
if nargout >= 4
    forcesParamMap = struct;
end

% Construct matrices where parametric elements are == 1
Aineq_temp = Aineq;
for i=1:numel(qcqpParams.Aineq)
    Aineq_temp(qcqpParams.Aineq(i).maps2index) = 1;
end
Aeq_temp = Aeq;
for i=1:numel(qcqpParams.Aeq)
    Aeq_temp(qcqpParams.Aeq(i).maps2index) = 1;
end
Q_temp = Q;
for i=1:numel(qcqpParams.Q)
    Q_temp{qcqpParams.Q(i).maps2mat}(qcqpParams.Q(i).maps2index) = 1;
end
l_temp = l;
for i=1:numel(qcqpParams.l)
    l_temp(qcqpParams.l(i).maps2index,qcqpParams.l(i).maps2mat) = 1;
end

% Go through all stages (and create them)
for i=1:G.n
    idx = G.vertices{i}; % sorting just to make sure code below works
    stages(i).dims.n = length(idx); % length of stage variable zi 
    stages(i).cost.H = H(idx,idx); % Hessian
    stages(i).cost.f = f(idx); % linear term
    
    % Is H a parameter? (TODO: H is a special case! Diagonal?!)
    relevantParams = findRelevantParams(idx,idx,size(H),qcqpParams.H);
    if ~isempty(relevantParams) % Yes!
        % Check if H is diagonal
        if all(stages(i).cost.H(~eye(length(idx))) == 0) && ... % nondiagonal entries are zero
                all(ismember([qcqpParams.H(relevantParams).maps2index], sub2ind(size(H),idx,idx))) % all parameters affect diagonal entries
            % We have a diagonal H --> create a diagonal parameter
            createDiagonalCostParameter(i, relevantParams, idx);
        else % H is not diagonal --> do the standard thing
            createParameter('H', i, relevantParams, size(H), idx, idx, 'cost.H');
        end
    end
    % Is f a parameter?
    relevantParams = findRelevantParams(idx,1,size(f),qcqpParams.f);
    if ~isempty(relevantParams)
        createParameter('f', i, relevantParams, size(f), idx, 1, 'cost.f');
    end
    
    % Select relevant equalities and set the constraints
    eq_idx = find(sum(Aeq_temp(:,idx)~=0, 2))';
    eq_idx = eq_idx(sum(Aeq_temp(eq_idx,setdiff(all_idx,[last_idx idx]))~=0,2)==0);
    eq_idx = sort(eq_idx);
    if i == 1
        d1IsSet = ~isempty(Aeq_temp(eq_idx,idx));
    end
    if d1IsSet
        if i > 1
            stages(i-1).eq.C = Aeq(eq_idx,last_idx);
        end
        stages(i).dims.r = length(eq_idx); % number of equality constraints
        stages(i).eq.c = beq(eq_idx);
        stages(i).eq.D = Aeq(eq_idx,idx);
    else
        if i > 1
            stages(i-1).dims.r = length(eq_idx); % number of equality constraints
            stages(i-1).eq.C = Aeq(eq_idx,last_idx);
            stages(i-1).eq.c = beq(eq_idx);
        end
        stages(i).eq.D = Aeq(eq_idx,idx);
    end

    % Is C a parameter?
    relevantParams = findRelevantParams(eq_idx,last_idx,size(Aeq),qcqpParams.Aeq);
    if ~isempty(relevantParams) && i > 1
        createParameter('Aeq', i-1, relevantParams, size(Aeq), eq_idx, last_idx, 'eq.C');
    end
    % Is D a parameter?
    relevantParams = findRelevantParams(eq_idx,idx,size(Aeq),qcqpParams.Aeq);
    if ~isempty(relevantParams)
        createParameter('Aeq', i, relevantParams, size(Aeq), eq_idx, idx, 'eq.D');
    end
    % Is c a parameter?
    relevantParams = findRelevantParams(eq_idx,1,size(beq),qcqpParams.beq);
    if ~isempty(relevantParams)
        if d1IsSet
            createParameter('beq', i, relevantParams, size(beq), eq_idx, 1, 'eq.c');
        else
            createParameter('beq', i-1, relevantParams, size(beq), eq_idx, 1, 'eq.c');
        end
    end
    
    % Lower bounds
    temp_lb = lb(idx);
    stages(i).dims.l = sum(temp_lb ~= -Inf); % number of lower bounds 
    stages(i).ineq.b.lbidx = find(temp_lb ~= -Inf)'; % index vector for lower bounds
    stages(i).ineq.b.lb = temp_lb(temp_lb ~= -Inf);    % lower bounds
    
    % Is lb a parameter?
    relevantParams = findRelevantParams(idx,1,size(lb),qcqpParams.lb);
    if ~isempty(relevantParams)
        createParameter('lb', i, relevantParams, size(lb), idx(stages(i).ineq.b.lbidx), 1, 'ineq.b.lb');
    end
    
    % Upper bounds
    temp_ub = ub(idx);
    stages(i).dims.u = sum(temp_ub ~= Inf); % number of upper bounds 
    stages(i).ineq.b.ubidx = find(temp_ub ~= Inf)'; % index vector for upper bounds
    stages(i).ineq.b.ub = temp_ub(temp_ub ~= Inf);    % upper bounds
    
    % Is ub a parameter?
    relevantParams = findRelevantParams(idx,1,size(ub),qcqpParams.ub);
    if ~isempty(relevantParams)
        createParameter('ub', i, relevantParams, size(ub), idx(stages(i).ineq.b.ubidx), 1, 'ineq.b.ub');
    end
    
    % Linear inequalities
    ineq_idx = find(sum(Aineq_temp(:,idx)~=0, 2))';
    stages(i).dims.p = length(ineq_idx); % number of polytopic constraints
    if ~isempty(ineq_idx)
        stages(i).ineq.p.A = Aineq(ineq_idx, idx); % Jacobian of linear inequality 
        stages(i).ineq.p.b = bineq(ineq_idx); % RHS of linear inequality

        % Is p.A a parameter?
        relevantParams = findRelevantParams(ineq_idx,idx,size(Aineq),qcqpParams.Aineq);
        if ~isempty(relevantParams)
            createParameter('Aineq', i, relevantParams, size(Aineq), ineq_idx, idx, 'ineq.p.A');
        end
        % Is p.b a parameter?
        relevantParams = findRelevantParams(ineq_idx,1,size(bineq),qcqpParams.bineq);
        if ~isempty(relevantParams)
            createParameter('bineq', i, relevantParams, size(bineq), ineq_idx, 1, 'ineq.p.b');
        end
    end
    
    
    % Quadratic constraints
    stages(i).dims.q = 0; % number of quadratic constraints
    if ~isempty(Q)
        stages(i).ineq.q.idx = {}; % index vectors 
        stages(i).ineq.q.Q = {}; % Hessians
        stages(i).ineq.q.l = {}; % linear terms
        stages(i).ineq.q.r = []; % RHSs 

        for k=1:numel(Q) % go through each quad. inequality
            subQ = Q_temp{k}(idx,idx);
            subL = l_temp(idx,k);
            if ~isempty([find(subQ,1) find(subL,1)]) % Inequality is relevant
                quad_idx = [];
                for s=1:size(subQ,1)
                    quad_idx = [quad_idx find(subQ(s,:)~=0)];
                    if subL(s)~=0
                        quad_idx(end+1) = s;
                    end
                end
                quad_idx = sort(unique(quad_idx));
                
                stages(i).ineq.q.idx{end+1} = quad_idx;
                stages(i).ineq.q.Q{end+1} = full(Q{k}(idx(quad_idx),idx(quad_idx)));
                stages(i).ineq.q.l{end+1} = full(l(idx(quad_idx),k)); % linear terms
                stages(i).ineq.q.r(end+1) = r(k); % RHSs 
                stages(i).dims.q = stages(i).dims.q + 1;
                    
                orig_idx = idx(quad_idx);
                
                % Is Q a param? (Q,l,r are special cases because of k)
                relevantParams = findRelevantParams(idx,idx,size(Q{k}),qcqpParams.Q,k);
                if ~isempty(relevantParams)
                    params(p) = newParam(sprintf('p_%u',p), i, 'ineq.q.Q',stages(i).dims.q);
                    standardParamValues.(sprintf('p_%u',p)) = stages(i).ineq.q.Q{stages(i).dims.q};
                    stages(i).ineq.q.Q{stages(i).dims.q} = [];
                    
                    forcesParamMap.(sprintf('p_%u',p)) = zeros(4,0);
                    paramSize = size(standardParamValues.(sprintf('p_%u',p)));
                    for j=relevantParams
                        for col=1:length(quad_idx)
                            for row=1:length(quad_idx)
                                value_idx = find(sub2ind(size(Q{k}),orig_idx(row),orig_idx(col)) == ...
                                            qcqpParams.Q(j).maps2index);
                                if length(value_idx) == 1
                                    forcesParamMap.(sprintf('p_%u',p))(:,end+1) = [sub2ind(paramSize,row,col);...
                                                                                   qcqpParams.Q(j).factor;...
                                                                                   yalmipParamMap(:,qcqpParams.Q(j).maps2origparam)];
                                elseif length(value_idx) > 1
                                    error('Mistake in the new stages formulation')
                                end
                            end
                        end
                    end
                    p = p + 1;
                end
                % Is l a param?
                relevantParams = findRelevantParams(idx,1,size(l(:,k)),qcqpParams.l,k);
                if ~isempty(relevantParams)
                    params(p) = newParam(sprintf('p_%u',p), i, 'ineq.q.l',stages(i).dims.q);
                    standardParamValues.(sprintf('p_%u',p)) = stages(i).ineq.q.l{stages(i).dims.q};
                    stages(i).ineq.q.l{stages(i).dims.q} = [];
                    
                    forcesParamMap.(sprintf('p_%u',p)) = zeros(4,0);
                    for j=relevantParams
                        for row=1:length(quad_idx)
                            value_idx = find(orig_idx(row) == qcqpParams.l(j).maps2index);
                            if length(value_idx) == 1
                                forcesParamMap.(sprintf('p_%u',p))(:,end+1) = [row;...
                                                                               qcqpParams.l(j).factor;...
                                                                               yalmipParamMap(:,qcqpParams.l(j).maps2origparam)];
                            elseif length(value_idx) > 1
                                error('Mistake in the new stages formulation')
                            end
                        end
                    end
                    p = p + 1;
                end
                % Is r a param?
                relevantParams = findRelevantParams(k,1,size(r),qcqpParams.r);
                if ~isempty(relevantParams)
                    assert(length(relevantParams)==1);
                    params(p) = newParam(sprintf('p_%u',p), i, 'ineq.q.r',stages(i).dims.q);
                    standardParamValues.(sprintf('p_%u',p)) = stages(i).ineq.q.r(stages(i).dims.q);
                    stages(i).ineq.q.r(stages(i).dims.q) = 0;
                    
                    forcesParamMap.(sprintf('p_%u',p)) = zeros(4,0);
                    for j=relevantParams
                        forcesParamMap.(sprintf('p_%u',p))(:,end+1) = [1;...
                                                                       qcqpParams.r(j).factor;...
                                                                       yalmipParamMap(:,qcqpParams.r(j).maps2origparam)];
                    end
                    p = p + 1;
                end
            end
        end
    end
    
    % Assign binary variables, if there are any
    bidx = find(ismember(idx,qcqpParams.bidx));
    if ~isempty(bidx)
        stages(i).bidx = bidx;
    end
    
    % Assign boundaries on binary variables if needed
    if isfield(stages(i),'bidx')
        for j=stages(i).bidx
            if ~ismember(j,stages(i).ineq.b.lbidx)
                stages(i).ineq.b.lbidx(end+1) = j;
                stages(i).ineq.b.lb(end+1) = 0;
                stages(i).dims.l = stages(i).dims.l + 1;
            end
            if ~ismember(j,stages(i).ineq.b.ubidx)
                stages(i).ineq.b.ubidx(end+1) = j;
                stages(i).ineq.b.ub(end+1) = 1;
                stages(i).dims.u = stages(i).dims.u + 1;
            end
        end
    end
    
    last_idx = idx;
end

    function createParameter(matrix, stage, relevantParams, matrixSize, row_idx, col_idx, name)
    % Helper function to create FORCES parameters

        params(p) = newParam(sprintf('p_%u',p), stage, name);
        var = regexp(name,'[.]','split');  % split name
        if length(var) == 2
            standardParamValues.(sprintf('p_%u',p)) = stages(stage).(genvarname(var{1})).(genvarname(var{2}));
            stages(stage).(genvarname(var{1})).(genvarname(var{2})) = [];
        elseif length(var) == 3
            standardParamValues.(sprintf('p_%u',p)) = stages(stage).(genvarname(var{1})).(genvarname(var{2})).(genvarname(var{3}));
            stages(stage).(genvarname(var{1})).(genvarname(var{2})).(genvarname(var{3})) = [];
        elseif length(var) == 4
            standardParamValues.(sprintf('p_%u',p)) = stages(stage).(genvarname(var{4}));
            stages(stage).(genvarname(var{1})).(genvarname(var{2})).(genvarname(var{3})).(genvarname(var{4})) = [];
        else
            error('This case exists?');
        end

        forcesParamMap.(sprintf('p_%u',p)) = zeros(4,0);
        paramSize = size(standardParamValues.(sprintf('p_%u',p)));
        for j=relevantParams
            for col=1:length(col_idx)
                for row=1:length(row_idx)
                    value_idx = find(sub2ind(matrixSize,row_idx(row),col_idx(col)) == ...
                            qcqpParams.(matrix)(j).maps2index);
                    if length(value_idx) == 1
                        forcesParamMap.(sprintf('p_%u',p))(:,end+1) = [sub2ind(paramSize,row,col);...
                                                                       qcqpParams.(matrix)(j).factor;...
                                                                       yalmipParamMap(:,qcqpParams.(matrix)(j).maps2origparam)];
                    elseif length(value_idx) > 1
                        error('Mistake in the new stages formulation')
                    end
                end
            end
        end
        p = p + 1;
    end

    function createDiagonalCostParameter(stage, relevantParams, element_idx)
    % Helper function to create a FORCES parameters for a diagonal cost
        params(p) = newParam(sprintf('p_%u',p), stage, 'cost.H', 'diag');
        standardParamValues.(sprintf('p_%u',p)) = stages(stage).cost.H(logical(eye(length(element_idx)))); % only use diagonal
        stages(stage).cost.H = [];

        forcesParamMap.(sprintf('p_%u',p)) = zeros(4,0);
        paramSize = size(standardParamValues.(sprintf('p_%u',p)));
        for j=relevantParams
            for element=1:length(element_idx)
                value_idx = find(sub2ind(size(H),element_idx(element),element_idx(element)) == ...
                        qcqpParams.H(j).maps2index);
                if length(value_idx) == 1
                    forcesParamMap.(sprintf('p_%u',p))(:,end+1) = [element;...
                                                                   qcqpParams.H(j).factor;...
                                                                   yalmipParamMap(:,qcqpParams.H(j).maps2origparam)];
                elseif length(value_idx) > 1
                    error('Mistake in the new stages formulation')
                end
            end
        end
        p = p + 1;
    end
            

end



