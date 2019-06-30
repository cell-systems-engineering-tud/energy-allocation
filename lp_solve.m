function [fopt,xopt,slack,status] = lp_solve(c, S, b, lb, ub, lecon, gecon, options, x0)
% lp_solve solves the linear programming problem
% min c'*x
%
%Input: 
%   c         objective function
%   S         stoichiometric matrix
%   b         Right hand sight
%   lb        lower variable bounds
%   ub        upper variable bounds 
%   lecon     indices of <= u constraints 
%   gecon     indices of >= RHS constraints
%   options   solver options
%   x0        start solution
% 
%Output:
%   fopt       objective function value
%   xopt       solution
%   slack      slack variables
%   status     solver status
%
    
    fopt = NaN;
    xopt = NaN;
    slack = NaN;
    status = NaN;
    
    solver = options.solver;
    options = rmfield(options, 'solver');
    
    options = optimoptions('linprog','Preprocess','none','Display','none');
    options = optimoptions('linprog','Algorithm','dual-simplex','Display','none'); 
    if ~exist('x0','var')
        x0 = [];
    end
    
        
    switch solver
      case 'cplex'   
        if isempty(gecon)
          if isempty(lecon)
             [xopt,fopt,exitflag,output] = linprog(c,[],[],S,b,lb,ub,x0,options);
          else
            [xopt,fopt,exitflag,output] = linprog(c,S(lecon,:), b(lecon), S(~lecon,:), b(~lecon),lb,ub, x0,options);
          end
        elseif isempty(lecon) 
          [xopt,fopt,exitflag,output] = linprog(c,-S(gecon,:), -b(gecon), S(~gecon,:), b(~gecon),lb,ub, x0,options);
        else
          econ = ~or(lecon,gecon);
          [xopt,fopt,exitflag,output] = linprog(c,[S(lecon,:); -S(gecon,:)], [b(lecon); -b(gecon)], S(econ,:), b(econ),lb,ub, x0,options);
        end
        status = exitflag;
      otherwise
          error('Unknown solver.');
    end
end