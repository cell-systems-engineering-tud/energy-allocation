function [ Res ] = getCompoundStartAmountVariabilities( Prob, lp_options )
%%getCompoundStartAmountVariabilities computes the variabilities of the compound start amounts
%
% Input:
%   Prob:        Problem structure
%   lp_options:  Solver option structure
%
% Output:        
%   Res:         Structure with fields
%     Values:      Matrix with first column as lower bounds and second column
%                  as upper bounds
%     Status:      Corresponding optimization status
%
  if ~exist('Prob', 'var'), error('cFBA:MissingParameter', 'Parameter ''Prob'' is missing.'); end
  if ~exist('lp_options', 'var'), error('cFBA:MissingParameter', 'Parameter ''lp_options'' is missing.'); end  

  nMLi = length(Prob.model.ImbalancedMets);
  Res = struct('Values', {NaN(nMLi,2)},'Status', NaN(nMLi,2));
  
  if exist('log', 'dir'), log_file = fopen(sprintf('log/%s.log',mfilename),'w'); else log_file=0; end
  
  if log_file, fprintf(log_file, '%s\tComputation started.\n\n',time); end
    
  OK_STATUS = [1 ];
  
  nVars = size(Prob.S,2);
 
  for iM=1:nMLi
    Objective = zeros(nVars,1);
    Objective(Prob.Vars.ImbMets.StartAmounts(iM)) = 1;
    if log_file, fprintf(log_file, '%s %17s (%2i): ',time,Prob.model.mets{Prob.model.ImbalancedMets(iM)},iM); end
    
    % maximize compound amount
    [~, xopt, ~, status] = lp_solve(-Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);
    
    if (status==1)
      % optimal status
      if log_file, fprintf(log_file, '.'); end
    elseif ismember(status, OK_STATUS)
      if log_file, fprintf(log_file, '[S:%i]', status); end
    else  
      % this probably happens, if our \alpha (compound fold change) is
      % slightly infeasible
      % to compute a better \alpha the so called 'Problem'-Objective is stored
      % A calling instance can catch this error and restart the \alpha
      % computation
      Objective=-Objective; if ~exist('tmp','dir'), mkdir('tmp'); end; save 'tmp/ProblemObjective.mat' Objective; 
      error('cFBA:ProblemObjectiveDetected', 'Unexpected status (%i).', status); 
    end
    Res.Values(iM,2) = xopt(Prob.Vars.ImbMets.StartAmounts(iM));
    Res.Status(iM,2) = status;
    
    % minimize compound amount
    if Res.Values(iM,2)==0
      % upper bound = zero implies the lower bound 
      Res.Values(iM,1) = Res.Values(iM,2);
      Res.Status(iM,1) = Res.Status(iM,2);
      if log_file, fprintf(log_file, '.'); end
    else
      [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);

      if (status==1)
        % optimal status
        if log_file, fprintf(log_file, '.'); end
      elseif ismember(status, OK_STATUS)
        if log_file, fprintf(log_file, '[S:%i]', status); end
      else 
        % this probably happens for the same reason as mentioned above
        if ~exist('tmp','dir'), mkdir('tmp'); end; save 'tmp/ProblemObjective.mat' Objective; 
        error('cFBA:ProblemObjectiveDetected', 'Unexpected status (%i).', status); 
      end
      Res.Values(iM,1) = xopt(Prob.Vars.ImbMets.StartAmounts(iM));
      Res.Status(iM,1) = status;
    end
    if log_file, fprintf(log_file, '\n'); end
  end
  
  if log_file
    fprintf(log_file, '\n%s\tComputation finished.\n',time);
    fclose(log_file);
  end
end

