function [Prob] = setMaximalCompoundFoldChange(Prob, MIN_PRECISION, lp_options)
%%setMaximalCompoundFoldChange computes the maximal compund fold change
%%(=\alpha)
%
% Input:
%   Prob:           LP problem structure
%   MIN_PRECISION:  Precision of the fold change which have to be achieved
%   lp_options:     solver option structure
%
% Output:
%   Prob:           Updated problem with maximal \alpha
%
  if exist('log','dir'), log_file = fopen(sprintf('log/%s.log',mfilename),'w'); else log_file=0; end
  if log_file>0, fprintf(log_file, '%s\tComputation started.\n',time); end
  
  PO = getSetOfProblemObjectives(Prob);
  nPO = size(PO,2);
  
  nVars = size(Prob.S,2);
  Objective = zeros(nVars,1);
  
  OK_STATUS = [1 ];
  EXPECTED_STATUS = [1 0 -2 -3 -4 -5 -7];
  
  % initiate Problem with non growth fold change
  Prob = setCompoundFoldChange(Prob,1);
  [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options);

  if ~ismember(status, OK_STATUS), error('Problem has no non-growth solution.'); end
  
  Feasible_FoldChanges = 1;
  Solutions = xopt;
  
  % double fold change until there is no feasible solution to the problem
  while ismember(status, OK_STATUS)
    Prob = setCompoundFoldChange(Prob,Feasible_FoldChanges(end)*2);
    if log_file>0, fprintf(log_file, '%s\tLevel 1 check FoldChange: %12.11g ... ', time, Prob.CompoundFoldChange); end
    [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options);
    
    if ismember(status, OK_STATUS)
      if log_file>0, fprintf(log_file, 'ok (%i).\n', status); end
      Feasible_FoldChanges = [Feasible_FoldChanges Prob.CompoundFoldChange];
      Solutions = [Solutions xopt];
    else
      if log_file>0, fprintf(log_file, 'failed (%i).\n', status); end
    end
  end
  
  Delta_FoldChange = Feasible_FoldChanges(end)/2;
  LastCheck2LevelOKIndex = 0;
  
  while Delta_FoldChange>MIN_PRECISION
    
    % Level 1 check
    % increase the precision step wise
    while Delta_FoldChange>MIN_PRECISION
      Prob = setCompoundFoldChange(Prob,Feasible_FoldChanges(end)+Delta_FoldChange);
      if log_file>0, fprintf(log_file, '%s\tLevel 1 check FoldChange: %12.11g, PRECISION: %g ... ', time, Prob.CompoundFoldChange, Delta_FoldChange); end
      [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options);

      if ismember(status, OK_STATUS)
        if log_file>0, fprintf(log_file, 'ok (%i).\n', status); end
        Feasible_FoldChanges = [Feasible_FoldChanges Prob.CompoundFoldChange];
        Solutions = [Solutions xopt];
      else
        if log_file>0, fprintf(log_file, 'failed (%i).\n', status); end
      end
      Delta_FoldChange = Delta_FoldChange/2;
    end
  
    % Perform level 2 check, backtrack if necessary
    % this check is necessary, because the level 1 check might yield
    % slightly infeasible FoldChanges (due to numerical artefacts of the solver)
    % by testing these feasible FoldChanges of the level 1 check with
    % problematic objectives, these slightly infeasible foldchanges are
    % hopefully discarded
    if nPO>0
      if numel(Feasible_FoldChanges)>LastCheck2LevelOKIndex
        for k=length(Feasible_FoldChanges):-1:(LastCheck2LevelOKIndex+1)
          Prob = setCompoundFoldChange(Prob,Feasible_FoldChanges(k));
          if log_file>0, fprintf(log_file, '%s\tLevel 2 check FoldChange: %12.11g [', time, Prob.CompoundFoldChange); end
          for iPO = 1:nPO    
            if log_file>0
              if iPO>1, fprintf(log_file, ', '); end
              fprintf(log_file, 'opt PO%i ... ', iPO);
            end

            [~, ~, ~, status] = lp_solve(PO(:,iPO), Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Solutions(:,k));
            if ~ismember(status,EXPECTED_STATUS), 
                error('Unexpected status (%i).', status); end
            if ~ismember(status,OK_STATUS), 
              if log_file>0, fprintf(log_file, 'failed (%i)', status); end
              break;
            else
              if log_file>0, fprintf(log_file, 'ok (%i)', status); end
            end
          end
          log_file = 1;
          fprintf(log_file, ']\n');
          log_file = 0;
          if ismember(status,OK_STATUS), break; end
        end

        if ismember(status,OK_STATUS) 
          if k==numel(Feasible_FoldChanges), 
            break;
          else
            Delta_FoldChange = (Feasible_FoldChanges(k+1)-Feasible_FoldChanges(k))/2;
            Feasible_FoldChanges = Feasible_FoldChanges(1:k);
            Solutions = Solutions(:,1:k);
            LastCheck2LevelOKIndex = k;
          end
        elseif LastCheck2LevelOKIndex==0
          error('Level 2 check failed for non growth setting.');
        else
          Delta_FoldChange = (Feasible_FoldChanges(LastCheck2LevelOKIndex+1)-Feasible_FoldChanges(LastCheck2LevelOKIndex))/2;
          Feasible_FoldChanges = Feasible_FoldChanges(1:LastCheck2LevelOKIndex);
          Solutions = Solutions(:,1:LastCheck2LevelOKIndex);
        end
      end
    else
      % no level 2 checkes given
      if log_file>0, fprintf(log_file, '%s\tLevel 2 check not performed (no ProblemObjective given).\n', time); end
    end
    
  end
  
  Prob = setCompoundFoldChange(Prob,Feasible_FoldChanges(end));  
  Prob.x0 = Solutions(:,end);
  if log_file>0, 
    fprintf(log_file, '%s\tResulting compound fold change: %12.11g\n', time, Prob.CompoundFoldChange); 
    fprintf(log_file, '\n%s\tComputation finished.\n',time); 
    fclose(log_file); 
  end
end
%
function [ PO ] = getSetOfProblemObjectives( Prob )
%%getSetOfProblemObjectives loads the set of 'Problem'-Objectives for the
%%specified problem
%
% Input: 
%   Prob:   problem structure
%
% Output:
%     PO:   Matrix in which each columns is a 'Problem'-Objective
%
  PROBLEM_OBJECTIVES_FILE = sprintf('tmp/ProblemObjectives_%s.mat', Prob.model.ID);
  nVar = size(Prob.S,2);

  if exist(PROBLEM_OBJECTIVES_FILE, 'file')
    load(PROBLEM_OBJECTIVES_FILE,'PO');
    if size(PO,1)~=nVar, error('cFBA:InvalidInput', 'Problem objectives do not match problem size.'); end
  else
    PO = zeros(nVar,0);
  end
end


