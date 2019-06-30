function [ Res ] = getFluxVariablities( Prob, lp_options )
%%getFluxVariablities computes the flux variabilities for the
%%specified problem
%
% Input:
%   Prob:        Problem structure
%   lp_options:  Solver option structure
%
% Output:        
%   Res:         Structure with fields
%     Values:      3D-matrix with the reactions in the first 
%                  dimension, the timepoints in the second dimension,
%                  and the lower/upper bounds as the first/second 
%                  element in the 3rd dimension 
%     Status:      3D-matrix with corresponding optimization status to
%                  Values matrix
%
  if ~exist('Prob', 'var'), error('cFBA:MissingParameter', 'Parameter ''Prob'' is missing.'); end
  if ~exist('lp_options', 'var'), error('cFBA:MissingParameter', 'Parameter ''lp_options'' is missing.'); end
  
  nR = size(Prob.model.S,2);
  nT = Prob.model.nT;
  nVars = size(Prob.S,2);
  
  if exist('log', 'dir'), log_file = fopen(sprintf('log/%s.log',mfilename),'w'); else log_file=0; end
   
  if log_file, fprintf(log_file, '%s\tComputation started.\n\n',time); end
  
  Res = struct();
  Res.Values = NaN(nR,nT,2);
  Res.Status = Res.Values;
  
  OK_STATUS = [1];
  
  % for each reaction
  for iR=1:nR
    if log_file, fprintf(log_file, '%s %19s (%2i): ', time, Prob.model.rxns{iR}, iR); end
    
    % for each interval
    for iT = 1:nT
      Objective = zeros(nVars,1);
      iVar = Prob.Vars.Fluxes(iT).Indices(iR);
      Objective(iVar) = 1;

      % maximize flux rate
      [~, xopt, ~, status] = lp_solve(-Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);

      if (status==1)
        % optimal status
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

      Res.Values(iR,iT,2) = xopt(iVar); 
      Res.Status(iR,iT,2) = status;

      if xopt(iVar)>Prob.lb(iVar)
        % Minimize flux rate
        [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);

        if (status==1)
          % optimal status
        elseif ismember(status, OK_STATUS)
          if log_file, fprintf(log_file, '[S:%i]', status); end
        else  
          % this probably happens for the same reason as mentioned above
          if ~exist('tmp','dir'), mkdir('tmp'); end; save 'tmp/ProblemObjective.mat' Objective; 
          error('cFBA:ProblemObjectiveDetected', 'Unexpected status (%i).', status); 
        end
      end

      Res.Values(iR,iT,1) = xopt(iVar);
      Res.Status(iR,iT,1) = status;
      
      if log_file
        % Show some progress in the log file
        if mod(iT,10)==0
          fprintf(log_file, '%i', mod(iT/10,10));
        else
          fprintf(log_file, '.');
        end
      end
    end
    if log_file, fprintf(log_file, '\n'); end 
  end
  
  if log_file
    fprintf(log_file, '\n%s\tComputation finished.\n',time);
    fclose(log_file);
  end
end

