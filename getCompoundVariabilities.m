function [ Res ] = getCompoundVariabilities( Prob, lp_options )
%%getCompoundVariabilities computes the amount variabilities of the
%%imbalanced compounds
%
% Input:
%   Prob:        Problem structure
%   lp_options:  Solver option structure
%
% Output:        
%   Res:         Structure with fields
%     Values:      3D-matrix with the compounds in the first 
%                  dimension, the timepoints in the second dimension,
%                  and the lower/upper bounds as the first/second 
%                  element in the 3rd dimension 
%     Status:      3D-matrix with corresponding optimization status to
%                  Values matrix
%
  if ~exist('Prob', 'var'), error('cFBA:MissingParameter', 'Parameter ''Prob'' is missing.'); end
  if ~exist('lp_options', 'var'), error('cFBA:MissingParameter', 'Parameter ''lp_options'' is missing.'); end
  
  nT = Prob.model.nT;
  nMLi = length(Prob.model.ImbalancedMets);
  nVars = size(Prob.S,2);
  
  
  if exist('log', 'dir'), log_file = fopen(sprintf('log/%s.log',mfilename),'w'); else log_file=0; end
   
  if log_file, fprintf(log_file, '%s\tComputation started.\n\n',time); end
  
  Res = struct();
  Res.Values = NaN(nMLi,nT,2);
  Res.Status = Res.Values;
  
  OK_STATUS = [1];
  
  % for each imbalanced compound
  for iM=1:nMLi
    metName = Prob.model.mets{Prob.model.ImbalancedMets(iM)};
    if log_file, fprintf(log_file, '%s %17s (%2i): ', time, metName, iM); end
    
    % for each time step but the start
    for iT = 1:nT
         
      Objective = zeros(nVars,1);
      Objective(1:Prob.Vars.ImbMets.Tind(iT)) = Prob.Vars.ImbMets.S(iM,1:Prob.Vars.ImbMets.Tind(iT));


      % Maximize compound amount
      [~, xopt, ~, status] = lp_solve(-Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);

      if (status==1)
        % ok
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


      Res.Values(iM,iT,2) = max(Prob.Vars.ImbMets.S(iM,1:Prob.Vars.ImbMets.Tind(iT))*xopt(1:Prob.Vars.ImbMets.Tind(iT)),0); 

      Res.Status(iM,iT,2) = status;

      if Res.Values(iM,iT,2)>0 
        % an upper bounds >0 does not imply the lower bound 
        % Minimize compound amount
        [~, xopt, ~, status] = lp_solve(Objective, Prob.S, Prob.b, Prob.lb, Prob.ub, Prob.lecon, [], lp_options, Prob.x0);

        if (status==1)
          % ok
        elseif ismember(status, OK_STATUS)
          if log_file, fprintf(log_file, '[S:%i]', status); end 
        else  
          % this probably happens for the same reason as mentioned above
          if ~exist('tmp','dir'), mkdir('tmp'); end; save 'tmp/ProblemObjective.mat' Objective; 
          error('cFBA:ProblemObjectiveDetected', 'Unexpected status (%i).', status); 
        end  
      end

      Res.Values(iM,iT,1) = max(Prob.Vars.ImbMets.S(iM,1:Prob.Vars.ImbMets.Tind(iT))*xopt(1:Prob.Vars.ImbMets.Tind(iT)),0);
      Res.Status(iM,iT,1) = status;
      if log_file
        % show some progress in log file
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
