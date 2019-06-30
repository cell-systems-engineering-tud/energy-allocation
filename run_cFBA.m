function [] = run_cFBA(model)
%% run_cFBA starts the simulation for specified model
% and stores the results in the subfolder results
%
% Input:
%   model:  model structure
%
  tic
  if ~exist('model', 'var'), error('cFBA:MissingParameter', 'Parameter ''model'' is missing.'); end
  
  if exist('log', 'dir'), log_file = fopen(sprintf('log/%s.log',mfilename),'w'); else log_file=0; end
  if log_file, fprintf(log_file, '%s\tComputation started.\n\n',time); end
  
  COMPOUND_FOLD_CHANGE_PRECISION = 1e-10;
  
  
  if ~exist('results','dir'), mkdir('results'); end
  
  res = struct();
  res.lp_options = initializeSolverOptions_cplex(); 
  
  res.Prob = create_cFBAProblem_FromModel(model);
  
  bolSimulationCompleted = false;
  
  while ~bolSimulationCompleted
    try
      if log_file, fprintf(log_file, '%s\tCompute maximal \\alpha ... \n',time); end
      [res.Prob] = setMaximalCompoundFoldChange(res.Prob,COMPOUND_FOLD_CHANGE_PRECISION,res.lp_options);

      
      % compute variabilities of the solution
      if log_file, fprintf(log_file, '%s\tCompute variabilities ... \n',time); end
      [res.CompoundStartVar] = getCompoundStartAmountVariabilities(res.Prob, res.lp_options);          
      [res.FluxVar] = getFluxVariablities(res.Prob, res.lp_options);                      
      [res.CompoundVar] = getCompoundVariabilities(res.Prob, res.lp_options);  
      bolSimulationCompleted = true;
    catch Ex
      if strcmp(Ex.identifier,'cFBA:ProblemObjectiveDetected')
        if log_file, fprintf(log_file, '%s\tProblem-Objective detected. Restart simulation ... \n',time); end
        load tmp/ProblemObjective.mat Objective;
        PROBLEM_OBJECTIVES_FILE = sprintf('tmp/ProblemObjectives_%s.mat', model.ID);
        if ~exist(PROBLEM_OBJECTIVES_FILE, 'file')
          PO = Objective;
        else
          load(PROBLEM_OBJECTIVES_FILE,'PO');
          if all(PO(:,end)==Objective), error('The last problem objective was reported again.'); end
          PO = [PO Objective];
        end
        save(PROBLEM_OBJECTIVES_FILE, 'PO');
      else
        rethrow(Ex);
      end
    end
  end
  
  save(sprintf('results/res_cFBA_%s.mat', model.ID), 'res');
  if log_file
    fprintf(log_file, '\n%s\tComputation finished.\n',time);
    fclose(log_file);
  end
  toc
end
%

function [options] = initializeSolverOptions_gurobi()
  options = struct();
  
  options.solver = 'gurobi';
  options.OutputFlag = 0;
  options.LogToConsole = 0;
  options.LogFile = '';
  options.NumericFocus = 3;
  options.OptimalityTol = 1e-9;
  options.FeasibilityTol = 1e-9; 
end
%
function [options] = initializeSolverOptions_cplex()
  options = struct();
  
  options.solver = 'cplex';
  
  
  options.tune.display = 0;
  options.simplex.display = 0;
  options.sifting.display = 0;
  options.conflict.display = 2; 
 
  options.emphasis.numerical = true; 
  options.simplex.tolerances.optimality = 1e-7; 
  options.simplex.tolerances.feasibility = 1e-7; 
      % if any(ismember(solver_configuration, 'NoSymmetryBreak')), options.preprocessing.symmetry = 0; end
  options.read.scale = -1; 
  options.options.preprocessing.reduce = 0; 
  
      % if any(ismember(solver_configuration, 'NoPresolve')), 
        options.preprocessing.numpass = 0; 
        options.preprocessing.presolve = 0;
        options.preprocessing.dual = -1;
      %end
end