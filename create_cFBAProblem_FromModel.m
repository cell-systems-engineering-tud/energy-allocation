function [ Prob ] = create_cFBAProblem_FromModel( model )
%% create_cFBAProblem_FromModel creates the cFBA LP for the specified model
%
%  Input:
%    model:  model structure to create the LP problem fLOADor
%
%  Output:
%    Prob:   Problem structure that contains the LP formulation.
%

  if ~is_model_valid(model), return; end
  t_model = transform_model(model); 
  
  [nMets, nRxns] = size(t_model.S);
  nT = length(t_model.dT);
  
  nMetsImb = length(t_model.ImbalancedMets);
  nConI = size(t_model.Constraints.ForAllIntervals,1);
  nConM = size(t_model.Constraints.Max,1);
  nConMin = size(t_model.Constraints.Min,1);
  nConeq0 = size(t_model.Constraints.Initial,1); 
  MaxStartT = t_model.Constraints.MaxStartT;
  MinStartT = t_model.Constraints.MinStartT;
  nCONM = size(t_model.Constraints.CONM,1);
  if isfield(t_model.Constraints, 'gamma')
    nConGamma = sum(~isnan(t_model.Constraints.gamma));
  else
    nConGamma = 0;
  end
  nConTP = size(t_model.Constraints.ForAllTimePoints,1);
  nCon0 = size(t_model.Constraints.ForStart,1);
  
  nRxnrev = sum(t_model.rev);
  
  cdT = cumsum(t_model.dT);
    
  % can the imbalanced compounds be consumed according to stoichiometric matrix
  bolCanImbMetBeConsumed = any([t_model.S(t_model.ImbalancedMets,:) -t_model.S(t_model.ImbalancedMets,t_model.rev)]<0,2);
  MetsImbConsumable = find(bolCanImbMetBeConsumed);
  nMetsImbConsumed = sum(bolCanImbMetBeConsumed);
  nMetsBal = nMets - nMetsImb;
  mets_balanced_log = true(nMets,1); 
  mets_balanced_log(t_model.ImbalancedMets) = false;
    
  nVars = nMetsImb + (nRxns + nRxnrev) * nT;
          
  nRows = nMetsBal*nT + ...          % balance constraints
          2*nRxnrev*nT + ...         % absolute fluxes
          (nConTP+nConI+nCONM)*nT + ...    % constraints for all intervall or timepoints
          nConM * (nT-MaxStartT) +...  % Max of acetate after x time, no last point
           nConMin * (nT-MinStartT) +...  %Min of po4 after x time
          nConGamma + ...            % constraints for varying light
          nCon0 + ...                % constraints for time point 0
          nConeq0 + ...              %Initial concentration, equality
          nMetsImbConsumed*nT + ...  % amount of imbalanced mets, which can be consumed) may never be negative
          nMetsImb + ...             % compound fold change constraint
          1;                         % biomass sum at the beginning
        
  
  % LP data structures are initialized
  S = sparse(nRows, nVars);
  lb = -Inf(nVars,1);
  ub = Inf(nVars,1);
  b = zeros(nRows,1);
  lecon = false(nRows,1);
  c = zeros(nVars,1);
  
  % structure for access on LP variables is initialized
  Vars = struct();
  Vars.Fluxes = struct(); Vars.Fluxes(nT).Indices = [];
  Vars.AbsFluxes = struct(); Vars.AbsFluxes(nT).Indices = [];
  Vars.ImbMets = struct();
  Vars.ImbMets.StartAmounts = 1:nMetsImb;
  % we have no explicit M^k (for k>0) variables here 
  % because they can be substituted by linear combinations of the other
  % variables
  Vars.ImbMets.S = repmat([t_model.S(t_model.ImbalancedMets,:) zeros(nMetsImb, nRxnrev)],1,nT);
  tmp = repmat(t_model.dT, 1, (nRxns+nRxnrev))';
  Vars.ImbMets.S = Vars.ImbMets.S .* repmat(tmp(:)',nMetsImb,1);
  clear tmp;
  Vars.ImbMets.S = [eye(nMetsImb) Vars.ImbMets.S];
  Vars.ImbMets.Tind = (nRxns+nRxnrev)*(1:nT)'+nMetsImb;
  
  for iT = 1:nT
    Vars.Fluxes(iT).Indices = nMetsImb + (nRxns + nRxnrev) * (iT-1) + (1:nRxns);
    Vars.AbsFluxes(iT).Indices = Vars.Fluxes(iT).Indices;
    Vars.AbsFluxes(iT).Indices(t_model.rev) = nMetsImb + (nRxns + nRxnrev) * (iT-1) + nRxns + (1:nRxnrev);
  end; clear iT;
  
  % structure for access on LP constraints is initialized
  Cons = struct;
  Cons.CompoundFoldChange = NaN(nMetsImb,1);
  
  last_used_row = 0;
  % insert the stoichiometric information
  for iT = 1:nT
    lb(Vars.Fluxes(iT).Indices) = t_model.lb;
    ub(Vars.Fluxes(iT).Indices) = t_model.ub;
    lb(Vars.AbsFluxes(iT).Indices) = max(0,t_model.lb);
    ub(Vars.AbsFluxes(iT).Indices) = max(abs(t_model.lb),abs(t_model.ub));
    
    S(last_used_row + (1:nMetsBal), Vars.Fluxes(iT).Indices) = t_model.S(mets_balanced_log,:);
        
    if nRxnrev>0
      % we need absolute fluxes for capacity constraints
      S(last_used_row + nMetsBal + (1:nRxnrev), Vars.Fluxes(iT).Indices(t_model.rev)) = eye(nRxnrev);
      S(last_used_row + nMetsBal + (1:nRxnrev), Vars.AbsFluxes(iT).Indices(t_model.rev)) = -eye(nRxnrev); 
      S(last_used_row + nMetsBal + nRxnrev + (1:nRxnrev), Vars.Fluxes(iT).Indices(t_model.rev)) = -eye(nRxnrev);
      S(last_used_row + nMetsBal + nRxnrev + (1:nRxnrev), Vars.AbsFluxes(iT).Indices(t_model.rev)) = -eye(nRxnrev); 
      lecon(last_used_row + nMetsBal + (1:2*nRxnrev)) = true;
    end
    last_used_row = last_used_row + nMetsBal + 2*nRxnrev;
  end; clear iT;
  
  %   insert constraints for all intervals 
  for iCon = 1:nConI
    for iT = 1:nT
      last_used_row = last_used_row + 1;
      lecon(last_used_row) = true;
      S(last_used_row, Vars.AbsFluxes(iT).Indices) = t_model.Constraints.ForAllIntervals(iCon,1:nRxns);
      if iT>1
        S(last_used_row, 1:Vars.ImbMets.Tind(iT-1)) = S(last_used_row, 1:Vars.ImbMets.Tind(iT-1)) + ...
                                       t_model.Constraints.ForAllIntervals(iCon,(nRxns+1):end) * Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(iT-1)); 
      else
        S(last_used_row, Vars.ImbMets.StartAmounts) = t_model.Constraints.ForAllIntervals(iCon,(nRxns+1):end);
      end
    end; clear iT;
  end; clear iCon;
  
  %DEGRADATION ACETATE
  start_ac = last_used_row +1;
    for iCon = 1:nCONM
    for iT = 1:nT
      last_used_row = last_used_row + 1;
      lecon(last_used_row) = true;
      S(last_used_row, Vars.AbsFluxes(iT).Indices) = t_model.Constraints.CONM(iCon,1:nRxns);
      if iT>1
        S(last_used_row, 1:Vars.ImbMets.Tind(iT-1)) = S(last_used_row, 1:Vars.ImbMets.Tind(iT-1)) + ...
                                       t_model.Constraints.CONM(iCon,(nRxns+1):end) * Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(iT-1)); 
      else
        S(last_used_row, Vars.ImbMets.StartAmounts) = t_model.Constraints.ForAllIntervals(iCon,(nRxns+1):end);
      end
    end; clear iT;
  end; clear iCon;
  end_ac = last_used_row;
  lecon(start_ac:end_ac) = true;
  %   insert constraints for all timepoints
  
  lecon(last_used_row + (1:nT*nConTP)) = true;
  for iT=1:nT
    S(last_used_row + (1:nConTP), 1:Vars.ImbMets.Tind(iT)) = t_model.Constraints.ForAllTimePoints * Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(iT));
    % b(nr + (1:nMetCon)) = - c_model.METCON_M * c_model.E0;
    % Cons.Mets(:, iT) = nr + (1:nMetCon)';
    last_used_row = last_used_row + nConTP;
  end; clear iT;
  
  %Max quota for acetate
  
  lecon(last_used_row + (1:(nT-MaxStartT-1)*nConM)) = true;
  for iT=MaxStartT:nT-1
    S(last_used_row + (1:nConM), 1:Vars.ImbMets.Tind(iT)) = -(t_model.Constraints.Max * Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(iT)));
    % b(nr + (1:nMetCon)) = - c_model.METCON_M * c_model.E0;
    % Cons.Mets(:, iT) = nr + (1:nMetCon)';
    last_used_row = last_used_row + nConM;
  end; clear iT;
  %MIn quota for P04 
  lecon(last_used_row + (1:(nT-MinStartT-1)*nConMin)) = true;
  for iT=MinStartT:nT-1
    S(last_used_row + (1:nConMin), 1:Vars.ImbMets.Tind(iT)) = t_model.Constraints.Min * Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(iT));
    % b(nr + (1:nMetCon)) = - c_model.METCON_M * c_model.E0;
    % Cons.Mets(:, iT) = nr + (1:nMetCon)';
    last_used_row = last_used_row + nConMin;
  end; clear iT;
  
  % insert constraints for simulation start
  if nCon0>0
    lecon(last_used_row + (1:nCon0)) = true;
    S(last_used_row+(1:nCon0),Vars.ImbMets.StartAmounts) = t_model.Constraints.ForStart;
    b(last_used_row + (1:nCon0)) = 0;
    last_used_row = last_used_row + nCon0;
  end
    % insert constraints for simulation start; equality
  if nConeq0>0
    lecon(last_used_row + (1:nConeq0)) = false;
    S(last_used_row+(1:nConeq0),Vars.ImbMets.StartAmounts) = t_model.Constraints.Initial(:,1:end-1);
    b(last_used_row + (1:nConeq0)) = t_model.Constraints.Initial(:,end);
    last_used_row = last_used_row + nConeq0;
  end
  
  if nConGamma>0
    lecon(last_used_row + (1:nConGamma)') = true;
    tmp_Pigment_log = ismember(t_model.mets(t_model.ImbalancedMets), 'Pigment');
    if ~any(tmp_Pigment_log), error('cFBA:InvalidInput', 'Compound ''Pigment'' not found.'); end
    tmp_Photon_Uptake_log = ismember(t_model.rxns, 'Photon_Uptake');
    if ~any(tmp_Photon_Uptake_log), error('cFBA:InvalidInput', 'Reaction ''Photon_Uptake'' not found.'); end
    for iT=find(~isnan(t_model.Constraints.gamma))'
      S(last_used_row + 1, Vars.Fluxes(iT).Indices(tmp_Photon_Uptake_log)) = t_model.Constraints.gamma(iT);
      if iT>1
        S(last_used_row + 1, 1:Vars.ImbMets.Tind(iT-1)) = S(last_used_row + 1, 1:Vars.ImbMets.Tind(iT-1)) + ...
                                                       -Vars.ImbMets.S(tmp_Pigment_log,1:Vars.ImbMets.Tind(iT-1));
      else
        S(last_used_row + 1, Vars.ImbMets.StartAmounts(tmp_Pigment_log)) = -1;
      end
      last_used_row = last_used_row + 1;
    end; clear iT;
  end
 
  
  % forbid non negative amounts of imbalanced metabolites
  for iMet = MetsImbConsumable'
    for iT = 1:nT
      last_used_row = last_used_row + 1;
      S(last_used_row, 1:Vars.ImbMets.Tind(iT)) = -Vars.ImbMets.S(iMet,1:Vars.ImbMets.Tind(iT));
      lecon(last_used_row) = true;
    end; clear iT;
  end; clear iMet;
  
  % consider varying bound information for intervals
  if isfield(t_model, 'ub_var')
    for cRxnName = fields(t_model.ub_var)'
      iRxn = find(ismember(t_model.rxns,cRxnName),1);
      strRxnName = cRxnName{1};
      for iT = 1:nT
        ub(Vars.Fluxes(iT).Indices(iRxn)) = t_model.ub_var.(strRxnName)(iT); 
      end; clear iT;
    end; clear cRxnName;
  end
  if isfield(t_model, 'lb_var')
    for cRxnName = fields(t_model.lb_var)'
      iRxn = find(ismember(t_model.rxns,cRxnName),1);
      strRxnName = cRxnName{1};
      for iT = 1:nT
        lb(Vars.Fluxes(iT).Indices(iRxn)) = t_model.lb_var.(strRxnName)(iT); 
      end; clear iT;
    end; clear cRxnName;
  end
  
  % Fix weighted sum of start amounts
  S(last_used_row + 1, Vars.ImbMets.StartAmounts) = t_model.Biomass;
  b(last_used_row + 1) = 1;
  last_used_row = last_used_row + 1;
  
  % compound fold change = 1;
  % reproduce startamounts constraint
  Cons.CompoundFoldChange(:) = last_used_row + (1:nMetsImb)';
  S(last_used_row+(1:nMetsImb), 1:Vars.ImbMets.Tind(nT))   = -Vars.ImbMets.S(:,1:Vars.ImbMets.Tind(nT));
  lecon(last_used_row+(1:nMetsImb)) = false;
  last_used_row = last_used_row + nMetsImb;
  
  % Set lb of startamounts
  lb(Vars.ImbMets.StartAmounts) = 0;
    
  Prob = struct();
  Prob.S = S;
  Prob.lb = lb;
  Prob.ub = ub;
  Prob.b = b;
  Prob.lecon = lecon;
  Prob.Vars = Vars;
  Prob.Cons = Cons;
  Prob.model = model;
end
%
function [res] = is_model_valid(model)
  % is_model_valid performs a short check on the model
  if ~all(model.ImbalancedMets(1:end-1)<model.ImbalancedMets(2:end))
    % this property is checked, because it might not cause a program stop
    % but a serious  miscalulation
    error('Property ''ImbalancedMets'' must have an increasing order.'); 
  end
  res = true;
end
%
function [model] = transform_model(model)
% transforms the model constraints in a different format
  nRxn = length(model.rxns);
  nMets = length(model.ImbalancedMets);
  CONI = zeros(0,nRxn + nMets);
  CONTP = zeros(0,nMets);
  CON0 = zeros(0,nMets);
  CONMax = zeros(0,nMets);
  CONMin = zeros(0,nMets);
  CONM = zeros(0,nRxn + nMets);
  CONeq0 = zeros(0,nMets+1);
  MaxStartT = 1;
  MinStartT = 1;
  if isfield(model, 'Constraints')
    if isfield(model.Constraints, 'Capacities')
      % Apply capacity constraints
      CONI = [CONI; [model.Constraints.Capacities.A -model.Constraints.Capacities.B]];
      if isfield(model.Constraints.Capacities, 'gamma')
        model.Constraints.gamma = model.Constraints.Capacities.gamma;
      end
      model.Constraints = rmfield(model.Constraints, 'Capacities'); 
    end
    if isfield(model.Constraints, 'Quotas')
      % apply quota constraints
      if isfield(model.Constraints.Quotas, 'ForAllTimePoints')
        CONTP = [CONTP; model.Constraints.Quotas.ForAllTimePoints.C*model.Biomass'-model.Constraints.Quotas.ForAllTimePoints.B];
      end
      if isfield(model.Constraints.Quotas, 'ForStart')
        CON0 = [CON0; model.Constraints.Quotas.ForStart.C*model.Biomass'-model.Constraints.Quotas.ForStart.B];
      end
      if isfield(model.Constraints.Quotas, 'Max')
        CONMax = [CONMax; model.Constraints.Quotas.Max.C*model.Biomass'-model.Constraints.Quotas.Max.B];
        MaxStartT = [model.Constraints.Quotas.Max.StartT];
      end
      if isfield(model.Constraints.Quotas, 'Min')
        CONMin = [CONMin; model.Constraints.Quotas.Min.C-model.Constraints.Quotas.Min.B];
        MinStartT = [model.Constraints.Quotas.Min.StartT];
      end
      if isfield(model.Constraints.Quotas, 'Initial')
        CONeq0 = [CONeq0; model.Constraints.Quotas.Initial.B model.Constraints.Quotas.Initial.C];
        
      end
      model.Constraints = rmfield(model.Constraints, 'Quotas');
    end 
    if isfield(model.Constraints, 'Maintenance')
      % apply maintenance constraints
      CONM = [CONM; [model.Constraints.Maintenance.C -model.Constraints.Maintenance.A]];
      model.Constraints = rmfield(model.Constraints, 'Maintenance');
    end
    if isfield(model.Constraints, 'Acdeg')
      % apply maintenance constraints
      CONI = [CONI; [-model.Constraints.Acdeg.A model.Constraints.Acdeg.C]];
      model.Constraints = rmfield(model.Constraints, 'Acdeg');
    end
  end
  model.Constraints.ForAllTimePoints = CONTP;
  model.Constraints.ForAllIntervals = CONI;
  model.Constraints.ForStart = CON0; 
  model.Constraints.Max = CONMax; 
  model.Constraints.Min = CONMin; 
  model.Constraints.Initial = CONeq0; 
  model.Constraints.MaxStartT = MaxStartT;
  model.Constraints.MinStartT = MinStartT;
  model.Constraints.CONM = CONM;
end
