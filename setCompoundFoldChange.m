function [ Prob ] = setCompoundFoldChange( Prob, COMPOUND_FOLD_CHANGE )
%%setCompoundFoldChange adapts the compound fold change constraints of the
%%specified problem to the specified value
%
% Input:
%   Prob:                   LP problem structure 
%   COMPOUND_FOLD_CHANGE:   new value for \alpha
%
% Output:
%   Prob:                   Updated problem
%
  if ~exist('COMPOUND_FOLD_CHANGE', 'var'), 
    error('COMPOUND_FOLD_CHANGE is missing.');
  elseif isnan(COMPOUND_FOLD_CHANGE) || isempty(COMPOUND_FOLD_CHANGE) || ~isnumeric(COMPOUND_FOLD_CHANGE) || COMPOUND_FOLD_CHANGE<=0
    error('COMPOUND_FOLD_CHANGE must be a positive number.');
  end
  nMLi = length(Prob.model.ImbalancedMets);
  Prob.S(Prob.Cons.CompoundFoldChange, Prob.Vars.ImbMets.StartAmounts) = (COMPOUND_FOLD_CHANGE-1)*eye(nMLi);
  Prob.lecon(Prob.Cons.CompoundFoldChange) = false;
  Prob.CompoundFoldChange = COMPOUND_FOLD_CHANGE;
end

