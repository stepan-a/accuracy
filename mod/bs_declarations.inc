// Borrowing constraint (buffer stock) model, shared by all the formulations.
// Taken from tests/occbin/model_borrcon of Dynare: b is debt, limited to M times
// income, lb the multiplier of that limit.

var b c ec lb y;
varexo u;

parameters RHO BETA M R SIGMA GAMMAC;
