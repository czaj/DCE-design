"""Compare predicted panel MXL standard errors with 600 independently fitted samples.

Run python tests/validate_panel_se.py with MATLAB on PATH, NumPy and SciPy.
The design has five binary tasks per respondent and one normally distributed
coefficient. All data are synthetic. This optional check uses Gaussian
quadrature to fit the model and the observed Hessian to estimate its covariance.
"""
import json
import subprocess
import tempfile
from itertools import product
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, roots_hermitenorm

root = Path(__file__).resolve().parents[1]
nodes, weights = roots_hermitenorm(42)
weights /= np.sqrt(2 * np.pi)
cards = np.arange(-2, 3)
y = np.array(list(product([0, 1], repeat=cards.size)))
truth = np.array([0.2, 0.45, 0.9])
n, reps = 1500, 600
rng = np.random.default_rng(291121)

def law(theta):
    v = theta[0] + (theta[1] + theta[2] * nodes[:, None]) * cards
    p = expit(v)
    weighted = np.exp(y @ v.T - np.logaddexp(0, v).sum(axis=1)) * weights
    mass = weighted.sum(axis=1)
    residual = y[:, None, :] - p
    gx = residual @ cards
    jac = np.column_stack([(weighted * residual.sum(axis=2)).sum(axis=1),
                           (weighted * gx).sum(axis=1), (weighted * gx * nodes).sum(axis=1)])
    return mass, jac

prob, deriv = law(truth)
assert np.isclose(prob.sum(), 1)
info = deriv.T @ (deriv / prob[:, None])
spec = {
    'rows': 5, 'alternatives': ['a', 'sq'], 'balance': 'none', 'dominance': False,
    'attributes': {'x': cards.tolist()}, 'fixed': {'sq': {'x': 0}},
    'model': {'type': 'rppanel', 'utility': {'a': 'asc + b_x*x', 'sq': 'b_x*x'},
              'priors': {'asc': truth[0], 'b_x': truth[1]}, 'random': {'b_x': 'normal(0.9)'}}
}
with tempfile.TemporaryDirectory() as temp:
    work = Path(temp)
    (work / 'spec.json').write_text(json.dumps(spec), encoding='utf-8')
    toolbox = root.as_posix().replace("'", "''")
    command = (
        f"restoredefaultpath; addpath('{toolbox}'); "
        "e=dce_evaluate('spec.json',(-2:2)','EvalPanelDraws',2048,'Seed',29,'Verbose',0); "
        "assert(isequal(e.models.parameters,{'asc';'b_x';'sd(b_x)'}) && e.models.identified); "
        "writematrix(e.models.information,'information.csv');"
    )
    subprocess.run(['matlab', '-batch', command], cwd=work, check=True)
    matlab = np.loadtxt(work / 'information.csv', delimiter=',')
relative = np.linalg.norm(info - matlab) / np.linalg.norm(info)
assert relative < 0.025, relative
pred = np.sqrt(np.diag(np.linalg.inv(n * matlab)))
estimates, observed_se = [], []
for counts in rng.multinomial(n, prob, size=reps):
    def objective(t):
        mass, jac = law(t)
        return -counts @ np.log(mass), -(counts / mass) @ jac
    fit = minimize(objective, truth, method='L-BFGS-B', jac=True,
                   bounds=[(-3,3),(-3,3),(0.03,3)], options={'ftol':1e-14,'gtol':1e-7})
    if not fit.success or np.linalg.norm(fit.jac) / n >= 1e-5:
        fit = minimize(objective, fit.x, method='BFGS', jac=True, options={'gtol':1e-6})
    assert np.linalg.norm(fit.jac) / n < 1e-5 and 0.03 < fit.x[2] < 3, (fit.message, fit.x, fit.jac)
    estimates.append(fit.x)
    steps = np.eye(3) * 1e-4
    h = np.column_stack([(objective(fit.x + step)[1] - objective(fit.x - step)[1]) / 2e-4 for step in steps])
    assert np.linalg.eigvalsh(h).min() > 0
    observed_se.append(np.sqrt(np.diag(np.linalg.inv((h+h.T)/2))))
empirical = np.std(estimates, axis=0, ddof=1)
print('relative information discrepancy', relative, flush=True)
print('parameter, predicted SE, empirical SD, ratio, mean observed-Hessian SE')
for row in zip(['asc','b_x','sd(b_x)'], pred, empirical, empirical/pred, np.mean(observed_se,axis=0)):
    print(*row, sep=',')
# A Monte Carlo standard deviation from 600 replications has roughly 3% relative sampling error.
assert np.all(np.abs(empirical/pred - 1) < 0.13)
