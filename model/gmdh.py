"""White-box GMDH (Group Method of Data Handling) regressor.

This module is a self-contained, dependency-light implementation of the
multilayer GMDH algorithm used in the project.  It is kept separate from the
training/validation pipeline (``model.py``) so the core algorithm can be read,
reused and unit-tested on its own.

Algorithm outline
-----------------
GMDH approximates an unknown response ``y = f(x)`` by a Kolmogorov--Gabor
functional series, which is approximated in turn by *partial quadratic
polynomials* over pairs of inputs.  Each polynomial is one "neuron"; neurons
are stacked in layers; only the neurons that improve a hold-out (validation)
split are carried forward.  The final prediction is the average of the selected
neurons in the last layer (the standard GMDH "p_average" output).

For ``d`` inputs and a single neuron built on a pair ``(v_i, v_j)``::

    u = a0 + a1*v_i + a2*v_j + a3*v_i*v_j + a4*v_i^2 + a5*v_j^2

The six coefficients ``a0..a5`` are fitted by ridge-regularised least squares
(Tikhonov) on the *training* portion of an internal random split; the neuron is
kept if its *validation* mean-squared error is within ``(1 + alpha)`` of the
best validation error in the layer.  ``alpha`` is the GMDH regularisation
(selection) coefficient; the project uses ``alpha = 0.6``.

Why p_average (single-layer stable)
-----------------------------------
Deeper layers multiply the model's degree (each layer can double the polynomial
degree), which overfits quickly on small geotechnical datasets.  Averaging the
selected neurons of the last layer acts as a model-average and keeps a
``max_layers = 1`` model numerically stable while still exposing nonlinear
pair interactions and squares.
"""

from __future__ import annotations

import numpy as np


class GMDH:
    """Multilayer GMDH with partial quadratic polynomials over input pairs.

    Parameters
    ----------
    alpha : float
        Regularisation / selection coefficient.  A neuron is retained when its
        validation MSE is ``<= (1 + alpha) * best_validation_MSE``.
    ridge : float
        Tikhonov ridge added to the diagonal of the normal-equation matrix
        (numerical stability for near-collinear pairs).
    max_layers : int
        Maximum number of hidden layers (a value of 1 gives the classic
        single-layer GMDH used throughout this project).
    val_frac : float
        Fraction of rows held out *internally* for neuron selection.
    max_neurons : int
        Cap on the number of neurons carried forward to the next layer.
    random_state : int
        Seed for the internal train/validation split.

    Attributes set by :meth:`fit`
    -----------------------------
    mu_, sig_ : ndarray
        Standardisation of the inputs (``z = (x - mu) / sig``).
    ym_, ys_ : float
        Standardisation of the target (``yw = (y - ym) / ys``).
    layers : list[list[dict]]
        One list per layer; each dict stores ``{"i", "j", "coef"}``.
    output_layer : int
        Index of the layer whose averaged neurons form the output.
    """

    def __init__(
        self,
        alpha: float = 0.6,
        ridge: float = 10.0,
        max_layers: int = 1,
        val_frac: float = 0.30,
        max_neurons: int = 8,
        random_state: int = 1,
        output_mode: str = "average",
    ):
        self.alpha = alpha
        self.ridge = ridge
        self.max_layers = max_layers
        self.val_frac = val_frac
        self.max_neurons = max_neurons
        self.random_state = random_state
        self.output_mode = output_mode  # "average" (p_average) or "best"

    # -- fitting ---------------------------------------------------------
    def fit(self, X, y) -> "GMDH":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float).ravel()
        n, d = X.shape

        # standardise inputs and target
        self.mu_ = X.mean(axis=0)
        self.sig_ = X.std(axis=0)
        self.sig_[self.sig_ == 0] = 1.0
        Z = (X - self.mu_) / self.sig_

        self.ym_ = float(y.mean())
        self.ys_ = float(y.std())
        if self.ys_ < 1e-9:
            self.ys_ = 1.0
        yw = (y - self.ym_) / self.ys_

        # internal train/validation split for neuron selection
        rng = np.random.RandomState(self.random_state)
        perm = rng.permutation(n)
        n_val = max(1, int(round(self.val_frac * n)))
        tr = perm[n_val:]
        va = perm[:n_val]
        Ztr, yw_tr = Z[tr], yw[tr]
        Zva, yw_va = Z[va], yw[va]

        nodes_tr = [Ztr[:, i] for i in range(d)]
        nodes_va = [Zva[:, i] for i in range(d)]

        self.layers: list[list[dict]] = []
        self.val_mse_per_layer: list[float] = []
        self.best: tuple | None = None  # (val_mse, layer_idx, neuron_idx)

        for _layer in range(1, self.max_layers + 1):
            m = len(nodes_tr)
            cands = []
            for i in range(m):
                for j in range(i, m):
                    zi, zj = nodes_tr[i], nodes_tr[j]
                    A = np.column_stack(
                        [np.ones_like(zi), zi, zj, zi * zj, zi**2, zj**2]
                    )
                    w = np.linalg.pinv(A.T @ A + self.ridge * np.eye(6)) @ A.T @ yw_tr

                    vi, vj = nodes_va[i], nodes_va[j]
                    Av = np.column_stack(
                        [np.ones_like(vi), vi, vj, vi * vj, vi**2, vj**2]
                    )
                    mse = float(np.mean((Av @ w - yw_va) ** 2))
                    cands.append(
                        {"i": i, "j": j, "coef": w, "mse": mse,
                         "pred_tr": A @ w, "pred_va": Av @ w}
                    )

            cands.sort(key=lambda c: c["mse"])
            best_mse = cands[0]["mse"]
            thresh = (1.0 + self.alpha) * best_mse
            sel = [c for c in cands if c["mse"] <= thresh]
            if len(sel) > self.max_neurons:
                sel = sel[: self.max_neurons]

            layer = [{"i": c["i"], "j": c["j"], "coef": c["coef"]} for c in sel]
            self.layers.append(layer)
            self.val_mse_per_layer.append(best_mse)

            for k, c in enumerate(sel):
                if self.best is None or c["mse"] < self.best[0]:
                    self.best = (c["mse"], len(self.layers) - 1, k)

            nodes_tr = [c["pred_tr"] for c in sel]
            nodes_va = [c["pred_va"] for c in sel]

            if len(sel) < 2:  # cannot form pairs -> stop growing
                break

        self.output_layer = len(self.layers) - 1
        return self

    # -- inference -------------------------------------------------------
    def _eval_layers(self, Z: np.ndarray) -> list[np.ndarray]:
        """Forward pass returning one (n, n_nodes) array per layer, plus Z."""
        layer_vals = [Z]
        cur = Z
        for layer in self.layers:
            nxt = []
            for neuron in layer:
                zi = cur[:, neuron["i"]]
                zj = cur[:, neuron["j"]]
                a0, a1, a2, a3, a4, a5 = neuron["coef"]
                v = a0 + a1 * zi + a2 * zj + a3 * zi * zj + a4 * zi**2 + a5 * zj**2
                nxt.append(v)
            cur = np.column_stack(nxt) if nxt else np.empty((Z.shape[0], 0))
            layer_vals.append(cur)
        return layer_vals

    def predict(self, X) -> np.ndarray:
        """Return predictions in the original target units."""
        X = np.asarray(X, dtype=float)
        Z = (X - self.mu_) / self.sig_
        vals = self._eval_layers(Z)
        if self.output_mode == "best" and self.best is not None:
            li, ni = self.best[1], self.best[2]
            out = vals[li + 1][:, ni]  # single best neuron, std-y units
        else:
            final = vals[self.output_layer + 1]  # (n, n_neurons) in std-y units
            out = final.mean(axis=1)
        return out * self.ys_ + self.ym_

    def _hidden_activations(self, X) -> np.ndarray:
        """Return the last layer's neuron outputs in *standardised-y* units.

        Used by the gradient and SHAP sensitivity routines, which operate on a
        differentiable scalar summary of the network.
        """
        X = np.asarray(X, dtype=float)
        Z = (X - self.mu_) / self.sig_
        vals = self._eval_layers(Z)
        return vals[self.output_layer + 1]  # (n, n_neurons)

    def predict_std(self, X) -> np.ndarray:
        """Return the averaged last-layer output in standardised-y units."""
        return self._hidden_activations(X).mean(axis=1)

    # -- symbolic expansion (for an explicit, human-readable equation) ---
    def raw_polynomial(self) -> dict[tuple, float]:
        """Expand the fitted network into a raw-unit polynomial.

        Returns a dict mapping an exponent tuple ``(e1, ..., ed)`` to a
        coefficient so that ``y = sum(coef * prod(x_k^e_k))``.  Keys with
        negligible coefficients are dropped.  This is what makes the model
        "white-box": the whole network collapses to one readable polynomial.
        """
        d = len(self.mu_)
        zero = (0,) * d

        def unit(i):
            e = [0] * d
            e[i] = 1
            return {zero: -self.mu_[i] / self.sig_[i], tuple(e): 1.0 / self.sig_[i]}

        def _add(p, q):
            r = dict(p)
            for k, v in q.items():
                r[k] = r.get(k, 0.0) + v
            return {k: v for k, v in r.items() if abs(v) > 1e-12}

        def _mul(p, q):
            r = {}
            for k1, v1 in p.items():
                for k2, v2 in q.items():
                    k = tuple(a + b for a, b in zip(k1, k2))
                    r[k] = r.get(k, 0.0) + v1 * v2
            return {k: v for k, v in r.items() if abs(v) > 1e-12}

        def _scale(p, c):
            return {k: v * c for k, v in p.items()}

        polys = [[unit(i) for i in range(d)]]
        cur = polys[0]
        for layer in self.layers:
            nxt = []
            for neuron in layer:
                pi, pj = cur[neuron["i"]], cur[neuron["j"]]
                a0, a1, a2, a3, a4, a5 = neuron["coef"]
                poly = {zero: a0}
                poly = _add(poly, _scale(pi, a1))
                poly = _add(poly, _scale(pj, a2))
                poly = _add(poly, _scale(_mul(pi, pj), a3))
                poly = _add(poly, _scale(_mul(pi, pi), a4))
                poly = _add(poly, _scale(_mul(pj, pj), a5))
                nxt.append(poly)
            cur = nxt
            polys.append(cur)

        final_polys = polys[self.output_layer + 1]
        avg = {}
        for p in final_polys:
            for k, v in p.items():
                avg[k] = avg.get(k, 0.0) + v
        avg = {k: v / len(final_polys) for k, v in avg.items()}

        # un-standardise: y = ys * yw + ym
        out = {k: v * self.ys_ for k, v in avg.items()}
        out[zero] = out.get(zero, 0.0) + self.ym_
        return out

    def polynomial_text(self, names: list[str]) -> str:
        """Render :meth:`raw_polynomial` as a readable multi-term equation."""
        poly = self.raw_polynomial()
        terms = []
        for exp, coef in sorted(poly.items(), key=lambda kv: -sum(kv[0])):
            if sum(exp) == 0:
                terms.append(f"{coef:.6g}")
                continue
            parts = [f"{coef:.6g}"]
            for e, nm in zip(exp, names):
                if e == 1:
                    parts.append(nm)
                elif e > 1:
                    parts.append(f"{nm}^{e}")
            terms.append("*".join(parts))
        return "y = " + " + ".join(terms)
