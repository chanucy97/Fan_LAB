"""Attention Cox: inner-only epoch selection; Breslow baseline from fitting rows."""
import random

import numpy as np
import torch
from torch import nn
from lifelines.utils import concordance_index

from analysis_core import HORIZON, DIMENSIONS, percentile

LEARNING_RATE = .02
WEIGHT_DECAY = .001


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_num_threads(1)


class FourScoreAttention(nn.Module):
    """Fixed orthogonal modality/position codes in the declared DIMENSIONS order."""
    def __init__(self):
        super().__init__()
        self.instance_embed = nn.Linear(1, 4)
        self.register_buffer("modality_position", torch.eye(4))
        self.attention_v = nn.Linear(4, 2)
        self.attention_u = nn.Linear(4, 2)
        self.attention_w = nn.Linear(2, 1)
        self.risk_head = nn.Linear(4, 1)

    def forward(self, scores):
        if scores.ndim != 2 or scores.shape[1] != 4:
            raise ValueError("Expected four ordered dimension scores")
        e = torch.relu(self.instance_embed(scores.unsqueeze(-1)) + self.modality_position.unsqueeze(0))
        w = torch.softmax(self.attention_w(torch.tanh(self.attention_v(e)) *
                          torch.sigmoid(self.attention_u(e))).squeeze(-1), dim=1)
        return self.risk_head((w.unsqueeze(-1) * e).sum(1)).squeeze(1), w


def cox_loss(risk, time, event):
    if event.sum().item() == 0:
        raise ValueError("No training events")
    loss = risk.new_zeros(())
    for t in torch.unique(time[event.bool()]):
        deaths = event.bool() & time.eq(t)
        loss -= risk[deaths].sum() - deaths.sum() * torch.logsumexp(risk[time.ge(t)], 0)
    return loss / event.sum()


def tensors(x, y):
    require_order(x)
    return (torch.tensor(x.to_numpy(), dtype=torch.float32),
            torch.tensor(y.rfs_time_months.to_numpy(), dtype=torch.float32),
            torch.tensor(y.rfs_event.to_numpy(), dtype=torch.float32))


def train(x, y, epochs, seed, tuning=None, patience=60):
    """tuning must be an inner heldout subset of the outer fitting patients."""
    seed_all(seed)
    model = FourScoreAttention()
    optim = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    a, time, event = tensors(x, y)
    best, best_epoch, remaining = -np.inf, 0, patience
    history = []
    for epoch in range(1, epochs + 1):
        model.train()
        optim.zero_grad(set_to_none=True)
        risk, _ = model(a)
        loss = cox_loss(risk, time, event)
        if not torch.isfinite(loss):
            raise ValueError("Nonfinite attention training loss")
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optim.step()
        if tuning is not None:
            tx, ty = tuning
            raw, _ = predict_raw(model, tx)
            score = float(concordance_index(ty.rfs_time_months, -raw, ty.rfs_event))
            history.append({"epoch": epoch, "inner_validation_c_index": score, "loss": float(loss.detach())})
            if score > best:
                best, best_epoch, remaining = score, epoch, patience
            else:
                remaining -= 1
                if remaining == 0:
                    break
    if tuning is not None:
        if best_epoch < 1:
            raise ValueError("No usable inner validation score")
        return best_epoch, history
    return model


def predict_raw(model, x):
    require_order(x)
    model.eval()
    with torch.no_grad():
        r, w = model(torch.tensor(x.to_numpy(), dtype=torch.float32))
    return r.numpy().astype(float), w.numpy().astype(float)


def require_order(x):
    if list(x.columns) != DIMENSIONS:
        raise ValueError(f"Attention columns must be exactly {DIMENSIONS}; never silently reorder")


def breslow_probability(fit_logrisk, y, apply_logrisk, horizon=HORIZON):
    """Breslow ties match the neural partial likelihood; no application outcomes."""
    shift = float(np.max(fit_logrisk))
    relative = np.exp(np.asarray(fit_logrisk) - shift)
    time, event = y.rfs_time_months.to_numpy(), y.rfs_event.to_numpy()
    baseline = 0.0
    for t in np.unique(time[(event == 1) & (time <= horizon)]):
        baseline += np.sum((time == t) & (event == 1)) / relative[time >= t].sum()
    cumulative = baseline * np.exp(np.clip(np.asarray(apply_logrisk) - shift, -700, 700))
    return -np.expm1(-cumulative), {"horizon": horizon, "shift": shift, "baseline_cumulative_hazard": float(baseline)}


def fit_predict(x, y, apply, epochs, seed):
    model = train(x, y, epochs, seed)
    reference, _ = predict_raw(model, x)
    raw, weights = predict_raw(model, apply)
    p, baseline = breslow_probability(reference, y, raw)
    return model, (raw, percentile(reference, raw), p), weights, reference, baseline
