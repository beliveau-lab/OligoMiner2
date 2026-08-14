#!/usr/bin/env python3
"""OM2 BiLSTM v2 — the production model. Tokenizer, condition encoder, network, train/predict.

    python om2_bilstm.py        # self-checks: tokenizer conformance, width-freedom, BCE path

ONE MODEL, NOT A ZOO. This study builds and characterises exactly one architecture. Variant
exploration finished in `20260727_c`; `duplex-paircol` v1 (the 96-wide soft-clip variant) is
DELIBERATELY ABSENT from this directory and must not be reintroduced here.

WHY v2 (tokens only, 48 wide) rather than v1 (96 wide with a soft-clip side channel):
  * the side channel was measured INFORMATION-FREE upstream -- soft clips are strictly terminal
    (0.000% internal), 93.6% of rows are left-clipped, and corr(soft0, pDup) = +0.0009
  * v1 delivered only `soft0`, ONE column broadcast as a constant bias, so 47 of its 48 flags were
    computed and discarded
  * `side_channels: 0` removes the token/side split, so a width mismatch becomes a loud shape error
    instead of confident nonsense -- a whole silent failure mode deleted
  * half the input width, and this path is encode-bound on GPU, so the tokenizer cost is real

  This was a DECISION, not a measured win: the head-to-head that would have settled it ran under a
  defective objective (`20260727_c` ISSUES #12) and never completed on two of three substrates. The
  reasons above are independent of that comparison. Recorded honestly rather than dressed up.

THE OBJECTIVE IS SOFT-TARGET BCE, and that is not a style choice. pDup is a PROBABILITY. Squared
error treats an error of 0.05 as equally costly at p = 0.02 and at p = 0.95; BCE punishes confident
wrongness near the extremes, which is exactly where pDup mass sits. `20260727_c` trained MSE by
mistake -- its effect size was never measured, so this study does not inherit a claim about it,
only the corrected default.

CONDITIONS. Channel ORDER is a contract: `ncond=k` uses the first k of (t_eff, length, log10_sodium)
and nothing else. A net fitted with length in channel 1 and served with temperature there is not
detectably broken -- it is just wrong. The spec is written into every artifact and asserted on load.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


# ---- tokenizer ----------------------------------------------------------------------------------
BASES, PAD = "ACGT-", 0
VOCAB = 1 + len(BASES) ** 2          # 26 = PAD + 25 ordered (probe base, target base) pair states
_IDX = {b: i for i, b in enumerate(BASES)}


def tokenize(df, width):
    """One token per aligned duplex COLUMN = the ORDERED base pair. token = 1 + 5*probe + target.

    THE ALIGNMENT IS GIVEN, SO DO NOT MAKE THE MODEL REDISCOVER IT. A cross-encoder over
    `[CLS] probe [SEP] target [SEP]` spends capacity learning which probe base sits opposite which
    target base. The CIGAR already answered that. Tokenizing column-wise hands the model base-PAIR
    states directly, so position in the sequence IS position in the duplex.

    A specific mismatch is its own token: G opposite T is not the same state as G opposite A, and a
    bulge on either strand is distinct again. 25 pair states + PAD.
    """
    tok = np.zeros((len(df), width), dtype=np.int64)
    for i, (pa, ta) in enumerate(zip(df.probe_aln.values, df.target_aln.values)):
        for j in range(min(len(pa), width)):
            p, t = _IDX.get(pa[j]), _IDX.get(ta[j])
            if p is not None and t is not None:
                tok[i, j] = 1 + p * len(BASES) + t
    return tok


# ---- conditions ---------------------------------------------------------------------------------
T_LO, T_HI = 7.0, 97.0            # C, effective temperature (T + 0.65 x %formamide)
L_LO, L_HI = 10, 80               # nt, probe length
NA_LO, NA_HI = 0.05, 1.00         # M sodium, normalised in log10 space
CHANNELS = ["t_eff", "length", "log10_sodium"]
CONDITION_SPEC = {
    "version": "1", "channels": CHANNELS,
    "normalization": "min-max to [0,1] over fixed ranges; sodium in log10 space",
    "ranges": {"t_eff": [T_LO, T_HI], "length": [L_LO, L_HI], "sodium_M": [NA_LO, NA_HI]},
    "note": "channel ORDER is part of the contract; ncond=k uses the first k channels",
}
THRESHOLD = 0.2                   # a binder is pDup >= 0.2, everywhere


def _mm(v, lo, hi):
    return np.clip((np.asarray(v, dtype=np.float32) - lo) / (hi - lo), 0.0, 1.0)


def encode_conditions(df, ncond):
    """(n, ncond) float32, channels in CONDITION_SPEC order. ncond=0 returns an (n, 0) array.

    Min-max, not `/100`: raw scaling leaves temperature at 0.07-0.97, length at 0.10-0.80 and
    log-salt NEGATIVE, so the channels reach the first layer at different scales. Salt enters as
    log10 because the nearest-neighbour salt correction is a ln[Na+] term -- handing the network a
    logarithm it would otherwise have to relearn.
    """
    if ncond == 0:
        return np.zeros((len(df), 0), dtype=np.float32)
    cols = [_mm(df["label_celsius"].values, T_LO, T_HI),
            _mm(df["length"].values, L_LO, L_HI),
            _mm(np.log10(np.clip(df["label_sodium"].values, 1e-6, None)),
                np.log10(NA_LO), np.log10(NA_HI))]
    return np.stack(cols[:ncond], axis=1).astype(np.float32)


# ---- the network --------------------------------------------------------------------------------
DEFAULT = {"embed_dim": 16, "hidden": 64, "layers": 2, "dropout": 0.3,
           "lr": 3e-4, "weight_decay": 1e-2, "batch_size": 512,
           "epochs": 30, "patience": 5, "val_frac": 0.1, "seed": 0, "loss": "bce"}


def _net(hp, ncond):
    import torch, torch.nn as nn

    class Net(nn.Module):
        """SINGLE-STREAM, COLUMN-LEVEL CROSS-ENCODER over base-pair tokens. NOT a Siamese net.

        The name `siamese-bilstm` in the older registry is a misnomer kept there for artifact
        provenance; a true two-tower net collapses to near-random on pDup (~0.43 PR-AUC vs 0.94).
        Nothing in this directory carries that name.

        MASKED, PACKED POOLING IS LOAD-BEARING -- and it is also why this architecture is WIDTH-FREE.
        An unmasked `mean(dim=1)` over a padded frame makes most of the pooled vector the LSTM's
        response to PAD, diluting signal and letting the model read length off the dilution factor
        as a spurious feature. Here: true lengths from the pad token, `pack_padded_sequence` so the
        recurrence never sees padding at all, masked mean+max pooling. Max rides alongside mean
        because one catastrophic column (a central bulge) is a max-like event a mean washes out.

        Consequence: an LSTM has no fixed input dimension and padding is never consumed, so a net
        fitted at width 42 serves width 94 natively. Width is a batching convention here, never a
        capability limit, and must not appear in a model card as one.

        CONDITION FUSION IS LATE, AND THAT IS AN OPEN QUESTION, NOT A SETTLED CHOICE. The condition
        vector is concatenated with the POOLED representation, so temperature can shift the output
        but cannot reshape how the sequence is read. On a pooled-temperature corpus the identical
        token sequence carries ten different labels, and a late-fused model hedges toward the middle
        -- the compression measured in `20260727_c`. Early fusion or FiLM-style modulation are the
        candidate fixes and are this study's one open modelling arm.
        """
        def __init__(s):
            super().__init__()
            s.emb = nn.Embedding(VOCAB, hp["embed_dim"], padding_idx=PAD)
            s.lstm = nn.LSTM(hp["embed_dim"], hp["hidden"], hp["layers"], batch_first=True,
                             bidirectional=True,
                             dropout=hp["dropout"] if hp["layers"] > 1 else 0.0)
            pooled = hp["hidden"] * 2 * 2 + ncond          # (mean, max) x bidirectional + conditions
            s.head = nn.Sequential(
                nn.Linear(pooled, 256), nn.ReLU(), nn.Dropout(hp["dropout"]),
                nn.Linear(256, 64), nn.ReLU(), nn.Dropout(hp["dropout"]),
                nn.Linear(64, 1), nn.Sigmoid())

        def forward(s, tok, cond):
            x = s.emb(tok)
            lens = (tok != PAD).sum(1).clamp(min=1)
            packed = nn.utils.rnn.pack_padded_sequence(x, lens.cpu(), batch_first=True,
                                                       enforce_sorted=False)
            out, _ = s.lstm(packed)
            out, _ = nn.utils.rnn.pad_packed_sequence(out, batch_first=True,
                                                      total_length=tok.shape[1])
            m = (tok != PAD).unsqueeze(-1).to(out.dtype)
            mean = (out * m).sum(1) / m.sum(1).clamp(min=1)
            mx = torch.nan_to_num(out.masked_fill(m == 0, float("-inf")).max(1).values, neginf=0.0)
            return s.head(torch.cat([mean, mx, cond], dim=1)).squeeze(-1)

    return Net()


class OM2BiLSTM:
    """The OM2 BiLSTM v2. `ncond=0` builds a condition-BLIND model — the ablation arm."""

    def __init__(self, params=None, ncond=2, device=None):
        import torch
        self.hp = dict(DEFAULT); self.hp.update(params or {})
        self.ncond = ncond
        self.device = torch.device(device) if device else torch.device(
            "cuda" if torch.cuda.is_available() else "cpu")

    def fit(self, tok, cond, y):
        import torch, torch.nn as nn
        from torch.utils.data import DataLoader, TensorDataset, Subset
        hp = self.hp
        torch.manual_seed(hp["seed"]); np.random.seed(hp["seed"])
        torch.cuda.manual_seed_all(hp["seed"])
        self.net = _net(hp, self.ncond).to(self.device)

        tt = torch.from_numpy(np.asarray(tok, np.int64))
        ct = torch.from_numpy(np.asarray(cond, np.float32).reshape(len(tok), self.ncond))
        yt = torch.from_numpy(np.asarray(y, np.float32))
        full = TensorDataset(tt, ct, yt)

        # INTERNAL validation split, carved from the TRAINING rows only. The corpus `val` split is
        # reserved for hyperparameter selection and `test` is the eval set; early stopping on either
        # would be the leak the probe-grouped split exists to prevent.
        g = torch.Generator().manual_seed(hp["seed"])
        perm = torch.randperm(len(yt), generator=g)
        nv = max(1, int(round(hp["val_frac"] * len(yt))))
        tr = DataLoader(Subset(full, perm[nv:].tolist()), batch_size=hp["batch_size"], shuffle=True)
        va = DataLoader(Subset(full, perm[:nv].tolist()), batch_size=hp["batch_size"])

        opt = torch.optim.AdamW(self.net.parameters(), lr=hp["lr"],
                                weight_decay=hp["weight_decay"])
        sched = torch.optim.lr_scheduler.OneCycleLR(
            opt, max_lr=hp["lr"], total_steps=max(1, hp["epochs"] * len(tr)), pct_start=0.05)
        # Soft-target BCE. NOTE its value has an IRREDUCIBLE FLOOR equal to the mean entropy of the
        # targets, so it is not comparable across corpora with different label distributions. It is
        # a training objective, never a reported metric. Kept out of autocast: BCELoss is explicitly
        # unsafe in fp16 (log of an underflowed probability is -inf).
        lossf = {"bce": nn.BCELoss(), "mse": nn.MSELoss()}[hp["loss"]]

        best, best_state, bad = float("inf"), None, 0
        self.history = []
        for _ in range(hp["epochs"]):
            self.net.train()
            for xb, cb, yb in tr:
                opt.zero_grad()
                lossf(self.net(xb.to(self.device), cb.to(self.device)),
                      yb.to(self.device)).backward()
                torch.nn.utils.clip_grad_norm_(self.net.parameters(), 1.0)
                opt.step(); sched.step()
            self.net.eval(); tot = n = 0.0
            with torch.no_grad():
                for xb, cb, yb in va:
                    p = self.net(xb.to(self.device), cb.to(self.device))
                    tot += float(lossf(p, yb.to(self.device))) * len(yb); n += len(yb)
            vl = tot / max(n, 1)
            self.history.append(round(vl, 6))
            if vl < best - 1e-6:
                best, bad = vl, 0
                best_state = {k: v.detach().clone() for k, v in self.net.state_dict().items()}
            else:
                bad += 1
                if bad >= hp["patience"]:
                    break
        if best_state:                              # BEST-validation weights, not the last ones
            self.net.load_state_dict(best_state)
        self.best_val, self.epochs_run = best, len(self.history)
        return self

    def predict(self, tok, cond, batch=4096):
        import torch
        self.net.eval()
        out = []
        with torch.no_grad():
            for i in range(0, len(tok), batch):
                t = torch.from_numpy(np.asarray(tok[i:i+batch], np.int64)).to(self.device)
                c = torch.from_numpy(np.asarray(cond[i:i+batch], np.float32)
                                     .reshape(len(t), self.ncond)).to(self.device)
                out.append(self.net(t, c).cpu().numpy())
        return np.clip(np.concatenate(out), 0.0, 1.0)

    def save(self, path):
        import torch
        torch.save({"state_dict": self.net.state_dict(), "hp": self.hp, "ncond": self.ncond,
                    "condition_spec": CONDITION_SPEC, "vocab": VOCAB, "pad": PAD,
                    "threshold": THRESHOLD}, str(path))

    @classmethod
    def load(cls, path, device=None):
        import torch
        ck = torch.load(str(path), map_location="cpu", weights_only=False)
        if ck.get("condition_spec") != CONDITION_SPEC:
            raise ValueError(f"{path} was fitted under a different condition encoder — serving it "
                             f"now would shift every channel silently.\n  artifact: "
                             f"{ck.get('condition_spec')}\n  current:  {CONDITION_SPEC}")
        obj = cls(params=ck["hp"], ncond=ck["ncond"], device=device)
        obj.net = _net(ck["hp"], ck["ncond"]).to(obj.device)
        obj.net.load_state_dict(ck["state_dict"])
        obj.net.eval()
        return obj
