'''
##########################
#below is for the shared trunk two-head model with the original 17 features (no new features)
###########################
import math
import os
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split

DATA = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_v1_dense_float32.npy"
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
BATCH_SIZE = 65536
EPOCHS = 100
LR = 1e-3
WEIGHT_DECAY = 1e-5
SEED = 42

#column order checking:
# columns
# 0 task_id
# 1 label
# 2 comp_id
# 3 log2_w
# 4 log2_h
# 5 scan_pos
# 6 scan_pos_norm_bucket
# 7 left_abs_bucket
# 8 top_abs_bucket
# 9 diag_abs_bucket
# 10 nnz_count
# 11 ctx_id
# 12 ctx_state
# 13 mps
# 14 range_before_norm
# 15 lps_range_norm
# 16 cabac_bitcost
# 17 poc

torch.manual_seed(SEED)
np.random.seed(SEED)

class CabacDataset(Dataset):
    def __init__(self, path):
        self.X = np.load(path, mmap_mode="r")

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        row = self.X[idx]
        sample = {
            "task_id": int(row[0]),
            "label": float(row[1]),
            "comp_id": int(row[2]),
            "log2_w": int(row[3]),
            "log2_h": int(row[4]),
            "scan_pos": int(row[5]),
            "scan_pos_norm_bucket": int(row[6]),
            "left_abs_bucket": int(row[7]),
            "top_abs_bucket": int(row[8]),
            "diag_abs_bucket": int(row[9]),
            "nnz_count": int(row[10]),
            "ctx_id": int(row[11]),
            "ctx_state": int(row[12]),
            "mps": int(row[13]),
            "range_before_norm": float(row[14]),
            "lps_range_norm": float(row[15]),
            "cabac_bitcost": float(row[16]),
        }
        return sample


def collate_fn(batch):
    out = {}
    for k in batch[0].keys():
        vals = [x[k] for x in batch]
        if k in ["label", "range_before_norm", "lps_range_norm", "cabac_bitcost"]:
            out[k] = torch.tensor(vals, dtype=torch.float32)
        else:
            out[k] = torch.tensor(vals, dtype=torch.long)
    return out


class SharedTrunkTwoHead(nn.Module):
    def __init__(self):
        super().__init__()

        self.emb_task = nn.Embedding(2, 2)
        self.emb_comp = nn.Embedding(3, 2)
        self.emb_log2w = nn.Embedding(8, 3)
        self.emb_log2h = nn.Embedding(8, 3)
        self.emb_scan = nn.Embedding(1024, 8)
        self.emb_scan_bucket = nn.Embedding(8, 3)
        self.emb_left = nn.Embedding(8, 3)
        self.emb_top = nn.Embedding(8, 3)
        self.emb_diag = nn.Embedding(8, 3)
        self.emb_nnz = nn.Embedding(8, 2)
        self.emb_ctxid = nn.Embedding(2048, 16)
        self.emb_ctxstate = nn.Embedding(256, 8)
        self.emb_mps = nn.Embedding(2, 2)

        emb_dim = 2 + 2 + 3 + 3 + 8 + 3 + 3 + 3 + 3 + 2 + 16 + 8 + 2
        cont_dim = 2
        in_dim = emb_dim + cont_dim

        self.trunk = nn.Sequential(
            nn.Linear(in_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
        )

        self.sig_head = nn.Linear(64, 1)
        self.gt1_head = nn.Linear(64, 1)

    def forward(self, batch):
        x = torch.cat([
            self.emb_task(batch["task_id"]),
            self.emb_comp(batch["comp_id"]),
            self.emb_log2w(batch["log2_w"].clamp(0, 7)),
            self.emb_log2h(batch["log2_h"].clamp(0, 7)),
            self.emb_scan(batch["scan_pos"].clamp(0, 1023)),
            self.emb_scan_bucket(batch["scan_pos_norm_bucket"].clamp(0, 7)),
            self.emb_left(batch["left_abs_bucket"].clamp(0, 7)),
            self.emb_top(batch["top_abs_bucket"].clamp(0, 7)),
            self.emb_diag(batch["diag_abs_bucket"].clamp(0, 7)),
            self.emb_nnz(batch["nnz_count"].clamp(0, 7)),
            self.emb_ctxid(batch["ctx_id"].clamp(0, 2047)),
            self.emb_ctxstate(batch["ctx_state"].clamp(0, 255)),
            self.emb_mps(batch["mps"].clamp(0, 1)),
            batch["range_before_norm"].unsqueeze(1),
            batch["lps_range_norm"].unsqueeze(1),
        ], dim=1)

        h = self.trunk(x)
        sig_logits = self.sig_head(h).squeeze(1)
        gt1_logits = self.gt1_head(h).squeeze(1)
        return sig_logits, gt1_logits


def bits_from_logits(logits, labels):
    # per-sample model bits in base-2
    probs = torch.sigmoid(logits).clamp(1e-6, 1 - 1e-6)
    bits = torch.where(labels > 0.5, -torch.log2(probs), -torch.log2(1.0 - probs))
    return bits


@torch.no_grad()
def evaluate(model, loader):
    model.eval()

    total_n = 0
    total_model_bits = 0.0
    total_cabac_bits = 0.0

    sig_n = gt1_n = 0
    sig_model = sig_cabac = 0.0
    gt1_model = gt1_cabac = 0.0

    for batch in loader:
        for k in batch:
            batch[k] = batch[k].to(DEVICE, non_blocking=True)

        sig_logits, gt1_logits = model(batch)

        labels = batch["label"]
        task_id = batch["task_id"]
        cabac_bits = batch["cabac_bitcost"]

        sig_mask = (task_id == 0)
        gt1_mask = (task_id == 1)

        model_bits = torch.zeros_like(labels)

        if sig_mask.any():
            model_bits[sig_mask] = bits_from_logits(sig_logits[sig_mask], labels[sig_mask])
            sig_n += int(sig_mask.sum())
            sig_model += model_bits[sig_mask].sum().item()
            sig_cabac += cabac_bits[sig_mask].sum().item()

        if gt1_mask.any():
            model_bits[gt1_mask] = bits_from_logits(gt1_logits[gt1_mask], labels[gt1_mask])
            gt1_n += int(gt1_mask.sum())
            gt1_model += model_bits[gt1_mask].sum().item()
            gt1_cabac += cabac_bits[gt1_mask].sum().item()

        total_n += labels.numel()
        total_model_bits += model_bits.sum().item()
        total_cabac_bits += cabac_bits.sum().item()

    metrics = {
        "overall_model_bits": total_model_bits / total_n,
        "overall_cabac_bits": total_cabac_bits / total_n,
        "overall_delta_bits": (total_model_bits - total_cabac_bits) / total_n,
        "sig_model_bits": sig_model / max(sig_n, 1),
        "sig_cabac_bits": sig_cabac / max(sig_n, 1),
        "sig_delta_bits": (sig_model - sig_cabac) / max(sig_n, 1),
        "gt1_model_bits": gt1_model / max(gt1_n, 1),
        "gt1_cabac_bits": gt1_cabac / max(gt1_n, 1),
        "gt1_delta_bits": (gt1_model - gt1_cabac) / max(gt1_n, 1),
        "n_total": total_n,
        "n_sig": sig_n,
        "n_gt1": gt1_n,
    }
    return metrics


def train():
    ds = CabacDataset(DATA)
    n = len(ds)
    n_train = int(0.9 * n)
    n_val = n - n_train

    train_ds, val_ds = random_split(ds, [n_train, n_val], generator=torch.Generator().manual_seed(SEED))

    train_loader = DataLoader(
        train_ds,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=4,
        pin_memory=True,
        collate_fn=collate_fn,
        drop_last=False,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=4,
        pin_memory=True,
        collate_fn=collate_fn,
        drop_last=False,
    )

    model = SharedTrunkTwoHead().to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    for epoch in range(1, EPOCHS + 1):
        model.train()
        running_loss = 0.0
        seen = 0

        for batch in train_loader:
            for k in batch:
                batch[k] = batch[k].to(DEVICE, non_blocking=True)

            optimizer.zero_grad()

            sig_logits, gt1_logits = model(batch)
            labels = batch["label"]
            task_id = batch["task_id"]

            sig_mask = (task_id == 0)
            gt1_mask = (task_id == 1)

            losses = []

            if sig_mask.any():
                sig_bits = bits_from_logits(sig_logits[sig_mask], labels[sig_mask])
                losses.append(sig_bits.mean())

            if gt1_mask.any():
                gt1_bits = bits_from_logits(gt1_logits[gt1_mask], labels[gt1_mask])
                losses.append(gt1_bits.mean())

            loss = sum(losses) / len(losses)
            loss.backward()
            optimizer.step()

            running_loss += loss.item() * labels.numel()
            seen += labels.numel()

        train_loss = running_loss / seen
        val_metrics = evaluate(model, val_loader)

        print(f"Epoch {epoch:02d}")
        print(f"  train_model_bits   = {train_loss:.6f}")
        print(f"  val_model_bits     = {val_metrics['overall_model_bits']:.6f}")
        print(f"  val_cabac_bits     = {val_metrics['overall_cabac_bits']:.6f}")
        print(f"  val_delta_bits     = {val_metrics['overall_delta_bits']:.6f}")
        print(f"  sig_model_bits     = {val_metrics['sig_model_bits']:.6f}")
        print(f"  sig_cabac_bits     = {val_metrics['sig_cabac_bits']:.6f}")
        print(f"  sig_delta_bits     = {val_metrics['sig_delta_bits']:.6f}")
        print(f"  gt1_model_bits     = {val_metrics['gt1_model_bits']:.6f}")
        print(f"  gt1_cabac_bits     = {val_metrics['gt1_cabac_bits']:.6f}")
        print(f"  gt1_delta_bits     = {val_metrics['gt1_delta_bits']:.6f}")
        print(f"  n_total={val_metrics['n_total']}  n_sig={val_metrics['n_sig']}  n_gt1={val_metrics['n_gt1']}")

    torch.save(model.state_dict(), "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_shared_trunk.pt")
    print("saved model.")


if __name__ == "__main__":
    train()
'''



################################
#below is for the shared trunk two-head model with five new features (ctx_off, left_nonzero, top_nonzero, diag_nonzero, sum_abs_neighbor_bucket)
################################
import os
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, random_split


DATA_PATH = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_v2_dense_float32.npy"
SAVE_DIR = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/result"
MODEL_PATH = os.path.join(SAVE_DIR, "cabac_sig_gt1_shared_trunk.pt")

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
BATCH_SIZE = 65536
EPOCHS = 50
LR = 1e-3
WEIGHT_DECAY = 1e-5
SEED = 42
NUM_WORKERS = 4

os.makedirs(SAVE_DIR, exist_ok=True)

torch.manual_seed(SEED)
np.random.seed(SEED)

#column order checking for new features:
# 0  task_id
# 1  label
# 2  comp_id
# 3  log2_w
# 4  log2_h
# 5  scan_pos
# 6  scan_pos_norm_bucket
# 7  left_abs_bucket
# 8  top_abs_bucket
# 9  diag_abs_bucket
# 10 nnz_count
# 11 ctx_off
# 12 left_nonzero
# 13 top_nonzero
# 14 diag_nonzero
# 15 sum_abs_neighbor_bucket
# 16 ctx_id
# 17 ctx_state
# 18 mps
# 19 range_before_norm
# 20 lps_range_norm
# 21 cabac_bitcost
# 22 poc

class CabacDataset(Dataset):
    def __init__(self, path: str):
        self.X = np.load(path, mmap_mode="r")

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        row = self.X[idx]
        return {
            "task_id": int(row[0]),
            "label": float(row[1]),
            "comp_id": int(row[2]),
            "log2_w": int(row[3]),
            "log2_h": int(row[4]),
            "scan_pos": int(row[5]),
            "scan_pos_norm_bucket": int(row[6]),
            "left_abs_bucket": int(row[7]),
            "top_abs_bucket": int(row[8]),
            "diag_abs_bucket": int(row[9]),
            "nnz_count": int(row[10]),
            "ctx_off": int(row[11]),
            "left_nonzero": int(row[12]),
            "top_nonzero": int(row[13]),
            "diag_nonzero": int(row[14]),
            "sum_abs_neighbor_bucket": int(row[15]),
            "ctx_id": int(row[16]),
            "ctx_state": int(row[17]),
            "mps": int(row[18]),
            "range_before_norm": float(row[19]),
            "lps_range_norm": float(row[20]),
            "cabac_bitcost": float(row[21]),
        }


def collate_fn(batch):
    out = {}
    float_keys = {"label", "range_before_norm", "lps_range_norm", "cabac_bitcost"}
    for k in batch[0].keys():
        vals = [x[k] for x in batch]
        if k in float_keys:
            out[k] = torch.tensor(vals, dtype=torch.float32)
        else:
            out[k] = torch.tensor(vals, dtype=torch.long)
    return out


class SharedTrunkTwoHead(nn.Module):
    def __init__(self):
        super().__init__()
        self.emb_task = nn.Embedding(2, 2)
        self.emb_comp = nn.Embedding(3, 2)
        self.emb_log2w = nn.Embedding(8, 3)
        self.emb_log2h = nn.Embedding(8, 3)
        self.emb_scan = nn.Embedding(1024, 8)
        self.emb_scan_bucket = nn.Embedding(8, 3)

        self.emb_left = nn.Embedding(8, 3)
        self.emb_top = nn.Embedding(8, 3)
        self.emb_diag = nn.Embedding(8, 3)
        self.emb_nnz = nn.Embedding(8, 2)

        self.emb_ctxoff = nn.Embedding(16, 4)
        self.emb_left_nz = nn.Embedding(2, 2)
        self.emb_top_nz = nn.Embedding(2, 2)
        self.emb_diag_nz = nn.Embedding(2, 2)
        self.emb_sum_abs = nn.Embedding(16, 4)

        self.emb_ctxid = nn.Embedding(2048, 16)
        self.emb_ctxstate = nn.Embedding(256, 8)
        self.emb_mps = nn.Embedding(2, 2)

        emb_dim = (
            2 + 2 + 3 + 3 + 8 + 3 +   
            3 + 3 + 3 + 2 +           
            4 + 2 + 2 + 2 + 4 +      
            16 + 8 + 2               
        )
        cont_dim = 2                  
        in_dim = emb_dim + cont_dim

        self.trunk = nn.Sequential(
            nn.Linear(in_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
        )

        self.sig_head = nn.Linear(64, 1)
        self.gt1_head = nn.Linear(64, 1)

    def forward(self, batch):
        x = torch.cat([
            self.emb_task(batch["task_id"].clamp(0, 1)),
            self.emb_comp(batch["comp_id"].clamp(0, 2)),
            self.emb_log2w(batch["log2_w"].clamp(0, 7)),
            self.emb_log2h(batch["log2_h"].clamp(0, 7)),
            self.emb_scan(batch["scan_pos"].clamp(0, 1023)),
            self.emb_scan_bucket(batch["scan_pos_norm_bucket"].clamp(0, 7)),

            self.emb_left(batch["left_abs_bucket"].clamp(0, 7)),
            self.emb_top(batch["top_abs_bucket"].clamp(0, 7)),
            self.emb_diag(batch["diag_abs_bucket"].clamp(0, 7)),
            self.emb_nnz(batch["nnz_count"].clamp(0, 7)),

            self.emb_ctxoff(batch["ctx_off"].clamp(0, 15)),
            self.emb_left_nz(batch["left_nonzero"].clamp(0, 1)),
            self.emb_top_nz(batch["top_nonzero"].clamp(0, 1)),
            self.emb_diag_nz(batch["diag_nonzero"].clamp(0, 1)),
            self.emb_sum_abs(batch["sum_abs_neighbor_bucket"].clamp(0, 15)),

            self.emb_ctxid(batch["ctx_id"].clamp(0, 2047)),
            self.emb_ctxstate(batch["ctx_state"].clamp(0, 255)),
            self.emb_mps(batch["mps"].clamp(0, 1)),

            batch["range_before_norm"].unsqueeze(1),
            batch["lps_range_norm"].unsqueeze(1),
        ], dim=1)

        h = self.trunk(x)
        sig_logits = self.sig_head(h).squeeze(1)
        gt1_logits = self.gt1_head(h).squeeze(1)
        return sig_logits, gt1_logits


def bits_from_logits(logits: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    probs = torch.sigmoid(logits).clamp(1e-6, 1 - 1e-6)
    return torch.where(labels > 0.5, -torch.log2(probs), -torch.log2(1.0 - probs))


@torch.no_grad()
def evaluate(model, loader):
    model.eval()

    total_n = 0
    total_model_bits = 0.0
    total_cabac_bits = 0.0

    sig_n = gt1_n = 0
    sig_model = sig_cabac = 0.0
    gt1_model = gt1_cabac = 0.0

    for batch in loader:
        for k in batch:
            batch[k] = batch[k].to(DEVICE, non_blocking=True)

        sig_logits, gt1_logits = model(batch)

        labels = batch["label"]
        task_id = batch["task_id"]
        cabac_bits = batch["cabac_bitcost"]

        sig_mask = (task_id == 0)
        gt1_mask = (task_id == 1)

        model_bits = torch.zeros_like(labels)

        if sig_mask.any():
            mb = bits_from_logits(sig_logits[sig_mask], labels[sig_mask])
            model_bits[sig_mask] = mb
            sig_n += int(sig_mask.sum())
            sig_model += mb.sum().item()
            sig_cabac += cabac_bits[sig_mask].sum().item()

        if gt1_mask.any():
            mb = bits_from_logits(gt1_logits[gt1_mask], labels[gt1_mask])
            model_bits[gt1_mask] = mb
            gt1_n += int(gt1_mask.sum())
            gt1_model += mb.sum().item()
            gt1_cabac += cabac_bits[gt1_mask].sum().item()

        total_n += labels.numel()
        total_model_bits += model_bits.sum().item()
        total_cabac_bits += cabac_bits.sum().item()

    metrics = {
        "overall_model_bits": total_model_bits / total_n,
        "overall_cabac_bits": total_cabac_bits / total_n,
        "overall_delta_bits": (total_model_bits - total_cabac_bits) / total_n,
        "sig_model_bits": sig_model / max(sig_n, 1),
        "sig_cabac_bits": sig_cabac / max(sig_n, 1),
        "sig_delta_bits": (sig_model - sig_cabac) / max(sig_n, 1),
        "gt1_model_bits": gt1_model / max(gt1_n, 1),
        "gt1_cabac_bits": gt1_cabac / max(gt1_n, 1),
        "gt1_delta_bits": (gt1_model - gt1_cabac) / max(gt1_n, 1),
        "n_total": total_n,
        "n_sig": sig_n,
        "n_gt1": gt1_n,
    }
    return metrics


def train():
    print(f"Loading dataset from: {DATA_PATH}")
    ds = CabacDataset(DATA_PATH)
    n = len(ds)
    n_train = int(0.9 * n)
    n_val = n - n_train

    train_ds, val_ds = random_split(
        ds, [n_train, n_val],
        generator=torch.Generator().manual_seed(SEED)
    )

    train_loader = DataLoader(
        train_ds,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=NUM_WORKERS,
        pin_memory=True,
        collate_fn=collate_fn,
        drop_last=False,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=NUM_WORKERS,
        pin_memory=True,
        collate_fn=collate_fn,
        drop_last=False,
    )

    model = SharedTrunkTwoHead().to(DEVICE)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    best_delta = float("inf")

    print(f"Device: {DEVICE}")
    print(f"Train size: {n_train}, Val size: {n_val}")

    for epoch in range(1, EPOCHS + 1):
        model.train()
        epoch_start = time.time()

        running_loss = 0.0
        seen = 0

        for batch in train_loader:
            for k in batch:
                batch[k] = batch[k].to(DEVICE, non_blocking=True)

            optimizer.zero_grad()

            sig_logits, gt1_logits = model(batch)
            labels = batch["label"]
            task_id = batch["task_id"]

            sig_mask = (task_id == 0)
            gt1_mask = (task_id == 1)

            losses = []

            if sig_mask.any():
                sig_bits = bits_from_logits(sig_logits[sig_mask], labels[sig_mask])
                losses.append(sig_bits.mean())

            if gt1_mask.any():
                gt1_bits = bits_from_logits(gt1_logits[gt1_mask], labels[gt1_mask])
                losses.append(gt1_bits.mean())

            loss = sum(losses) / len(losses)
            loss.backward()
            optimizer.step()

            running_loss += loss.item() * labels.numel()
            seen += labels.numel()

        train_model_bits = running_loss / seen
        val_metrics = evaluate(model, val_loader)

        if val_metrics["overall_delta_bits"] < best_delta:
            best_delta = val_metrics["overall_delta_bits"]
            torch.save(model.state_dict(), MODEL_PATH)

        print(f"Epoch {epoch:02d} | time={time.time()-epoch_start:.2f}s")
        print(f"  train_model_bits   = {train_model_bits:.6f}")
        print(f"  val_model_bits     = {val_metrics['overall_model_bits']:.6f}")
        print(f"  val_cabac_bits     = {val_metrics['overall_cabac_bits']:.6f}")
        print(f"  val_delta_bits     = {val_metrics['overall_delta_bits']:.6f}")
        print(f"  sig_model_bits     = {val_metrics['sig_model_bits']:.6f}")
        print(f"  sig_cabac_bits     = {val_metrics['sig_cabac_bits']:.6f}")
        print(f"  sig_delta_bits     = {val_metrics['sig_delta_bits']:.6f}")
        print(f"  gt1_model_bits     = {val_metrics['gt1_model_bits']:.6f}")
        print(f"  gt1_cabac_bits     = {val_metrics['gt1_cabac_bits']:.6f}")
        print(f"  gt1_delta_bits     = {val_metrics['gt1_delta_bits']:.6f}")
        print(f"  n_total={val_metrics['n_total']}  n_sig={val_metrics['n_sig']}  n_gt1={val_metrics['n_gt1']}")
        print(f"  best_delta_bits    = {best_delta:.6f}")

    print(f"Saved best model to: {MODEL_PATH}")


if __name__ == "__main__":
    train()