import numpy as np

src = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_v0.bin"
out = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_v2_dense_float32.npy"

dtype = np.dtype([
    ("task_id", "u1"),
    ("label", "u1"),
    ("comp_id", "u1"),
    ("log2_w", "u1"),
    ("log2_h", "u1"),
    ("scan_pos", "<u2"),
    ("left_abs_bucket", "u1"),
    ("top_abs_bucket", "u1"),
    ("diag_abs_bucket", "u1"),
    ("nnz_count", "u1"),
    ("ctx_off", "u1"),
    ("left_nonzero", "u1"),
    ("top_nonzero", "u1"),
    ("diag_nonzero", "u1"),
    ("sum_abs_neighbor_bucket", "u1"),
    ("ctx_id", "<u2"),
    ("ctx_state", "u1"),
    ("mps", "u1"),
    ("range_before", "<u2"),
    ("lps_range", "<u2"),
    ("cabac_bitcost", "<f4"),
    ("poc", "<i4"),
], align=False)

arr = np.fromfile(src, dtype=dtype)

num_coeff = (2 ** arr["log2_w"].astype(np.int32)) * (2 ** arr["log2_h"].astype(np.int32))
scan_pos_norm_bucket = np.floor(8.0 * arr["scan_pos"] / np.maximum(num_coeff, 1)).astype(np.int32)
scan_pos_norm_bucket = np.clip(scan_pos_norm_bucket, 0, 7)

dense = np.column_stack([
    arr["task_id"].astype(np.float32),                  # 0
    arr["label"].astype(np.float32),                    # 1
    arr["comp_id"].astype(np.float32),                  # 2
    arr["log2_w"].astype(np.float32),                   # 3
    arr["log2_h"].astype(np.float32),                   # 4
    arr["scan_pos"].astype(np.float32),                 # 5
    scan_pos_norm_bucket.astype(np.float32),            # 6
    arr["left_abs_bucket"].astype(np.float32),          # 7
    arr["top_abs_bucket"].astype(np.float32),           # 8
    arr["diag_abs_bucket"].astype(np.float32),          # 9
    arr["nnz_count"].astype(np.float32),                # 10
    arr["ctx_off"].astype(np.float32),                  # 11
    arr["left_nonzero"].astype(np.float32),             # 12
    arr["top_nonzero"].astype(np.float32),              # 13
    arr["diag_nonzero"].astype(np.float32),             # 14
    arr["sum_abs_neighbor_bucket"].astype(np.float32),  # 15
    arr["ctx_id"].astype(np.float32),                   # 16
    arr["ctx_state"].astype(np.float32),                # 17
    arr["mps"].astype(np.float32),                      # 18
    (arr["range_before"].astype(np.float32) / 512.0),   # 19
    (arr["lps_range"].astype(np.float32) / 512.0),      # 20
    arr["cabac_bitcost"].astype(np.float32),            # 21
    arr["poc"].astype(np.float32),                      # 22
]).astype(np.float32)

np.save(out, dense)

print("saved:", out)
print("shape:", dense.shape)
print("dtype:", dense.dtype)