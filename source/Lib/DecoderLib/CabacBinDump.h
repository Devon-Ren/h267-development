// source/Lib/DecoderLib/CabacBinDump.h
#pragma once

#include <cstdint>

struct CabacFeatCtx
{
  uint8_t  valid = 0;      // 1 means the next decodeBin() should dump one record
  uint8_t  task_id = 0;    // 0 = SIG_COEFF_MAP_FLAG, 1 = GT1_FLAG
  uint8_t  comp_id = 0;    // 0:Y 1:Cb 2:Cr
  uint8_t  log2_w = 0;
  uint8_t  log2_h = 0;

  uint16_t scan_pos = 0;

  uint8_t  left_abs_bucket = 0;
  uint8_t  top_abs_bucket  = 0;
  uint8_t  diag_abs_bucket = 0;
  uint8_t  nnz_count       = 0;   // among left/top/diag

  int32_t  poc = -1;              // optional, set if easy to get

  uint8_t  ctx_off = 0;
  uint8_t  left_nonzero = 0;
  uint8_t  top_nonzero  = 0;
  uint8_t  diag_nonzero = 0;
  uint8_t  sum_abs_neighbor_bucket = 0;
};

#pragma pack(push, 1)
struct CabacBinRecord
{
  uint8_t  task_id;
  uint8_t  label;
  uint8_t  comp_id;
  uint8_t  log2_w;
  uint8_t  log2_h;

  uint16_t scan_pos;

  uint8_t  left_abs_bucket;
  uint8_t  top_abs_bucket;
  uint8_t  diag_abs_bucket;
  uint8_t  nnz_count;

  uint8_t  ctx_off;
  uint8_t  left_nonzero;
  uint8_t  top_nonzero;
  uint8_t  diag_nonzero;
  uint8_t  sum_abs_neighbor_bucket;

  uint16_t ctx_id;
  uint8_t  ctx_state;
  uint8_t  mps;

  uint16_t range_before;
  uint16_t lps_range;

  float    cabac_bitcost;

  int32_t  poc;
};
#pragma pack(pop)
extern thread_local CabacFeatCtx g_cabacFeatCtx;

void dumpCabacBinRecord(
    uint16_t ctx_id,
    uint8_t  ctx_state,
    uint8_t  mps,
    uint8_t  label,
    uint16_t range_before,
    uint16_t lps_range,
    float    cabac_bitcost
);