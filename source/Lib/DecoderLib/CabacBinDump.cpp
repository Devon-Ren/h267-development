// source/Lib/DecoderLib/CabacBinDump.cpp
#include "CabacBinDump.h"

#include <cstdio>
#include <cstdlib>

thread_local CabacFeatCtx g_cabacFeatCtx;

namespace
{
  static FILE* g_fp = nullptr;

  static void closeDump()
  {
    if (g_fp)
    {
      std::fclose(g_fp);
      g_fp = nullptr;
    }
  }

  static void ensureOpen()
  {
    if (g_fp) return;

    const char* path = std::getenv("VTM_CABAC_DUMP_PATH");
    if (!path || !path[0])
    {
      path = "/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/stats/cabac_sig_gt1_v0.bin";
    }

    g_fp = std::fopen(path, "ab");
    if (!g_fp)
    {
      return;
    }

    std::atexit(closeDump);
  }
}

void dumpCabacBinRecord(
    uint16_t ctx_id,
    uint8_t  ctx_state,
    uint8_t  mps,
    uint8_t  label,
    uint16_t range_before,
    uint16_t lps_range,
    float    cabac_bitcost
)
{
  if (!g_cabacFeatCtx.valid)
  {
    return;
  }

  ensureOpen();
  if (!g_fp)
  {
    g_cabacFeatCtx.valid = 0;
    return;
  }

  CabacBinRecord rec{};
  rec.task_id         = g_cabacFeatCtx.task_id;
  rec.label           = label;
  rec.comp_id         = g_cabacFeatCtx.comp_id;
  rec.log2_w          = g_cabacFeatCtx.log2_w;
  rec.log2_h          = g_cabacFeatCtx.log2_h;
  rec.scan_pos        = g_cabacFeatCtx.scan_pos;
  rec.left_abs_bucket = g_cabacFeatCtx.left_abs_bucket;
  rec.top_abs_bucket  = g_cabacFeatCtx.top_abs_bucket;
  rec.diag_abs_bucket = g_cabacFeatCtx.diag_abs_bucket;
  rec.nnz_count       = g_cabacFeatCtx.nnz_count;
  rec.ctx_id          = ctx_id;
  rec.ctx_state       = ctx_state;
  rec.mps             = mps;
  rec.poc             = g_cabacFeatCtx.poc;
  rec.range_before    = range_before;
  rec.lps_range       = lps_range;
  rec.cabac_bitcost   = cabac_bitcost;
  rec.ctx_off                 = g_cabacFeatCtx.ctx_off;
  rec.left_nonzero            = g_cabacFeatCtx.left_nonzero;
  rec.top_nonzero             = g_cabacFeatCtx.top_nonzero;
  rec.diag_nonzero            = g_cabacFeatCtx.diag_nonzero;
  rec.sum_abs_neighbor_bucket = g_cabacFeatCtx.sum_abs_neighbor_bucket;

  std::fwrite(&rec, sizeof(CabacBinRecord), 1, g_fp);

  // very important: consume exactly one pending context
  g_cabacFeatCtx.valid = 0;
}