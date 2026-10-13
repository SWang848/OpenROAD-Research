// SPDX-License-Identifier: BSD-3-Clause
//
// lgista_test -- minimal standalone driver for the local graph ISTA engine.
//
// Usage:
//   lgista_test <design.odb>  <lib1.lib> [lib2.lib ...] [opts]
//   lgista_test <design.v>    <lib1.lib> [lib2.lib ...] --lef <f.lef> [--lef ...] --top <name> [opts]
//
//   opts: [--inst <inst_name>] [--test t1|t2a|t2b|all]
//
// Loading paths:
//   * .odb : read an OpenDB database directly.
//   * .v   : read LEF (tech + std cells) + Verilog, then link the design
//            (no placement / RC needed -- timing comes from Liberty).
//
// Then it either runs the correctness tests (--test) or the analyze demo.
// A self-contained fixture lives in test/chain.v (Nangate45):
//   lgista_test test/chain.v test/Nangate45/Nangate45_typ.lib \
//       --lef test/Nangate45/Nangate45.lef --top chain --test all

#include <tcl.h>

#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

#include "odb/db.h"
#include "odb/lefin.h"
#include "utl/Logger.h"

#include "db_sta/dbReadVerilog.hh"
#include "db_sta/dbSta.hh"
#include "sta/MinMax.hh"
#include "sta/Sta.hh"
#include "sta/VerilogReader.hh"

#include "lgista/LocalGraphIsta.hh"
#include "LgistaTests.hh"

namespace {
bool endsWith(const std::string& s, const char* suffix)
{
  size_t n = std::strlen(suffix);
  return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}
}  // namespace

int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: %s <design.odb|design.v> <lib1.lib> [lib2.lib ...] "
                 "[--lef <f.lef>]... [--top <name>] [--inst <name>] "
                 "[--test t1|t2a|t2b|all]\n",
                 argv[0]);
    return 1;
  }

  const std::string design_path = argv[1];
  std::vector<std::string> lib_paths;
  std::vector<std::string> lef_paths;
  std::string top = "chain";
  std::string single_inst;
  std::string test_mode;
  for (int i = 2; i < argc; ++i) {
    if (std::strcmp(argv[i], "--inst") == 0 && i + 1 < argc) {
      single_inst = argv[++i];
    } else if (std::strcmp(argv[i], "--test") == 0 && i + 1 < argc) {
      test_mode = argv[++i];
    } else if (std::strcmp(argv[i], "--lef") == 0 && i + 1 < argc) {
      lef_paths.emplace_back(argv[++i]);
    } else if (std::strcmp(argv[i], "--top") == 0 && i + 1 < argc) {
      top = argv[++i];
    } else {
      lib_paths.emplace_back(argv[i]);
    }
  }

  try {
    // --- 1. base components -------------------------------------------------
    Tcl_Interp* interp = Tcl_CreateInterp();
    utl::Logger* logger = new utl::Logger();
    odb::dbDatabase* db = odb::dbDatabase::create();
    db->setLogger(logger);

    // --- 2. dbSta (registers itself as a db observer) -----------------------
    sta::dbSta* sta = new sta::dbSta(interp, db, logger);

    const bool verilog_flow = endsWith(design_path, ".v");

    // --- 3a. LEF (verilog flow only) ---------------------------------------
    if (verilog_flow) {
      if (lef_paths.empty()) {
        std::fprintf(stderr, "error: .v design needs at least one --lef\n");
        return 1;
      }
      odb::lefin lef_reader(db, logger, /*ignore_non_routing_layers=*/false);
      odb::dbLib* lib = nullptr;
      for (size_t i = 0; i < lef_paths.size(); ++i) {
        std::printf("Reading LEF: %s\n", lef_paths[i].c_str());
        if (i == 0)
          lib = lef_reader.createTechAndLib("lgista_tech", "lgista_lib",
                                            lef_paths[i].c_str());
        else
          lef_reader.updateLib(lib, lef_paths[i].c_str());
      }
    }

    // --- 3b. Liberty -------------------------------------------------------
    for (const std::string& lib : lib_paths) {
      std::printf("Reading liberty: %s\n", lib.c_str());
      // Read max only -> a single dcalc analysis point (index 0), so the
      // PtGraph ap (LocalSta uses MinMax::max) and the global graph share one
      // consistent index. (Reading MinMaxAll::all() yields two aps and the
      // PtGraph's local ap storage doesn't line up with the global index.)
      sta->readLiberty(lib.c_str(),
                       sta->cmdScene(),
                       sta::MinMaxAll::max(),
                       /*infer_latches=*/false);
    }

    // --- 3c. design (.odb or .v) -------------------------------------------
    if (verilog_flow) {
      std::printf("Reading verilog: %s (top=%s)\n", design_path.c_str(),
                  top.c_str());
      ord::dbVerilogNetwork vnet(sta);
      sta::VerilogReader reader(&vnet);
      ord::setDbNetworkLinkFunc(&vnet, &reader);
      reader.read(design_path.c_str());
      bool ok = ord::dbLinkDesign(top.c_str(), &vnet, db, logger,
                                  /*hierarchy=*/false);
      if (!ok) {
        std::fprintf(stderr, "error: link of top '%s' failed\n", top.c_str());
        return 1;
      }
      db->triggerPostReadDb();
    } else {
      std::printf("Reading design: %s\n", design_path.c_str());
      std::ifstream stream(design_path, std::ios::binary);
      if (!stream.good()) {
        std::fprintf(stderr, "error: cannot open odb file %s\n",
                     design_path.c_str());
        return 1;
      }
      db->read(stream);
    }

    // --- 4. either run correctness tests, or the analyze demo ---------------
    if (!test_mode.empty()) {
      int rc = lgista::runTest(sta, test_mode);
      std::printf("=== test '%s' -> %s ===\n", test_mode.c_str(),
                  rc == 0 ? "PASS" : "FAIL");
      return rc;
    }

    lgista::LocalGraphIsta timer(sta);
    timer.prepare();

    if (!single_inst.empty()) {
      std::printf("=== local timing for instance '%s' ===\n",
                  single_inst.c_str());
      timer.analyzeInstance(single_inst, /*print_report=*/true);
    } else {
      std::printf("=== local timing for all instances ===\n");
      size_t n = timer.analyzeAll(/*print_report=*/true);
      std::printf("built %zu local graphs\n", n);
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }

  return 0;
}
