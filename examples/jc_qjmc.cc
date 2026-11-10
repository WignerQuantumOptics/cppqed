// Copyright András Vukics 2024. Distributed under the Boost Software License, Version 1.0.
//
// Highest-level executable: Jaynes-Cummings QJMC simulation
// axis 0 = qubit, axis 1 = mode

#include "Composite.h"
#include "JaynesCummings.h"
#include "Mode.h"
#include "Qbit.h"

#include "QuantumJumpMonteCarlo.h"

#include "Pars.h"
#include "pcg_random.hpp"

using namespace cppqedutils;
using namespace quantumdata;
using namespace quantumtrajectory;
using namespace quantumdynamics;


int main(int argc, const char* argv[])
{
  // ─── Parameter registration ───────────────────────────────────────────
  // All subsystem and trajectory parameters registered into a single op

  auto op = optionParser();

  cppqedutils::trajectory::Pars<qjmc::Pars<qjmc::Algorithm::integrating, pcg64>> pt{op};
  qbit::Pars pq{op, "q"};
  mode::Pars pm{op, "m"};
  double g;
  add(op, parameters::_("g", "Jaynes-Cummings coupling strength", 1., g));

  parse(op, argc, argv);

  // ─── Initial state ────────────────────────────────────────────────────
  // rank-2 tensor product: psi_ij = psiQ_i * psiM_j

  auto psiQ = qbit::init(pq);
  auto psiM = mode::init(pm);

  StateVector<2> psi{ psiQ % psiM };

  // ─── Broadcasters ─────────────────────────────────────────────────────

  Broadcaster<2,0> bQ{psi.extents};  // qubit on axis 0
  Broadcaster<2,1> bM{psi.extents};  // mode on axis 1

  // ─── Unpack subsystem descriptors ─────────────────────────────────────

  auto [cutoffQ, haQ, liQ, evQ, descrQ] = qbit::make(pq);
  auto [cutoffM, haM, liM, evM, descrM] = mode::make(pm);

  // ─── Compose ──────────────────────────────────────────────────────────

  auto ha = hana::make_tuple(
    bQ(haQ),
    bM(haM),
    jaynescummings::hamiltonian<0,1>(pm.cutoff, g)
  );

  auto li = bQ(liQ);
  for (auto& l : bM(liM)) li.push_back(std::move(l));

  auto ev = hana::make_tuple( bQ(evQ), bM(evM) );

  auto descr = json::object{
    {"system", "Jaynes-Cummings"},
    {"qubit",  descrQ},
    {"mode",   descrM},
    {"g",      g}
  };

  // ─── Run ──────────────────────────────────────────────────────────────

  quantumtrajectory::run<qjmc::Algorithm::integrating, ODE_EngineBoost>( ha, li, ev, descr, std::move(psi), pt, trajectory::observerNoOp );

}
