// Generated MadGraph amplitude process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"

#include <utility>

namespace gra::amplitude {

// Compute all generated MG5 amplitude processes by value
std::vector<Process> AllProcesses() {
  return {
      Process{"DURHAM",
       "gg_ddbar",
       "d d~",
       "d d~",
       "d,d~",
       "g g > d d~ QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{1, -1},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ssbar",
       "s s~",
       "s s~",
       "s,s~",
       "g g > s s~ QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{3, -3},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ccbar",
       "c c~",
       "c c~",
       "c,c~",
       "g g > c c~ QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{4, -4},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_bbbar",
       "b b~",
       "b b~",
       "b,b~",
       "g g > b b~ QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{5, -5},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_uubar",
       "u u~",
       "u u~",
       "u,u~",
       "g g > u u~ QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{2, -2},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_uubarg",
       "u u~ g",
       "u u~ g",
       "u,u~,g",
       "g g > u u~ g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{2, -2, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ddbarg",
       "d d~ g",
       "d d~ g",
       "d,d~,g",
       "g g > d d~ g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{1, -1, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ssbarg",
       "s s~ g",
       "s s~ g",
       "s,s~,g",
       "g g > s s~ g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{3, -3, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ccbarg",
       "c c~ g",
       "c c~ g",
       "c,c~,g",
       "g g > c c~ g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{4, -4, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_bbbarg",
       "b b~ g",
       "b b~ g",
       "b,b~,g",
       "g g > b b~ g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{5, -5, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_gg",
       "gg",
       "g g",
       "g,g",
       "g g > g g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{21, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_ggg",
       "ggg",
       "g g g",
       "g,g,g",
       "g g > g g g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{21, 21, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"PHOTON",
       "yy_ll",
       "e+ e-",
       "e+ e-",
       "e+,e-",
       "a a > e+ e- QED=2 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"PHOTON",
       "yy_uubarg",
       "u u~ g",
       "u u~ g",
       "u,u~,g",
       "a a > u u~ g QED=2 QCD=1",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{2, -2, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"DURHAM",
       "gg_gggg",
       "gggg",
       "g g g g",
       "g,g,g,g",
       "g g > g g g g QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{21, 21, 21, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"PHOTON",
       "yy_uubargg",
       "u u~ g g",
       "u u~ g g",
       "u,u~,g,g",
       "a a > u u~ g g QED=2 QCD=2",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{2, -2, 21, 21},
       DecayStructure{DecayType::Full, true},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_epem",
       "e+ e-",
       "e+ e-",
       "e+,e-",
       "p p > e+ e- / h",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_mupmum",
       "mu+ mu-",
       "mu+ mu-",
       "mu+,mu-",
       "p p > mu+ mu- / h",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_taptam",
       "tau+ tau-",
       "tau+ tau-",
       "ta+,ta-",
       "p p > ta+ ta- / h",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_decay_epem",
       "e+ e-",
       "Z > {e+ e-}",
       "z(e+,e-)",
       "p p > z, z > e+ e-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-11, 11},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_decay_mupmum",
       "mu+ mu-",
       "Z > {mu+ mu-}",
       "z(mu+,mu-)",
       "p p > z, z > mu+ mu-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-13, 13},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_Z",
       "pp_z_decay_taptam",
       "tau+ tau-",
       "Z > {tau+ tau-}",
       "z(ta+,ta-)",
       "p p > z, z > ta+ ta-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-15, 15},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_ZJ",
       "pp_zj_epem",
       "e+ e- j",
       "Z > {e+ e-} j",
       "z(e+,e-),j",
       "p p > z j, z > e+ e-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 0},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::FlavourSet},
      Process{"MG5_PP_ZJ",
       "pp_zj_mupmum",
       "mu+ mu- j",
       "Z > {mu+ mu-} j",
       "z(mu+,mu-),j",
       "p p > z j, z > mu+ mu-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 0},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::FlavourSet},
      Process{"MG5_PP_ZJ",
       "pp_zj_taptam",
       "tau+ tau- j",
       "Z > {tau+ tau-} j",
       "z(ta+,ta-),j",
       "p p > z j, z > ta+ ta-",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 0},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::FlavourSet},
      Process{"MG5_PP_JJ",
       "pp_jj",
       "j j",
       "j j",
       "j,j",
       "p p > j j QED=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{21}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{0, 0},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::FlavourSet},
      Process{"MG5_PP_W",
       "pp_w_epve",
       "e+ ve",
       "e+ ve",
       "e+,ve",
       "p p > e+ ve",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 12},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_W",
       "pp_w_emvebar",
       "e- ve~",
       "e- ve~",
       "e-,ve~",
       "p p > e- ve~",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{11, -12},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_W",
       "pp_w_mupvm",
       "mu+ vm",
       "mu+ vm",
       "mu+,vm",
       "p p > mu+ vm",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 14},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_W",
       "pp_w_mumvmbar",
       "mu- vm~",
       "mu- vm~",
       "mu-,vm~",
       "p p > mu- vm~",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{13, -14},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_W",
       "pp_w_tapvt",
       "tau+ vt",
       "tau+ vt",
       "ta+,vt",
       "p p > ta+ vt",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 16},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_PP_W",
       "pp_w_tamvtbar",
       "tau- vt~",
       "tau- vt~",
       "ta-,vt~",
       "p p > ta- vt~",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{15, -16},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_JJ",
       "yy_jj",
       "j j",
       "j j",
       "j,j",
       "a a > j j QED=2 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}, AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{0, 0},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::FlavourSet},
      Process{"MG5_YY_WW",
       "yy_ww_epveemvebar",
       "e+ ve e- ve~",
       "W+ > {e+ ve} W- > {e- ve~}",
       "w+(e+,ve),w-(e-,ve~)",
       "a a > w+ w-, w+ > e+ ve, w- > e- ve~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-11, 12, 11, -12},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_epvemumvmbar",
       "e+ ve mu- vm~",
       "W+ > {e+ ve} W- > {mu- vm~}",
       "w+(e+,ve),w-(mu-,vm~)",
       "a a > w+ w-, w+ > e+ ve, w- > mu- vm~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-11, 12, 13, -14},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_epvetamvtbar",
       "e+ ve tau- vt~",
       "W+ > {e+ ve} W- > {tau- vt~}",
       "w+(e+,ve),w-(ta-,vt~)",
       "a a > w+ w-, w+ > e+ ve, w- > ta- vt~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{12}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-11, 12, 15, -16},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_mupvmemvebar",
       "mu+ vm e- ve~",
       "W+ > {mu+ vm} W- > {e- ve~}",
       "w+(mu+,vm),w-(e-,ve~)",
       "a a > w+ w-, w+ > mu+ vm, w- > e- ve~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-13, 14, 11, -12},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_mupvmmumvmbar",
       "mu+ vm mu- vm~",
       "W+ > {mu+ vm} W- > {mu- vm~}",
       "w+(mu+,vm),w-(mu-,vm~)",
       "a a > w+ w-, w+ > mu+ vm, w- > mu- vm~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-13, 14, 13, -14},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_mupvmtamvtbar",
       "mu+ vm tau- vt~",
       "W+ > {mu+ vm} W- > {tau- vt~}",
       "w+(mu+,vm),w-(ta-,vt~)",
       "a a > w+ w-, w+ > mu+ vm, w- > ta- vt~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{14}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-13, 14, 15, -16},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_tapvtemvebar",
       "tau+ vt e- ve~",
       "W+ > {tau+ vt} W- > {e- ve~}",
       "w+(ta+,vt),w-(e-,ve~)",
       "a a > w+ w-, w+ > ta+ vt, w- > e- ve~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-12}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-15, 16, 11, -12},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_tapvtmumvmbar",
       "tau+ vt mu- vm~",
       "W+ > {tau+ vt} W- > {mu- vm~}",
       "w+(ta+,vt),w-(mu-,vm~)",
       "a a > w+ w-, w+ > ta+ vt, w- > mu- vm~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-14}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-15, 16, 13, -14},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_WW",
       "yy_ww_tapvttamvtbar",
       "tau+ vt tau- vt~",
       "W+ > {tau+ vt} W- > {tau- vt~}",
       "w+(ta+,vt),w-(ta-,vt~)",
       "a a > w+ w-, w+ > ta+ vt, w- > ta- vt~ QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{16}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{-24}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-16}, std::vector<AmplitudeTopologyNode>{}}}}}},
       std::vector<int>{-15, 16, 15, -16},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_epemuubar",
       "e+ e- u u~",
       "Z > {e+ e-} u u~",
       "z(e+,e-),u,u~",
       "a a > z u u~, z > e+ e- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 2, -2},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_mupmumuubar",
       "mu+ mu- u u~",
       "Z > {mu+ mu-} u u~",
       "z(mu+,mu-),u,u~",
       "a a > z u u~, z > mu+ mu- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 2, -2},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_taptamuubar",
       "tau+ tau- u u~",
       "Z > {tau+ tau-} u u~",
       "z(ta+,ta-),u,u~",
       "a a > z u u~, z > ta+ ta- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{2}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-2}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 2, -2},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_epemddbar",
       "e+ e- d d~",
       "Z > {e+ e-} d d~",
       "z(e+,e-),d,d~",
       "a a > z d d~, z > e+ e- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 1, -1},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_mupmumddbar",
       "mu+ mu- d d~",
       "Z > {mu+ mu-} d d~",
       "z(mu+,mu-),d,d~",
       "a a > z d d~, z > mu+ mu- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 1, -1},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_taptamddbar",
       "tau+ tau- d d~",
       "Z > {tau+ tau-} d d~",
       "z(ta+,ta-),d,d~",
       "a a > z d d~, z > ta+ ta- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{1}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-1}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 1, -1},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_epemssbar",
       "e+ e- s s~",
       "Z > {e+ e-} s s~",
       "z(e+,e-),s,s~",
       "a a > z s s~, z > e+ e- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 3, -3},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_mupmumssbar",
       "mu+ mu- s s~",
       "Z > {mu+ mu-} s s~",
       "z(mu+,mu-),s,s~",
       "a a > z s s~, z > mu+ mu- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 3, -3},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_taptamssbar",
       "tau+ tau- s s~",
       "Z > {tau+ tau-} s s~",
       "z(ta+,ta-),s,s~",
       "a a > z s s~, z > ta+ ta- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{3}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-3}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 3, -3},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_epemccbar",
       "e+ e- c c~",
       "Z > {e+ e-} c c~",
       "z(e+,e-),c,c~",
       "a a > z c c~, z > e+ e- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 4, -4},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_mupmumccbar",
       "mu+ mu- c c~",
       "Z > {mu+ mu-} c c~",
       "z(mu+,mu-),c,c~",
       "a a > z c c~, z > mu+ mu- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 4, -4},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_taptamccbar",
       "tau+ tau- c c~",
       "Z > {tau+ tau-} c c~",
       "z(ta+,ta-),c,c~",
       "a a > z c c~, z > ta+ ta- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{4}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-4}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 4, -4},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_epembbbar",
       "e+ e- b b~",
       "Z > {e+ e-} b b~",
       "z(e+,e-),b,b~",
       "a a > z b b~, z > e+ e- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-11}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{11}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-11, 11, 5, -5},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_mupmumbbbar",
       "mu+ mu- b b~",
       "Z > {mu+ mu-} b b~",
       "z(mu+,mu-),b,b~",
       "a a > z b b~, z > mu+ mu- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-13}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{13}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-13, 13, 5, -5},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact},
      Process{"MG5_YY_ZJJ",
       "yy_zjj_taptambbbar",
       "tau+ tau- b b~",
       "Z > {tau+ tau-} b b~",
       "z(ta+,ta-),b,b~",
       "a a > z b b~, z > ta+ ta- QED=4 QCD=0",
       AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}},
       std::vector<AmplitudeTopology>{AmplitudeTopology{AmplitudeTopologyNode{std::vector<int>{23}, std::vector<AmplitudeTopologyNode>{AmplitudeTopologyNode{std::vector<int>{-15}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{15}, std::vector<AmplitudeTopologyNode>{}}}}, AmplitudeTopologyNode{std::vector<int>{5}, std::vector<AmplitudeTopologyNode>{}}, AmplitudeTopologyNode{std::vector<int>{-5}, std::vector<AmplitudeTopologyNode>{}}}},
       std::vector<int>{-15, 15, 5, -5},
       DecayStructure{DecayType::Full, false},
       MatrixElementForm::Generated,
       TopologyMode::Exact}
  };
}

// Compute processes in one generated amplitude family by value
std::vector<Process> Processes(
    const std::string &process_family) {
  std::vector<Process> selected;
  for (auto process : AllProcesses()) {
    if (process.process_family == process_family) {
      selected.push_back(std::move(process));
    }
  }
  return selected;
}

// Compute one exact generated process as a vector
std::vector<Process> Processes(
    const std::string &process_family, const std::string &process_name) {
  const auto process =
      FindProcess(process_family, process_name);
  if (!process.has_value()) { return {}; }
  return {*process};
}

// Find one generated process by its family and process names
std::optional<Process> FindProcess(
    const std::string &process_family, const std::string &process_name) {
  for (auto process : AllProcesses()) {
    if (process.process_family == process_family &&
        process.process_name == process_name) {
      return process;
    }
  }
  return std::nullopt;
}

// Compute the parameter card owned by one exact generated process
std::optional<std::string> ParameterCard(
    const std::string &process_family, const std::string &process_name) {
  if (process_family == "DURHAM" &&
      process_name == "gg_ddbar") {
    return "MG5cards/Durham/gg_ddbar/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ssbar") {
    return "MG5cards/Durham/gg_ssbar/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ccbar") {
    return "MG5cards/Durham/gg_ccbar/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_bbbar") {
    return "MG5cards/Durham/gg_bbbar/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_uubar") {
    return "MG5cards/Durham/gg_uubar/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_uubarg") {
    return "MG5cards/Durham/gg_uubarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ddbarg") {
    return "MG5cards/Durham/gg_ddbarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ssbarg") {
    return "MG5cards/Durham/gg_ssbarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ccbarg") {
    return "MG5cards/Durham/gg_ccbarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_bbbarg") {
    return "MG5cards/Durham/gg_bbbarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_gg") {
    return "MG5cards/Durham/gg_gg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_ggg") {
    return "MG5cards/Durham/gg_ggg/param_card.dat";
  }
  if (process_family == "PHOTON" &&
      process_name == "yy_ll") {
    return "MG5cards/Photon/yy_ll/param_card.dat";
  }
  if (process_family == "PHOTON" &&
      process_name == "yy_uubarg") {
    return "MG5cards/Photon/yy_uubarg/param_card.dat";
  }
  if (process_family == "DURHAM" &&
      process_name == "gg_gggg") {
    return "MG5cards/Durham/gg_gggg/param_card.dat";
  }
  if (process_family == "PHOTON" &&
      process_name == "yy_uubargg") {
    return "MG5cards/Photon/yy_uubargg/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_epem") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_mupmum") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_taptam") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_decay_epem") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_decay_mupmum") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_Z" &&
      process_name == "pp_z_decay_taptam") {
    return "MG5cards/Parton/MG5_PP_Z/param_card.dat";
  }
  if (process_family == "MG5_PP_ZJ" &&
      process_name == "pp_zj_epem") {
    return "MG5cards/Parton/MG5_PP_ZJ/param_card.dat";
  }
  if (process_family == "MG5_PP_ZJ" &&
      process_name == "pp_zj_mupmum") {
    return "MG5cards/Parton/MG5_PP_ZJ/param_card.dat";
  }
  if (process_family == "MG5_PP_ZJ" &&
      process_name == "pp_zj_taptam") {
    return "MG5cards/Parton/MG5_PP_ZJ/param_card.dat";
  }
  if (process_family == "MG5_PP_JJ" &&
      process_name == "pp_jj") {
    return "MG5cards/Parton/MG5_PP_JJ/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_epve") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_emvebar") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_mupvm") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_mumvmbar") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_tapvt") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_PP_W" &&
      process_name == "pp_w_tamvtbar") {
    return "MG5cards/Parton/MG5_PP_W/param_card.dat";
  }
  if (process_family == "MG5_YY_JJ" &&
      process_name == "yy_jj") {
    return "MG5cards/Photon/MG5_YY_JJ/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_epveemvebar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_epvemumvmbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_epvetamvtbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_mupvmemvebar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_mupvmmumvmbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_mupvmtamvtbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_tapvtemvebar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_tapvtmumvmbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_WW" &&
      process_name == "yy_ww_tapvttamvtbar") {
    return "MG5cards/Photon/MG5_YY_WW/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_epemuubar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_mupmumuubar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_taptamuubar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_epemddbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_mupmumddbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_taptamddbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_epemssbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_mupmumssbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_taptamssbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_epemccbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_mupmumccbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_taptamccbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_epembbbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_mupmumbbbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  if (process_family == "MG5_YY_ZJJ" &&
      process_name == "yy_zjj_taptambbbar") {
    return "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat";
  }
  return std::nullopt;
}

// Find one generated process by typed topology and ordered stable PDGs
std::optional<Process> FindProcess(
    const std::string &process_family, const AmplitudeTopology &topology,
    const std::vector<int> &stable_pdgs) {
  for (auto process : AllProcesses()) {
    if (process.process_family == process_family &&
        process.topology == topology &&
        process.stable_pdgs == stable_pdgs) {
      return process;
    }
  }
  return std::nullopt;
}

}  // namespace gra::amplitude
