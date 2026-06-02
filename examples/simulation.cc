// example1.cc - constructing an automaton, then dumping it

#include "mata/nfa/nfa.hh"
#include "mata/nfa/simulation.hh"
#include <cmath>
#include <random>
#include "mata/utils/sparse-set.hh"
#include "mata/utils/partition-relation-pair.hh"
#include "mata/nfa/delta.hh"
#include "mata/nfa/nfa.hh"
#include "mata/nfa/strings.hh"
#include "mata/nfa/builder.hh"
#include "mata/nfa/plumbing.hh"
#include "mata/nfa/algorithms.hh"
#include "mata/parser/re2parser.hh"
#include <chrono>

#include <iostream>
#include <fstream>

using namespace mata::nfa;
using namespace mata::utils;
using namespace mata::parser;
using namespace mata::nfa::algorithms;
using SimlibRel = Simlib::Util::BinaryRelation;
using NfaSim = simulation::Simulation;
using PRP = PartitionRelationPair;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Input file missing\n";
        return EXIT_FAILURE;
    }

    std::string filename = argv[1];

    std::fstream fs(filename, std::ios::in);
    if (!fs) {
        std::cerr << "Could not open file \'" << filename << "'\n";
        return EXIT_FAILURE;
    }

    mata::parser::Parsed parsed;
    Nfa aut;
    try {
        parsed = mata::parser::parse_mf(fs, true);
        fs.close();

        std::vector<mata::IntermediateAut> inter_auts = mata::IntermediateAut::parse_from_mf(parsed);

        if (inter_auts[0].is_nfa())
            aut = mata::nfa::builder::construct(inter_auts[0]);
    }
    catch (const std::exception& ex) {
        fs.close();
        std::cerr << "libMATA error: " << ex.what() << "\n";
        return EXIT_FAILURE;
    }

    
    auto simulation = simulation::Simulation(aut).compute_simulation();
    std::cout << simulation;

    assert(simulation.related_states(0, 0));
    assert(simulation.related_blocks(0, 0));

    return EXIT_SUCCESS;
}
