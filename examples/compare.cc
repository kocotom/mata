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

bool sim_eq(PRP first, SimlibRel second) {
    if(second.size() != first.partition.num_of_states()) { return false; }
    size_t num_of_states = second.size();
    for(size_t state_row = 0; state_row < num_of_states; ++state_row) {
        for(size_t state_col = 0; state_col < num_of_states; ++state_col) {
            size_t block_idx_row = first.partition.get_block_idx_of_state(state_row);
            size_t block_idx_col = first.partition.get_block_idx_of_state(state_col);
            if(first.relation.get(block_idx_row, block_idx_col) !=
                second.get(state_row, state_col)) { return false; }
        }
    }
    return true;
}

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

    std::cout << "Filename: " << argv[1] << std::endl;

    using clock = std::chrono::high_resolution_clock;

    auto new_start = clock::now();
    auto new_res = NfaSim(aut).compute_simulation();
    auto new_end = clock::now();

    auto old_start = clock::now();
    auto old_res = compute_relation(aut);
    auto old_end = clock::now();

    bool correctness = sim_eq(new_res, old_res);
    std::cout << "Correctness: " << correctness << std::endl;

    std::cout << "New simulation: " << std::chrono::duration_cast<std::chrono::microseconds>(new_end - new_start).count() << std::endl;
    std::cout << "Old simulation: " << std::chrono::duration_cast<std::chrono::microseconds>(old_end - old_start).count() << std::endl;
    std::cout << std::endl;


    return EXIT_SUCCESS;
}
