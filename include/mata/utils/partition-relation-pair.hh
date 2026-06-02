/** @file partitionrelation_pair.hh
 *  @brief Definition of a partition-relation pair.
 *
 *  This file contains definition of partition-relation pair and operations
 *  which allow us to manipulate with it.
 *
 *  Description:
 *
 *  A partition-relation pair is a tuple (P, Rel). It is an efficient
 *  representation of a preorder/quasiorder R, which is a reflexive and
 *  transitive binary relation.
 *  In this context, we consider a carrier set S which contains all
 *  natural numbers from 0 to |S|-1. These numbers are called states.
 *  P is a partition of S which corresponds to an equivalence relation
 *  induced by the preorder R.
 *  Rel is a partial order over P.
 *  Thus, (P, Rel) corresponds to a preorder relation R over states S.
 *  
 *  This file provides implementation of a partition P and defines the
 *  ExtendableSquareMatrix class which can be used to represent 
 *  the binary relation Rel. 
 *  These classes can be combined to represent the preorder R.
 *
 *  @author Tomáš Kocourek
 */

#ifndef _PARTITION_RELATION_PAIR_HH_
#define _PARTITION_RELATION_PAIR_HH_

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cassert>
#include <memory>

#include "mata/utils/partition.hh"
#include "mata/utils/extendable-square-matrix.hh"

namespace mata::utils {

/************************************************************************
*
*
*
*                         PARTITION-RELATION PAIR
*
*
*
*************************************************************************/

using Relation = CascadeSquareMatrix<char>;

class PartitionRelationPair {
    public:
    
        Partition partition{};
        Relation relation{1, 1};
    
    private:
            
        Relation pairs_to_relation(
            size_t num_of_states, std::vector<size_t> pairs) {
            
            assert(!(pairs.size() % 2) && 
                "The pairs vector has to contain even number of indices");
            size_t num_of_pairs = pairs.size() / 2;
            size_t max_block = 0;
            for(size_t value : pairs) {
                max_block = value > max_block ? value : max_block;
            }
            assert(max_block < num_of_states &&
                "Index bigger than the number of states occured.");
            Relation rel = Relation(num_of_states, max_block+1);        
            for(size_t idx = 0; idx < num_of_pairs; ++idx) {
                rel.set(pairs[2 * idx], pairs[2 * idx + 1], 1);
            }
            return rel;
        }
        
    public:

        PartitionRelationPair(size_t num_of_states) {    
            partition = Partition(num_of_states);
            relation = Relation(num_of_states, 1);
        }
        
        PartitionRelationPair(
            Partition part, Relation rel) {
            assert(part.num_of_blocks() == rel.size());
            partition = std::move(part);
            relation = std::move(rel);
            }
            
        PartitionRelationPair(
            size_t num_of_states, StateBlocks& part,
            const std::vector<size_t>& pairs) {
            
            partition = Partition(num_of_states, part);
            relation = pairs_to_relation(num_of_states, pairs);
            assert(partition.num_of_blocks() == relation.size());   
        }
            
        PartitionRelationPair(const PartitionRelationPair &prp) {
            *this = prp;
        }
        
        
        PartitionRelationPair& operator=(const PartitionRelationPair &prp) {
            this->partition = prp.partition;
            this->relation = prp.relation;
            return *this;
        }

        size_t num_of_states(void) const {
            return partition.num_of_states();}
        
        size_t num_of_blocks(void) const {
            return partition.num_of_blocks();}
            
        bool related_states(State fst, State snd) const {
            assert(fst < partition.num_of_states());
            assert(snd < partition.num_of_states());    
            size_t fst_block_idx = partition.get_block_idx_of_state(fst);
            size_t snd_block_idx = partition.get_block_idx_of_state(snd);
            return relation.get(fst_block_idx, snd_block_idx);
        }
        
        bool related_blocks(size_t fst, size_t snd) const {
            assert(fst < partition.num_of_blocks());
            assert(snd < partition.num_of_blocks());    
            return relation.get(fst, snd);
        }

        void relate_blocks(size_t fst, size_t snd) {
            assert(fst < partition.num_of_blocks());
            assert(snd < partition.num_of_blocks());    
            relation.set(fst, snd, true);
        }
        
        void unrelate_blocks(size_t fst, size_t snd) {
            assert(fst < partition.num_of_blocks());
            assert(snd < partition.num_of_blocks());    
            relation.set(fst, snd, false);
        }
        
        std::vector<State> partition_block(State state) const {
            assert(state < partition.num_of_states());
            return partition.states_in_same_block(state);
        }
        
        bool is_preorder(void) const {
            return relation.is_reflexive() 
                && relation.is_antisymetric() 
                && relation.is_transitive();
        }
        
        // operators
        
        bool operator==(const PartitionRelationPair& other) const {
            if(partition != other.partition) { return false; }
            if(relation.size() != other.relation.size()) { return false; }
            size_t rel_size = relation.size();
            for(size_t row = 0; row < rel_size; ++row) {
                for(size_t col = 0; col < rel_size; ++col) {
                    State row_repr = partition.get_block(row).first().state();
                    State col_repr = partition.get_block(col).first().state();
                    size_t other_row = other.partition.get_block_idx_of_state(row_repr);
                    size_t other_col = other.partition.get_block_idx_of_state(col_repr);
                    if(relation.get(row, col) != 
                        other.relation.get(other_row, other_col)) { 
                        
                        return false;
                    }
                }
            }
            return true;
        }
        
        friend std::ostream& operator<<(
            std::ostream& os, const PartitionRelationPair& prp) {
        
            os << "PARTITION:" << std::endl << prp.partition << std::endl;
            os << "RELATION:" << std::endl << prp.relation
               << std::endl;
            return os;
        }
        
}; // PartitionRelationPair



}

#endif
