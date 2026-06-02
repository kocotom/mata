#include <cstdlib>
#include <vector>
#include <cassert>

#include "mata/nfa/simulation.hh"
#include "mata/nfa/delta.hh"
#include "mata/nfa/nfa.hh"
#include "mata/utils/generation-sparse-set.hh"
#include "mata/utils/partition-relation-pair.hh"
#include "mata/utils/partition.hh"
#include "mata/utils/sparse-set.hh"

using Partition = mata::utils::Partition;
using Relation = mata::utils::Relation;
using PartitionRelationPair = mata::utils::PartitionRelationPair;

/***************************************
*                                      *
* SIMULATION AUXILIARY DATA STRUCTURES *
*                                      *
***************************************/

namespace mata::nfa::simulation::detail {

/***************************************
*                                      *
*             FLATTEN NFA              *
*                                      *
***************************************/

FlattenNfa::FlattenNfa(const Nfa& aut, size_t num_of_symbols) : num_of_symbols_(num_of_symbols) {
    num_of_states_ = aut.num_of_states();
    targets_.resize(num_of_symbols_);
    offsets_.resize(num_of_symbols_);
    state_hashers_.resize(num_of_symbols);
    rev_targets_.resize(num_of_symbols_);
    rev_offsets_.resize(num_of_symbols_);
    hashers_.resize(num_of_symbols_);
    out_symbols_.resize(num_of_states_);
    valid_symbols_.reserve(num_of_symbols_);
    std::vector<size_t> rev_counts(num_of_states_+1, 0);
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        offsets_[symbol].reserve(num_of_states_ + 1);
        rev_offsets_[symbol].resize(num_of_states_ + 1, 0);
        state_hashers_[symbol].resize(num_of_states_);
    }
    std::vector<size_t> sym_idx(num_of_symbols_, 0);
    for(size_t state = 0; state < num_of_states_; ++state) {
        for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
            offsets_[symbol].push_back(sym_idx[symbol]);
        }
        for(const SymbolPost& post : aut.delta[state]) {
            out_symbols_[state].push_back(post.symbol);
            for(size_t target : post.targets) {
                targets_[post.symbol].push_back(target);
                hashers_[post.symbol].add_transition(uint32_t(state), uint32_t(target));
                state_hashers_[post.symbol][target].add_transition(uint32_t(target), uint32_t(state));
                ++rev_offsets_[post.symbol][target];
                ++sym_idx[post.symbol];
            }
        }
    }
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        offsets_[symbol].push_back(sym_idx[symbol]);
        if(targets_[symbol].size()) { valid_symbols_.insert(symbol); }
    }
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        size_t cumsum = 0;
        std::swap(rev_counts, rev_offsets_[symbol]);
        for(size_t source = 0; source < num_of_states_; ++source) {
            cumsum += rev_counts[source];
            rev_offsets_[symbol][source + 1] = cumsum;
        }
        rev_targets_[symbol].resize(cumsum);
        for(size_t source = 0; source < num_of_states_; ++source) {
            for(size_t target : successors(source, symbol)) {
                rev_targets_[symbol][rev_offsets_[symbol][target + 1] - rev_counts[target]] = source;
                --rev_counts[target];
            }
        }
    }
    if(num_of_symbols_ == 1) { return; }
    size_t num_of_buckets = get_bucket_count();
    std::vector<std::vector<size_t>> buckets(num_of_buckets);
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        hash::HashTuple hash_tuple = hashers_[symbol].get_hashes();
        buckets[get_bucket_idx(hash_tuple.hash1, hash_tuple.hash2, num_of_buckets)].push_back(symbol);
    }
    utils::SparseSet<size_t> bucket_collision{};
    for(size_t bucket_idx = 0; bucket_idx < num_of_buckets; ++bucket_idx) {
        size_t bucket_size = buckets[bucket_idx].size();
        if(bucket_size > 1) {
            bucket_collision.clear();
            for(size_t i = 0; i < bucket_size; ++i) {
                for(size_t j = i + 1; j < bucket_size; ++j) {
                    if(bucket_collision[j]) { continue; }
                    size_t symbol1 = buckets[bucket_idx][i];
                    size_t symbol2 = buckets[bucket_idx][j];
                    hash::HashTuple hash_tuple1 = hashers_[symbol1].get_hashes();
                    hash::HashTuple hash_tuple2 = hashers_[symbol2].get_hashes();
                    if(hash_tuple1.hash1 != hash_tuple2.hash1 || hash_tuple1.hash2 != hash_tuple2.hash2) { continue; }
                    if(are_interchangeable(symbol1, symbol2)) {
                        bucket_collision.insert(j);
                        valid_symbols_.erase(symbol2);
                    }
                }
            }
        }
    }
}

bool FlattenNfa::are_interchangeable(size_t symbol1, size_t symbol2) {
    assert(symbol1 < num_of_symbols_ && symbol2 < num_of_symbols_ && "Nonexisting symbols occured.");
    if(symbol1 == symbol2) { return true; }
    if(targets_[symbol1].size() != targets_[symbol2].size()) { return false; }
    for(size_t source = 0; source < num_of_states_; ++source) {
        if(offsets_[symbol1][source] != offsets_[symbol2][source] ) { return false; }
    }
    size_t num_of_targets = targets_[symbol1].size();
    for(size_t target = 0; target < num_of_targets; ++target) {
        if(targets_[symbol1][target] != targets_[symbol2][target]) {
            return false;
        }
    }
    return true;
}

bool FlattenNfa::is_covered_by(size_t symbol1, size_t symbol2) {
    assert(symbol1 < num_of_symbols_ && symbol2 < num_of_symbols_ && "Nonexisting symbols occured.");
    if(symbol1 == symbol2) { return true; }
    for(size_t source = 0; source < num_of_states_; ++source) {
        size_t from1 = offsets_[symbol1][source];
        size_t to1 = offsets_[symbol1][source + 1];
        size_t from2 = offsets_[symbol2][source];
        size_t to2 = offsets_[symbol2][source + 1];
        auto goals = std::span(targets_[symbol1]).subspan(from1, to1 - from1);
        for(size_t candidate : std::span(targets_[symbol2]).subspan(from2, to2 - from2)) {
            if(goals.empty()) { break; }
            if(candidate > goals.front()) { return false; }
            if(goals.front() == candidate) { goals = goals.subspan(1); }
            
        }
        if(!goals.empty()) { return false; }
    }
    return true;
}

bool FlattenNfa::is_reversed_correctly(void) {
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        for(size_t source = 0; source < num_of_states_; ++source) {
            for(size_t target : successors(source, symbol)) {
                bool success = true;
                for(size_t predecessor : predecessors(target, symbol)) {
                    success = false;
                    if(predecessor == source) {
                        success = true;
                        break;
                    }
                }
                if(!success) {return false; }
            }
        }
        for(size_t target = 0; target < num_of_states_; ++target) {
            for(size_t source : predecessors(target, symbol)) {
                bool success = true;
                for(size_t successor : successors(source, symbol)) {
                    success = false;
                    if(successor == target) {
                        success = true;
                        break;
                    }
                }
                if(!success) {return false; }
            }
        }
    }
    return true;
}

/***************************************
*                                      *
*          MERGED TRANSITIONS          *
*                                      *
***************************************/


MergedTransitions::MergedTransitions(FlattenNfa& aut, utils::Partition& part, size_t symbol) {
    symbol_ = symbol;
    num_of_states_ = aut.num_of_states();
    blocks_.reserve(num_of_states_);
    nodes_.reserve(2 * num_of_states_ + 1);
    blocks_.resize(part.num_of_blocks());
    nodes_.resize(part.num_of_nodes());
    utils::SparseSet<size_t> targets{};
    utils::SparseSet<size_t> new_targets{};
    targets.reserve(num_of_states_);
    new_targets.reserve(num_of_states_);
    size_t max_ancestor_idx = 0;
    while(max_ancestor_idx < part.num_of_nodes() && !part.get_node(max_ancestor_idx).parent()) {
        std::queue<size_t> ancestors_worklist{};
        ancestors_worklist.push(max_ancestor_idx);
        while(!ancestors_worklist.empty()) {
            size_t current_ancestor = ancestors_worklist.front();
            ancestors_worklist.pop();
            size_t current_node_idx = current_ancestor;
            while(part.get_node(current_node_idx).left_child()) {
                current_node_idx = *part.get_node(current_node_idx).left_child();
            }
            new_targets.clear();
            targets.clear();
            size_t first_idx = data_.size();
            while(true){
                for(auto& block_item : part.get_node(current_node_idx)) {
                    for(size_t state : aut.predecessors(block_item.state(), symbol_)) {
                        if(!targets[state]) { 
                            targets.insert(state);
                            new_targets.insert(state);
                        }
                    }
                }
                for(size_t state : new_targets) {
                    data_.push_back(state);
                }
                new_targets.clear();
                nodes_[current_node_idx] = Block{first_idx, data_.size()};
                if(!part.get_node(current_node_idx).left_child()) {
                    blocks_[part.get_node(current_node_idx).first().block().idx()] = Block{first_idx, data_.size()};
                } else {
                    ancestors_worklist.push(*part.get_node(current_node_idx).right_child());
                }
                if(current_node_idx == current_ancestor) { break; }
                current_node_idx = *part.get_node(current_node_idx).parent();
            }
        }
        ++max_ancestor_idx;
    }
}

void MergedTransitions::split_block(FlattenNfa& aut, utils::Partition& part, size_t old_block_idx) {
    size_t old_first = blocks_[old_block_idx].first_;
    size_t old_last = blocks_[old_block_idx].last_;
    size_t new_block_idx = blocks_.size();
    blocks_.emplace_back(Block{data_.size(), data_.size()});
    size_t source = old_first;
    while(source < old_last) {
        bool old_found = false;
        bool new_found = false;
        for(size_t target : aut.successors(data_[source], symbol_)) {
            size_t target_block_idx = part.get_block_idx_of_state(target);
            if(target_block_idx == old_block_idx) {
                old_found = true;
            } else if(target_block_idx == new_block_idx) {
                new_found = true;
            }
        }
        if(new_found) {
            blocks_[new_block_idx].last_++;
            data_.push_back(data_[source]);
        }
        if(!old_found) {
            old_last--;
            blocks_[old_block_idx].last_--;
            std::swap(data_[old_last], data_[source]);
        } else {
            source++;
        }
    }
    if(part.get_node(part.num_of_nodes() - 1).first().block().idx() ==
        old_block_idx) {
        nodes_.push_back(blocks_[new_block_idx]);
        nodes_.push_back(blocks_[old_block_idx]);
    } else {
        nodes_.push_back(blocks_[old_block_idx]);
        nodes_.push_back(blocks_[new_block_idx]);
    }
}

/***************************************
*                                      *
*           SPARSE RELATION            *
*                                      *
***************************************/


SparseRelation::SparseRelation(size_t degree) : degree_(degree) {
    assert(degree > 0);
    for(size_t i = 0; i < degree; ++i) {
        pool.emplace_back(0, 0, 0, 4 * i + 1, 0, 0);
        pool.emplace_back(0, 0, 4 * i, 0, 0, 0);
        pool.emplace_back(0, 0, 0, 0, 0, 4 * i + 3);
        pool.emplace_back(0, 0, 0, 0, 4 * i + 2, 0);
        rows_first.push_back(4 * i);
        rows_last.push_back(4 * i + 1);
        cols_first.push_back(4 * i + 2);
        cols_last.push_back(4 * i + 3);
    }
}


void SparseRelation::extend_and_copy(size_t row, size_t col) {
    assert(row < degree_);
    assert(col < degree_);

    size_t new_row_first = pool.size();
    rows_first.emplace_back(pool.size());
    pool.emplace_back();

    size_t new_col_first = pool.size();
    cols_first.emplace_back(pool.size());
    pool.emplace_back();

    size_t new_row_last = pool.size();
    rows_last.emplace_back(pool.size());
    pool.emplace_back();

    size_t new_col_last = pool.size();
    cols_last.emplace_back(pool.size());
    pool.emplace_back();

    pool[new_row_first].right = new_row_last;
    pool[new_row_last].left = new_row_first;
    pool[new_col_first].down = new_col_last;
    pool[new_col_last].up = new_col_first;

    size_t prev_old_cell = rows_first[row];
    size_t next_old_cell = pool[prev_old_cell].right;
    size_t prev_new_cell = new_row_first;

    size_t old_row_last = rows_last[row];

    while(next_old_cell != old_row_last) {
        size_t current_col = pool[next_old_cell].col;
        size_t new_down = cols_last[current_col];
        size_t new_up = pool[new_down].up;
        size_t id = pool.size();
        pool.emplace_back(degree_, current_col, prev_new_cell, new_row_last, new_up, new_down);
        pool[prev_new_cell].right = id;
        pool[new_down].up = id;
        pool[new_up].down = id;
        pool[new_row_last].left = id;
        prev_new_cell = id;
        prev_old_cell = next_old_cell;
        next_old_cell = pool[next_old_cell].right;
    }

    prev_old_cell = cols_first[col];
    next_old_cell = pool[prev_old_cell].down;
    prev_new_cell = new_col_first;

    size_t old_col_last = cols_last[col];

    while(next_old_cell != old_col_last) {
        size_t current_row = pool[next_old_cell].row;
        size_t new_right = rows_last[current_row];
        size_t new_left = pool[new_right].left;
        size_t id = pool.size();
        pool.emplace_back(current_row, degree_, new_left, new_right, prev_new_cell, new_col_last);
        pool[prev_new_cell].down = id;
        pool[new_right].left = id;
        pool[new_left].right = id;
        pool[new_col_last].up = id;
        prev_new_cell = id;
        prev_old_cell = next_old_cell;
        next_old_cell = pool[next_old_cell].down;
    }

    auto last_diag = pool[pool[rows_last.back()].left];
    if(!(last_diag.row == degree_ && last_diag.col == degree_)) {
        size_t id = pool.size();
        pool.emplace_back(degree_, degree_, pool[new_row_last].left, new_row_last, pool[new_col_last].up, new_col_last);
        pool[pool[new_row_last].left].right = id;
        pool[pool[new_col_last].up].down = id;
        pool[new_row_last].left = id;
        pool[new_col_last].up = id;
    }
    ++degree_;
}

bool SparseRelation::erase(size_t row_idx, size_t col_idx) {
    if(row_idx >= degree_ || col_idx >= degree_) {
        return false;
    }
    auto row_range = row(row_idx);
    for(auto it = row_range.begin(); it != row_range.end(); ++it) {
        if(*it == col_idx) { 
            unlink(it.cell()); 
            return true;
        }
    }
    return false;
}

void SparseRelation::set(size_t row_idx, size_t col_idx) {
    size_t row_elem = rows_first[row_idx];
    size_t row_stop = rows_last[row_idx];

    while(pool[row_elem].right != row_stop) {
        row_elem = pool[row_elem].right;
        auto& current_elem = pool[row_elem];
        if(current_elem.col == col_idx) { return; }
        if(current_elem.col > col_idx) { break; }
    }

    size_t col_elem = cols_first[col_idx];
    size_t col_stop = cols_last[col_idx];

    while(pool[col_elem].down != col_stop) {
        auto& current_elem = pool[pool[col_elem].down];
        if(current_elem.row > row_idx) { break; }
        col_elem = current_elem.down;
    }

    pool.push_back(Cell(row_idx, col_idx, row_elem, pool[row_elem].right, col_elem, pool[col_elem].down));
    size_t new_idx = pool.size() - 1;
    auto& current_row_elem = pool[row_elem];
    auto& current_col_elem = pool[col_elem];
    pool[current_row_elem.right].left = new_idx;
    pool[current_col_elem.down].up = new_idx;
    current_row_elem.right = new_idx;
    current_col_elem.down = new_idx;
}

void SparseRelation::print_sr(void) {
    std::cout << "DEGREE: " << degree_ << std::endl << std::endl;
    std::cout << "ROWS: " << std::endl;

    for(size_t row_idx = 0; row_idx < degree_; ++row_idx) {
        auto& start_border = pool[rows_first[row_idx]];
        auto& end_border = pool[rows_last[row_idx]];
        
        std::cout << rows_first[row_idx] << " : " << start_border.right << " => || "; 

        auto current = start_border.right;
        auto last = rows_last[row_idx];
        while(current != last) {
            auto& c = pool[current];
            std::cout << "<" << current << " : (" << c.row << " " << c.col << ")[" << c.left << " " << c.right << " " << c.up << " " << c.down << "] "; 
            current = c.right;
        }

        std::cout << " || <= " << rows_last[row_idx] << " : " << end_border.left << std::endl;
    }
}

void SparseRelation::print_borders(void) {
    std::cout << "ROW START: ";
    for(size_t i = 0; i < degree_; ++i) {
        std::cout << rows_first[i] << " ";
    }
    std::cout << std::endl << "ROW END: ";
    for(size_t i = 0; i < degree_; ++i) {
        std::cout << rows_last[i] << " ";
    }
    std::cout << std::endl << "COL START: ";
    for(size_t i = 0; i < degree_; ++i) {
        std::cout << cols_first[i] << " ";
    }
    std::cout << std::endl << "COL END: ";
    for(size_t i = 0; i < degree_; ++i) {
        std::cout << cols_last[i] << " ";
    }
}

}

namespace mata::nfa::simulation{

/***************************************
*                                      *
*              SIMULATION              *
*                                      *
***************************************/


/** This function performs splitting of blocks with respect to a given
* marker idx. If the marker is omitted, the default marker 0 is used.
* The function calls split_update_data method which performs neccessary
* updates of various data structures.
* @brief performs a blocks splitting
* @param symbol given marker idx
*/
void Simulation::split_blocks() {
    // iterates through all blocks which contain at least one marked
    // state with respect to the given marker index
    std::vector<size_t> unprocessed_blocks{};
    for(size_t block_idx : part_.touched_blocks()) { 
        if(!part_.is_fully_marked(block_idx)) {
            unprocessed_blocks.emplace_back(block_idx);
        } 
    };
    for(size_t block_idx : unprocessed_blocks) {
        // performs splitting of a given block with respect to the
        // given marker. It also updates touched flags in context
        // of all markers
        std::vector<utils::SplitPair> split_pairs =
            part_.split_blocks_inner(block_idx);
        // if the block was split, data structures will be updated
        if(split_pairs.size()) {
            size_t idx_c = split_pairs.front().first_block_idx;
            size_t idx_d = split_pairs.front().second_block_idx;
            split_update_data(idx_c, idx_d);
        }
    }
}

mata::utils::PartitionRelationPair Simulation::relation_from_desymbolized_nfa(void) {
    std::vector<size_t> valid_blocks_idx{};
    size_t num_of_blocks = part_.num_of_blocks();
    for(size_t i = 0; i < num_of_blocks; ++i) {
        if(part_.get_block_repr_of_block_idx(i) < init_nfa_.num_of_states()) {
            valid_blocks_idx.push_back(i);
        }
    }

    std::vector<std::vector<size_t>> valid_blocks{};
    for(size_t i : valid_blocks_idx) {
        valid_blocks.emplace_back();
        for(auto &block_item : part_.get_block(i)) {
            valid_blocks.back().emplace_back(block_item.state());
        }
    }

    utils::Partition final_part{init_nfa_.num_of_states(), valid_blocks};
    utils::Relation final_rel{init_nfa_.num_of_states(), valid_blocks.size()};

    size_t num_of_final_blocks = valid_blocks.size();
    for(size_t i = 0; i < num_of_final_blocks; ++i) {
        for(size_t j = i; j < num_of_final_blocks; ++j) {
            size_t former_i_idx = part_.get_block_idx_of_state(final_part.get_block_repr_of_block_idx(i));
            size_t former_j_idx = part_.get_block_idx_of_state(final_part.get_block_repr_of_block_idx(j));
            if(sparse_rel_.get(former_i_idx, former_j_idx)) {
                final_rel.set(i, j, true);
            }
            if(sparse_rel_.get(former_j_idx, former_i_idx)) {
                final_rel.set(j, i, true);
            }
        }
    }

    return PartitionRelationPair(final_part, final_rel);
}

PartitionRelationPair Simulation::compute_in_desymbolized_simulation(void) {

   // automaton with no states produces a trivial preorder
    if(!init_nfa_.num_of_states()) { 
        return PartitionRelationPair(1);
    }

    size_t num_of_states{init_nfa_.num_of_states()};
    size_t num_of_symbols{0};
    size_t new_state = num_of_states;
    std::vector<std::vector<size_t>> symbol_classes{};
    Nfa desymbolized_nfa(num_of_states);
    for(size_t state = 0; state < num_of_states; ++state) {
        for(const SymbolPost& post : init_nfa_.delta[state]) {
            desymbolized_nfa.delta.add(state, 0, new_state);
            if(post.symbol >= num_of_symbols ) { 
                num_of_symbols = post.symbol + 1; 
                symbol_classes.resize(num_of_symbols);
            }
            symbol_classes[post.symbol].push_back(new_state);
            for(size_t target : post.targets) {
                desymbolized_nfa.delta.add(new_state, 0, target);
            }
            ++new_state;
        }
        if(init_nfa_.final[state]) { desymbolized_nfa.final.insert(state); }
    }

    if(num_of_symbols <= 1) { return compute_simulation(); }

    std::vector<std::vector<size_t>> used_symbol_classes{};
    for(auto i : symbol_classes) {
        if(!std::empty(i)) { used_symbol_classes.push_back(i); }
    }

    state_structures_reset(desymbolized_nfa);
    num_of_symbols_ = 1;

    //constructing simpler automata representation for fast search
    nfa_ = FlattenNfa(desymbolized_nfa, 1);
    size_t num_of_blocks = used_symbol_classes.size() + 1;

    part_ = Partition(new_state, used_symbol_classes);
    rel_ = Relation(new_state, num_of_blocks);
    sparse_rel_ = SparseRelation(num_of_blocks);
    for(size_t i = 0; i < num_of_blocks; ++i) { sparse_rel_.set(i, i); }

    sim_init();

    while(!refiner_1.empty()) {
        sim_update_data();
        size_t blocks;
        do{
        blocks = part_.num_of_blocks();
            split_first(); 
            split_second();
        } while(blocks != part_.num_of_blocks());
        refine();
    }

    return relation_from_desymbolized_nfa();
}

/** This function performs initial splitting of blocks with respect
* to their final/nonfinal property. After the split, no final state
* will share a block with a non-final state and vice versa.
* All necessary data structures are updated in place.
* @brief performs initial block splitting using non/final states condition
*/
void Simulation::split_by_final_states(void) {
    // if there is no final state or if all states are final, no
    // split is needed to be done
    if(num_of_states_ != init_nfa_.final.size() && init_nfa_.final.size()) {
        part_.clear_markings();
        size_t num_of_blocks = part_.num_of_blocks();

        // marking all final states
        for(const auto& i : init_nfa_.final) { part_.mark(i); }
        
        // iterating through blocks which contain at least one final state
        std::vector<size_t> touched_blocks;
        for(size_t idx : part_.touched_blocks()) {
            touched_blocks.push_back(idx);
        }

        for(size_t idx : touched_blocks) {

            // splitting
            std::vector<utils::SplitPair> split_pairs =
            part_.split_blocks_inner(idx);            
            
            // if the blocks were split, the simulation 
            // data structures have to be updated
            if(split_pairs.size()) {
                size_t idx_c = split_pairs.front().first_block_idx;
                size_t idx_d = split_pairs.front().second_block_idx;
                num_of_blocks = part_.num_of_blocks();
                sparse_rel_.extend_and_copy(idx_c, idx_c);
                pred_repr.emplace_back();
                in.emplace_back();
                in.emplace_back();

                if(part_.is_touched(idx_d)) {
                    for(size_t i = 0; i < num_of_blocks; ++i) {
                        if(!part_.is_touched(i)) {
                            sparse_rel_.erase(idx_d, i);
                        }
                    }
                } else {
                    for(size_t i = 0; i < num_of_blocks; ++i) {
                        if(part_.is_touched(i)) {
                            sparse_rel_.erase(i, idx_d);
                        }
                    }
                }

                // updating MergedTransitions
                for(size_t symbol : nfa_.get_valid_symbols()) {
                    mt_[symbol].split_block(nfa_, part_, idx_c);
                }
            }
        }
        
        part_.clear_markings();
    }
}

/** This function performs initialization of all data structures
* so the computation won't be corrupted if the Simulation object
* will be reused.
*/
void Simulation::state_structures_reset(const Nfa &nfa) {
    num_of_states_ = nfa.num_of_states();

    mt_.clear();
    in.clear();
    out.clear();
    repr_.clear();
    split_counts_.clear();
    pred_tmp.clear();
    used_symbols.clear();
    touched_1.clear();
    touched_2.clear();
    not_rel_1.clear();
    not_rel_2.clear();
    refiner_1.clear();
    refiner_2.clear();


    // limited by maximal number of nodes
    not_rel_1.reserve(2 * num_of_states_  + 1);
    not_rel_2.reserve(2 * num_of_states_  + 1);

    // limited by maximal number of blocks
    touched_1.extend_domain(2 * num_of_states_ + 1);
    touched_2.extend_domain(2 * num_of_states_ + 1);

    num_of_symbols_ = get_num_of_symbols();

    in = std::vector<std::vector<size_t>>{};

    // limited by maximal number of nodes
    refiner_1 = GenerationSparseSet<size_t>(2 * num_of_states_ + 1);
    refiner_2 = GenerationSparseSet<size_t>(2 * num_of_states_ + 1);

    // limited by maximal number of blocks
    pred_tmp = std::vector<size_t>{};
    pred_tmp.reserve(num_of_states_);

    pred_repr = std::vector<std::vector<std::vector<size_t>>>();
    pred_repr.resize(num_of_symbols_);

    noticed = GenerationSparseSet<size_t>(num_of_states_ + num_of_symbols_);

    smart_counters_ = SmartCounters(num_of_symbols_);

    // reusage of the the initial partition-relation pair
    part_ = init_part_;
    sparse_rel_ = init_sparse_rel_;
}


/** Key function which computes simulation over a non-deterministic
* automaton.
* @brief computes a simulation
* @return returns simulation over an automaton
*/
PartitionRelationPair Simulation::compute_simulation(void) {
    // automaton with no states produces a trivial preorder
    // TODO: Empty partition-relation set
    if(!init_nfa_.num_of_states()) { 
        return PartitionRelationPair(1);
    }

    state_structures_reset(init_nfa_);

    // constructing simpler automata representation for fast search
    nfa_ = FlattenNfa(init_nfa_, num_of_symbols_);
    sim_init();
    refiner_2.clear();
    while(!refiner_1.empty()) {
        sim_update_data();
        while(true) {
            size_t former_blocks = part_.num_of_blocks();
            if(split_first()) {
                while(true) {
                    size_t processed = split_first_correction(former_blocks);
                    if(part_.num_of_blocks() == processed) { break; }
                    former_blocks = processed;
                }
            }
            if(!split_second()) { break; }
        }
        refine();
    }

    size_t num_of_blocks = part_.num_of_blocks();
    rel_ = Relation(num_of_blocks, num_of_blocks);
    for(size_t i = 0; i < num_of_blocks; ++i) {
        auto row = sparse_rel_.row(i);
        for(auto j = row.begin(); j != row.end(); ++j) {
            rel_.set(i, *j, true);
        }
    }

    return PartitionRelationPair(part_, rel_);
}

/** Initialization of key data structures for computing
* simulation relation
* @brief initializes data structures
* @return returns simulation over an automaton
*/
void Simulation::sim_init(void) {
    
    // merged transition structure init
    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
        if(nfa_.get_valid_symbols()[symbol]) {
            mt_.push_back(MergedTransitions(nfa_, part_, symbol));
        } else { mt_.emplace_back(); }
    }
    
    split_by_final_states();

    size_t num_of_blocks = part_.num_of_blocks();
    size_t num_of_nodes = part_.num_of_nodes();
    
    for(size_t symbol : nfa_.get_valid_symbols()) {
        pred_repr[symbol].resize(num_of_blocks);
    }

    for(size_t row = 1; row < num_of_blocks; ++row) {
        smart_counters_.copy_row(0);
    }
    smart_counters_.resize_map(num_of_blocks);
    for(size_t i = 0; i < num_of_blocks; ++i) {
        size_t repr = part_.get_block_repr_of_block_idx(i);
        for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
            if(nfa_.get_valid_symbols()[symbol] && nfa_.has_successor(repr, symbol)) {
                smart_counters_.map_insert(i, symbol);
            }
        }
    }
    smart_counters_.introduce_keys();

    // initialization of block info
    for(size_t block_idx = 0; block_idx < num_of_blocks; ++block_idx) {
        split_counts_.push_back(0);
        
        // choosing a representative state
        repr_.push_back(part_.get_block_repr_of_block_idx(block_idx));
        
        // creating a list of representative states which
        // are predecessors of a processed block
        // TODO: aviod duplicate representative predecessors
        for(size_t symbol : nfa_.get_valid_symbols()) {
            noticed.clear();
            for(size_t state : nfa_.successors(repr_.back(), symbol)) {
                size_t block_idx = part_.get_block_idx_of_state(state);
                if(!noticed[block_idx]) {
                    pred_repr[symbol][block_idx].emplace_back(repr_.back());
                    noticed.insert(block_idx);
                }
                
            }
        }
        noticed.clear();
    }
    
    // initialization of node info
    for(size_t node_idx = 0; node_idx < num_of_nodes; ++node_idx) {
        in.emplace_back();
        not_rel_1.emplace_back();
        not_rel_2.emplace_back();
    }
    
    for(size_t symbol : nfa_.get_valid_symbols()) {
        for(size_t node_idx = 0; node_idx < num_of_nodes; ++node_idx) {
            if(!mt_[symbol].get_node(node_idx).empty()) {
                in[node_idx].push_back(symbol);
            }
        }
    }

    // collecting states which have at least one successor
    // with respect to a given symbol
    for(size_t symbol : nfa_.get_valid_symbols()) {
        
        for(State from = 0; from < num_of_states_; ++from) {
            if(!nfa_.successors(from, symbol).empty()) {
                part_.mark(from);
            }
        }
        // state with no successor cannot simulate state with a successor  
        split_blocks();
        
        num_of_blocks = part_.num_of_blocks();
        for(size_t idx_d = 0; idx_d < num_of_blocks; ++idx_d) {
            if(part_.touched_blocks()[idx_d]) { continue; }

            auto sparse_col = sparse_rel_.col(idx_d);
            for(auto it = sparse_col.begin(); it != sparse_col.end(); ) {
                if(part_.touched_blocks()[*it]) {
                    it = sparse_col.erase(it);
                } else {
                    ++it;
                }

            }
        }
        part_.clear_markings();
    }

    // initialization of counters
    for(size_t idx_c = 0; idx_c < num_of_blocks; ++idx_c) {
        size_t idx_c_node = part_.get_node_idx_of_block_idx(idx_c);
        for(size_t idx_d = 0; idx_d < num_of_blocks; ++idx_d) {
            for(size_t symbol : nfa_.get_valid_symbols()) {
                //counters_[symbol].set(idx_c, idx_d, 0);
                smart_counters_.set(symbol, idx_c, idx_d, 0);
            }
            // if d cannot simulate c, then Notrel(c, d)
            // c becomes a refiner node
            if(!sparse_rel_.get(idx_c, idx_d)) {
                size_t idx_d_node = part_.get_node_idx_of_block_idx(idx_d);
                not_rel_1[idx_c_node].emplace_back(idx_d_node);
                refiner_1.insert(idx_c_node);
            }
        }
    }
 
    noticed.clear();
    // updating counters
    for(size_t idx_c = 0; idx_c < num_of_blocks; ++idx_c) {
        for(const auto& symbol : in[part_.get_node_idx_of_block_idx(idx_c)]) {
            for(auto& to : mt_[symbol].get_block(idx_c)) {
                size_t idx_to = part_.get_block_idx_of_state(to);
                if(to != repr_[idx_to]) { continue; }
                noticed.insert(idx_to);
            }

            for(size_t idx_e : noticed) {
                for(size_t idx_d = 0; idx_d < num_of_blocks; ++idx_d) {
                    smart_counters_.inc(symbol, idx_d, idx_e);
                }
            }
            noticed.clear();
        }
    }

    noticed.clear();
}


/** Preparation for data structures updates
* after splitting one block into two blocks
* @brief data structures update
* @param idx_c index of the former block
* @param idx_d index of the new block
*/
void Simulation::split_update_data_prepare(size_t idx_c, size_t idx_d) {

    for(size_t symbol : nfa_.get_valid_symbols()) {
        mt_[symbol].split_block(nfa_, part_, idx_c);
        pred_repr[symbol].emplace_back();
    }
    split_counts_.emplace_back(0);
    in.emplace_back();
    in.emplace_back();
    sparse_rel_.extend_and_copy(idx_c, idx_c);
    repr_.push_back(repr_[idx_c]);

    smart_counters_.copy_row(idx_c);
    smart_counters_.resize_map(part_.num_of_blocks());
    for(size_t symbol : nfa_.get_valid_symbols()) {
        if(nfa_.has_successor(part_.get_block_repr_of_block_idx(idx_d), symbol)) {
            smart_counters_.map_insert(idx_d, symbol);
        }
    }
    smart_counters_.introduce_keys();
}

/** Updating data after splitting one block into two pieces
* @brief data structures update
* @param idx_c index of the former block
* @param idx_d index of the new block
*/
void Simulation::split_update_data(size_t idx_c, size_t idx_d) {
    size_t idx_x = idx_d;
    size_t idx_c_node = part_.get_node_idx_of_block_idx(idx_c);
    size_t idx_d_node = part_.get_node_idx_of_block_idx(idx_d);
    
    split_update_data_prepare(idx_c, idx_d);
    
    // updating representative predecessors after split
    for(size_t symbol : nfa_.get_valid_symbols()) {
        pred_tmp.clear();
        for(size_t state : pred_repr[symbol][idx_c]) {
            bool old_found = false;
            bool new_found = false;
            for(size_t s : nfa_.successors(state, symbol)) {
                if(part_.get_block_idx_of_state(s) == idx_c) {
                    old_found = true;
                } else if(part_.get_block_idx_of_state(s) == idx_d) {
                    new_found = true;
                }
                if(old_found && new_found) { break; }
            }
            if(old_found) { pred_tmp.push_back(state); }
            if(new_found) { pred_repr[symbol][idx_d].push_back(state); }
        }
        std::swap(pred_tmp, pred_repr[symbol][idx_c]);
        pred_tmp.clear();
    }
    noticed.clear();
    
    // updating counters after split
    for(size_t symbol : in[part_.get_node(idx_c_node).parent().value()]) {
        for(auto e : pred_repr[symbol][idx_c]) {
            touched_1.insert(part_.get_block_idx_of_state(e));
        }

        for(auto e : pred_repr[symbol][idx_d]) {
            size_t idx_e = part_.get_block_idx_of_state(e);
            if(touched_1[idx_e]) { touched_2.insert(idx_e); }
        }

        for(const size_t& idx_e : touched_2) {
            
            auto col_range = sparse_rel_.col(idx_c);
            for(auto it = col_range.begin(); it != col_range.end(); ++it) {
                smart_counters_.inc(symbol, *it, idx_e);
            } 
            
        }
        touched_1.clear();
        touched_2.clear();
    }
    
    repr_[idx_x] = part_.get_block_repr_of_block_idx(idx_x);
    touched_1.clear();
    noticed.clear();


    for(size_t symbol : nfa_.out_symbols(repr_[idx_x])) {
        if(!nfa_.get_valid_symbols()[symbol]) { continue; }
        for(const State& e : nfa_.successors(repr_[idx_x], symbol)) {
            touched_1.insert(part_.get_block_idx_of_state(e));
        }
        for(const size_t& idx_e : touched_1) {

            auto col_range = sparse_rel_.col(idx_e);
            for(auto it = col_range.begin(); it != col_range.end(); ++it) {
                smart_counters_.inc(symbol, *it, idx_x);
            }
            
        }
        touched_1.clear();
    }

    for(size_t symbol : in[part_.get_node(idx_c_node).parent().value()]) {
        if(!mt_[symbol].get_node(idx_c_node).empty()) {
            in[idx_c_node].push_back(symbol);
        }
        if(!mt_[symbol].get_node(idx_d_node).empty()) {
            in[idx_d_node].push_back(symbol);
        }
    }

    for(size_t symbol : nfa_.out_symbols(repr_.back())) {
        if(!nfa_.get_valid_symbols()[symbol]) { continue; }
        noticed.clear();
        for(size_t state : nfa_.successors(repr_.back(), symbol)) {
            size_t block_idx = part_.get_block_idx_of_state(state);
            if(!noticed[block_idx]) {
                pred_repr[symbol][block_idx].push_back(repr_.back());
                noticed.insert(block_idx);
            }
        }    
    }

    noticed.clear();

    not_rel_1.emplace_back();
    not_rel_1.emplace_back();
    not_rel_2.emplace_back();
    not_rel_2.emplace_back();
}

/** Determines whether a split of the first type is possible
* @brief checks the possibility of a split of the first type
* @param idx_b index of a processed node
* @param symbol processed symbol
* @return flag which says whether a split of the first type is possible
*/
bool Simulation::at_least_one_split(size_t idx_b, size_t symbol) {
    size_t block_b_idx = part_.get_node(idx_b).first().block().idx();
    for (State target : mt_[symbol].get_node(idx_b)) {

            if(!smart_counters_.get(symbol, block_b_idx, part_.get_block_idx_of_state(target))) {
                return true;
            }
    }

    return false;
}


/** Performs a split of the first type. Looks for a splitter
* transtitions of the first type
* @brief performs a split of the first type
*/
bool Simulation::split_first(void) {
    size_t num_of_blocks = part_.num_of_blocks();
    for(const size_t& idx_b : refiner_1) {
        const auto& block_b = part_.get_node(idx_b).first().block();
        for(size_t symbol : in[idx_b]) {
            // TODO: Are there any symbols related to the others which are possible to be skipped?
            if(!at_least_one_split(idx_b, symbol)) { continue; }        

            auto row_range = sparse_rel_.row(block_b.idx());
            for(auto it = row_range.begin(); it != row_range.end(); ++it) {
                for(const State& e : mt_[symbol].get_block(*it)) {
                    part_.mark(e);
                }
            }

            split_blocks();
            part_.clear_markings();
        }
    }
    return num_of_blocks != part_.num_of_blocks();
}

size_t Simulation::split_first_correction(size_t former_num_of_blocks) {
    size_t current_num_of_blocks = part_.num_of_blocks();
    for(size_t new_block = former_num_of_blocks; new_block < current_num_of_blocks; ++new_block) {
        for(const size_t& idx_b : refiner_1) {
            for(size_t symbol : in[idx_b]) {
                const auto& block_b = part_.get_node(idx_b).first().block();
                const size_t block_b_idx = block_b.idx();
                if(smart_counters_.get(symbol, block_b_idx, new_block)) { continue; }
                bool split = false;
                for(size_t state : mt_[symbol].get_node(idx_b)) {
                    if(!split && part_.get_block_idx_of_state(state) == new_block) {
                        split = true;
                    }
                    part_.mark(state);
                }
                if(!split) { 
                    part_.clear_markings();
                    continue;
                }
                auto row_range = sparse_rel_.row(block_b_idx);
                for(auto it = row_range.begin(); it != row_range.end(); ++it) {
                    if(part_.node_contains_block(idx_b, *it)) { continue; }
                    for(const State& e : mt_[symbol].get_block(*it)) {
                        part_.mark(e);
                    }
                }
            split_blocks();
            part_.clear_markings();
            }
        }
    }
    return current_num_of_blocks;
}

/** Performs a split of the second type. Looks for a splitter
* transtitions of the second type
* @brief performs a split of the second type
*/
bool Simulation::split_second(void) {
    size_t num_of_blocks = part_.num_of_blocks();
    for(const size_t& idx_a : refiner_1) {
        const auto& node_a = part_.get_node(idx_a);
        const auto& block_a = node_a.first().block();
        for(const auto& symbol : in[idx_a]) {

            for(size_t idx_b : part_.get_blocks_in_node(idx_a)) {
                for(size_t state : pred_repr[symbol][idx_b]) {
                    size_t idx_e = part_.get_block_idx_of_state(state);
                    noticed.insert(idx_e);
                }
                for(const size_t& noticed_idx : noticed) {
                    split_counts_[noticed_idx]++;
                    touched_1.insert(noticed_idx);
                }
                noticed.clear();
            }
            
            for(const size_t& idx_e : touched_1) {
                if(smart_counters_.get(symbol, block_a.idx(), idx_e) == split_counts_[idx_e]) {
                    touched_2.insert(idx_e);
                }

                split_counts_[idx_e] = 0;
            }
            for(const State& e : mt_[symbol].get_node(idx_a)) {
                if(touched_2[part_.get_block_idx_of_state(e)]) {
                    part_.mark(e);
                }
            }
            touched_1.clear();
            touched_2.clear();

            split_blocks();
            part_.clear_markings();

        }
    }
    return num_of_blocks != part_.num_of_blocks();
}

/** Refine first
* @brief refine first TODO
*/
void Simulation::refine_notice(size_t idx_a, size_t symbol) {
    const auto& block = part_.get_node(idx_a).first().block();
    size_t block_idx = block.idx();
    for(const size_t& idx_b : not_rel_1[idx_a]) {
        for(const State& e : mt_[symbol].get_node(idx_b)) {

                if(!smart_counters_.get(symbol, block_idx, part_.get_block_idx_of_state(e))) {
                    noticed.insert(part_.get_block_idx_of_state(e));
                }
        }
    }
}

/** Refine second
* @brief refine second TODO
* Is it possible to process the same block multiple times?
*/
void Simulation::refine_remove(size_t idx_a, size_t symbol) {
    touched_1.clear();
    for(const State& c : mt_[symbol].get_node(idx_a)) {
        size_t bl_idx_c = part_.get_block_idx_of_state(c);
        touched_1.insert(bl_idx_c);   
    }

    for(size_t bl_idx_c : touched_1) {
        auto row_range = sparse_rel_.row(bl_idx_c);
        for(auto it = row_range.begin(); it != row_range.end(); ) {
            if(!noticed[*it]) { 
                ++it;
                continue; 
            }
            size_t n_idx_c = part_.get_node_idx_of_block_idx(bl_idx_c);
            size_t n_idx_d = part_.get_node_idx_of_block_idx(*it);
            it = row_range.erase(it);
            not_rel_2[n_idx_c].emplace_back(n_idx_d);
            refiner_2.insert(n_idx_c);
        }
    }
    touched_1.clear();
}

/** Refine third
* @brief refine third TODO
*/
void Simulation::refine(void) {
    for(const size_t& idx_a : refiner_1) {
        for(const size_t& symbol : in[idx_a]) {
            refine_notice(idx_a, symbol);           
            refine_remove(idx_a, symbol);
            noticed.clear();    
        } 
        not_rel_1[idx_a].clear();
    }
    refiner_1.clear();

    std::swap(not_rel_1, not_rel_2);
    std::swap(refiner_1, refiner_2);
}

/** Sim update
* @brief sim update TODO
*/
void Simulation::sim_update_data(void) {
    for(const size_t& idx_b : refiner_1) {
        const auto& block_b = part_.get_node(idx_b).first().block();
        size_t block_b_idx = block_b.idx();
        for(const auto& symbol : in[idx_b]) {        
            for(const size_t& idx_e : not_rel_1[idx_b]) {
                size_t idx_block_e = part_.get_node(idx_e).first().block().idx();
                for(size_t state : pred_repr[symbol][idx_block_e]) {
                    smart_counters_.dec(symbol, block_b_idx,part_.get_block_idx_of_state(state));
                }
                noticed.clear();
            }
        }
    }      
}

}