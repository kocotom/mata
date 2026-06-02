#ifndef MATA_NFA_SIMULATION_HH
#define MATA_NFA_SIMULATION_HH

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <vector>
#include <cassert>
#include <queue>

#include "mata/nfa/delta.hh"
#include "mata/utils/partition.hh"
#include "mata/utils/sparse-set.hh"
#include "mata/utils/generation-sparse-set.hh"
#include "mata/utils/partition-relation-pair.hh"
#include "mata/nfa/nfa.hh"

namespace mata::nfa::simulation::detail {

    namespace hash {
        static constexpr uint64_t ADDITIVE_CONST = 0x9e3779b97f4a7c15ull;
        static constexpr uint64_t MULTIPLICATIVE_CONST_ONE  = 0xbf58476d1ce4e5b9ull;
        static constexpr uint64_t MULTIPLICATIVE_CONST_TWO  = 0x94d049bb133111ebull;

        static constexpr unsigned SHIFT_ONE = 30u;
        static constexpr unsigned SHIFT_TWO = 27u;
        static constexpr unsigned SHIFT_THREE = 31u;
        static constexpr unsigned SHIFT_FOUR = 32u;

        static constexpr uint64_t DEFAULT_ACC_ONE = 0x243f6a8885a308d3ull;
        static constexpr uint64_t DEFAULT_ACC_TWO = 0x9e3779b97f4a7c15ull;

        static constexpr uint64_t DEFAULT_SEED_ONE = 0x1bd11bdaa9fc1a22ull;
        static constexpr uint64_t DEFAULT_SEED_TWO = 0x57fefc2d2d0a3f4bull;

        constexpr uint64_t splitmix64(uint64_t val) noexcept {
            val += ADDITIVE_CONST;
            val = (val ^ (val >> SHIFT_ONE)) * MULTIPLICATIVE_CONST_ONE;
            val = (val ^ (val >> SHIFT_TWO)) * MULTIPLICATIVE_CONST_TWO;
            return val ^ (val >> SHIFT_THREE);
        }

        constexpr uint64_t encode_transition(uint32_t source, uint32_t target, uint64_t seed) noexcept {
            return splitmix64((uint64_t(source) << SHIFT_FOUR) ^ uint64_t(target) ^ seed);
        }

        struct HashTuple { std::uint64_t hash1, hash2; };

        class NfaHasher {
            private:
                uint64_t acc1_ = DEFAULT_ACC_ONE;
                uint64_t acc2_ = DEFAULT_ACC_TWO;
                uint64_t seed1_ = DEFAULT_SEED_ONE;
                uint64_t seed2_ = DEFAULT_SEED_TWO;
                uint32_t num_of_transitions_ = 0;

            public:
                constexpr NfaHasher(uint64_t seed1 = DEFAULT_SEED_ONE, uint64_t seed2 = DEFAULT_SEED_TWO) noexcept :
                    seed1_(seed1), seed2_(seed2) {}
                
                constexpr void reset() noexcept {
                    acc1_ = DEFAULT_ACC_ONE;
                    acc2_ = DEFAULT_ACC_TWO;
                    num_of_transitions_ = 0;
                }

                inline void add_transition(uint32_t source, uint32_t target) noexcept {
                    ++num_of_transitions_;
                    const uint64_t code1 = encode_transition(source, target, seed1_);
                    const uint64_t code2 = encode_transition(source, target, seed2_);
                    acc1_ ^= code1;
                    acc2_ ^= code2;
                }

                constexpr HashTuple get_hashes(void) const noexcept {
                    return {
                        acc1_ ^ num_of_transitions_ * MULTIPLICATIVE_CONST_ONE,
                        acc2_ ^ num_of_transitions_ * MULTIPLICATIVE_CONST_TWO
                    };
                }
        };

        class StateTransitionsHasher {
            private:
                uint64_t hash_ = 0;
                uint64_t seed_ = DEFAULT_SEED_ONE;

            public:
                constexpr StateTransitionsHasher(uint64_t seed = DEFAULT_SEED_ONE) noexcept : seed_(seed) {}

                constexpr void reset() noexcept {
                    hash_ = 0;
                    seed_ = DEFAULT_SEED_ONE;
                }

                inline void add_transition(uint32_t source, uint32_t target) noexcept {
                    hash_ ^= encode_transition(source, target, seed_);
                    hash_ |= 1ULL << 63;
                }

                inline void set_artificial_value(uint32_t value) { hash_ = value; }

                constexpr uint64_t get_hash(void) const noexcept { return hash_; }

        };

    }

    class CompactIndexMap{
    
        private:

            static constexpr int32_t empty_ = -1;

            std::vector<int32_t> map_{};
            int32_t num_of_keys_{0};

        public:

            CompactIndexMap() = default;
            void extend(size_t num_of_elements) {
                assert(num_of_elements <= map_.max_size() - map_.size());
                map_.resize(map_.size() + num_of_elements, empty_);
            }

            void add_element(size_t element) {
                if(element >= map_.size()) { map_.resize(element + 1, empty_); };
                if(map_[element] != empty_) { return; }
                assert(num_of_keys_ != std::numeric_limits<int32_t>::max());
                map_[element] = num_of_keys_++;            
            }

            [[nodiscard]] inline int32_t get_safe(size_t idx) const noexcept {
                if(idx >= map_.size()) { return empty_; }
                return map_[idx];
            }

            [[nodiscard]] inline int32_t get(size_t idx) const noexcept { 
                assert(idx < map_.size());
                return map_[idx]; 
            }

            bool has_key(size_t idx) const noexcept {
                if(idx >= map_.size()) { return false; }
                return map_[idx] != empty_;
            }

            inline int32_t get_num_of_keys(void) const noexcept { return num_of_keys_; }
            inline int32_t empty(void) const noexcept { return empty_; }
            inline size_t size(void) const noexcept { return map_.size(); }
    };


    class CompactTable{

        private:

            static constexpr int32_t row_size = 32;

            typedef struct SharedVec{
                std::vector<size_t> data_{};
                size_t participants_{0};

                SharedVec(size_t participants = 1) : participants_(participants) {
                    data_.resize(row_size, 0);
                }

            }SharedVector;

            
            std::vector<std::vector<size_t>> rows_{};
            std::vector<SharedVector> data_pool_{};
            size_t used_keys_{0};

        public:

            CompactTable() {
                rows_.emplace_back();
                data_pool_.emplace_back(1);
                rows_[0].emplace_back(0);
            }

            CompactTable(size_t num_of_rows, size_t num_of_keys) {
                rows_.resize(num_of_rows);

                size_t num_of_vectors = (num_of_keys + row_size - 1) / row_size;
                data_pool_.emplace_back(num_of_rows * num_of_vectors);

                for(size_t row = 0; row < num_of_rows; ++row) { 
                    rows_[row].resize(num_of_vectors, 0);
                }

                used_keys_ = num_of_keys;
            }

            size_t get(size_t row_idx, size_t key) const noexcept {
                assert(row_idx < rows_.size());
                return data_pool_[rows_[row_idx][key / row_size]].data_[key % row_size];
            }

            size_t num_of_used_keys(void) const noexcept { return used_keys_; }
            size_t num_of_rows(void) const noexcept { return rows_.size(); }
            size_t num_of_vecs(void) const noexcept { return data_pool_.size(); }
            double ratio_of_saved_cells(void) {
                size_t saved = 0;
                size_t all = 0;
                size_t num_of_pools = data_pool_.size();
                for(size_t i = 0; i < num_of_pools; ++i) {
                    saved += (data_pool_[i].participants_ - 1) * row_size;
                    all += data_pool_[i].participants_ * row_size;
                }
                return double(saved) / double(all);
            }


        private:

            size_t separate_vec(size_t vec_idx) {
                assert(vec_idx < data_pool_.size());
                if(data_pool_[vec_idx].participants_ == 1) { return vec_idx; }
                data_pool_.push_back(data_pool_[vec_idx]);
                --data_pool_[vec_idx].participants_;
                data_pool_[data_pool_.size() - 1].participants_ = 1;
                return data_pool_.size() - 1;
            }

        public:

            void set(size_t row_idx, size_t key, size_t value) noexcept {
                assert(row_idx < rows_.size());
                if(value == get(row_idx, key)) { return; }
                size_t vec_idx = rows_[row_idx][key / row_size];
                size_t new_vec_idx = separate_vec(vec_idx);
                if(vec_idx != new_vec_idx) {
                    rows_[row_idx][key / row_size] = new_vec_idx;
                }
                data_pool_[new_vec_idx].data_[key % row_size] = value;
            }

            void inc(size_t row_idx, size_t key) noexcept {
                assert(row_idx < rows_.size());
                size_t vec_idx = rows_[row_idx][key / row_size];
                size_t new_vec_idx = separate_vec(vec_idx);
                if(vec_idx != new_vec_idx) {
                    rows_[row_idx][key / row_size] = new_vec_idx;
                }
                ++data_pool_[new_vec_idx].data_[key % row_size];
            }

            void dec(size_t row_idx, size_t key) noexcept {
                assert(row_idx < rows_.size());
                size_t vec_idx = rows_[row_idx][key / row_size];
                size_t new_vec_idx = separate_vec(vec_idx);
                if(vec_idx != new_vec_idx) {
                    rows_[row_idx][key / row_size] = new_vec_idx;
                }
                --data_pool_[new_vec_idx].data_[key % row_size];
            }

            void copy_row(size_t row_idx) noexcept {
                assert(row_idx < rows_.size());
                rows_.push_back(rows_[row_idx]);
                size_t num_of_vecs = rows_.back().size();
                for(size_t i = 0; i < num_of_vecs; ++i) {
                    ++data_pool_[rows_.back()[i]].participants_;
                }
            }

            void introduce_keys(size_t num_of_keys = 1) noexcept {
                assert(used_keys_ <= std::numeric_limits<size_t>::max() - num_of_keys);
                if(num_of_keys < 1) { return; }
                const size_t old_vecs = (used_keys_ + row_size - 1) / row_size;
                used_keys_ += num_of_keys;
                const size_t new_vecs = (used_keys_ + row_size - 1) / row_size;
                size_t num_of_added_vecs = new_vecs - old_vecs;
                if(num_of_added_vecs) {
                    size_t num_of_rows = rows_.size();
                    size_t new_vec_idx = data_pool_.size();
                    data_pool_.emplace_back(num_of_added_vecs * num_of_rows);                
                    for(size_t row = 0; row < num_of_rows; ++row) {
                        rows_[row].resize(rows_[row].size() + num_of_added_vecs, new_vec_idx);
                    }
                }
            }
    };


    class SmartCounters {
        private:
            CompactTable data_{};
            CompactIndexMap map_{};

            size_t num_of_keys_{0};
            size_t num_of_active_keys_{0};
            size_t num_of_symbols_{0};

        public: 
            
            SmartCounters() = default;
            SmartCounters(size_t num_of_symbols) : num_of_symbols_(num_of_symbols) {};

            inline size_t get_num_of_keys(void) const noexcept { return num_of_keys_; }
            inline size_t get_num_of_active_keys(void) const noexcept { return num_of_active_keys_; }
            inline size_t get_num_of_rows(void) const noexcept { return data_.num_of_rows(); }

            inline void resize_map(size_t num_of_blocks) {
                map_.extend(num_of_blocks * num_of_symbols_ - map_.size());
            }

            inline void map_insert(size_t num_of_block, size_t symbol) {
                size_t key = num_of_block * num_of_symbols_ + symbol;
                if(map_.has_key(key)) { return; }
                map_.add_element(key);
                ++num_of_keys_;
            }

            inline void introduce_keys(void) {
                data_.introduce_keys(num_of_keys_ - num_of_active_keys_);
                num_of_active_keys_ = num_of_keys_;
            }

            inline void copy_row(size_t row_idx) { data_.copy_row(row_idx); }

            size_t get(size_t symbol, size_t row_idx, size_t col_idx) {
                size_t key = col_idx * num_of_symbols_ + symbol; 
                if(!map_.has_key(key)) { return 0; }
                return data_.get(row_idx, size_t(map_.get(key)));
            }

            void set(size_t symbol, size_t row_idx, size_t col_idx, size_t value) {
                size_t key = col_idx * num_of_symbols_ + symbol; 
                if(!map_.has_key(key)) { return; }
                data_.set(row_idx, size_t(map_.get(key)), value);
            }

            void inc(size_t symbol, size_t row_idx, size_t col_idx) {
                size_t key = col_idx * num_of_symbols_ + symbol; 
                if(!map_.has_key(key)) { return; }
                data_.inc(row_idx, size_t(map_.get(key)));
            }

            void dec(size_t symbol, size_t row_idx, size_t col_idx) {
                size_t key = col_idx * num_of_symbols_ + symbol; 
                if(!map_.has_key(key)) { return; }
                data_.dec(row_idx, size_t(map_.get(key)));
            }

            double ratio_of_saved_cells(void) {
                return data_.ratio_of_saved_cells();
            }
    };
    


    class SparseRelation{
        public:

            class Cell{
                public:

                    size_t row{0};
                    size_t col{0};
                    
                    size_t left{0};
                    size_t right{0};
                    size_t up{0};
                    size_t down{0}; 
                    
                    Cell() = default;
                    Cell(size_t row_, size_t col_, size_t left_, size_t right_, size_t up_, size_t down_) :
                    row(row_), col(col_), left(left_), right(right_), up(up_), down(down_) { }
                };


        private:

            std::vector<size_t> rows_first{};
            std::vector<size_t> cols_first{};
            std::vector<size_t> rows_last{};
            std::vector<size_t> cols_last{};
            std::vector<Cell> pool{};

            size_t degree_{};

        public:

            SparseRelation(size_t degree = 1);
            void extend_and_copy(size_t row, size_t col);

            class RowIterator {
                private:
                
                    friend class SparseRelation;
                    SparseRelation *rel_{};
                    size_t current_cell_{};

                public:

                    RowIterator(SparseRelation* rel, size_t current_cell) : rel_(rel), current_cell_(current_cell) {}

                    using iterator_category = std::forward_iterator_tag;
                    using value_type = size_t;
                    using difference_type = std::ptrdiff_t;

                    size_t operator*() const { return rel_->pool[current_cell_].col; }

                    RowIterator& operator++() { 
                        current_cell_ = rel_->pool[current_cell_].right; 
                        return *this; 
                    }

                    bool operator==(const RowIterator& other) const { return current_cell_ == other.current_cell_; }
                    bool operator!=(const RowIterator& other) const { return current_cell_ != other.current_cell_; }

                    size_t cell() const noexcept { return current_cell_; }
            };

            class RowRange {
                private:

                    friend class SparseRelation;
                    SparseRelation *rel_{};
                    size_t row_{0};
                    size_t row_start_{0};
                    size_t row_end_{0};

                    RowRange(SparseRelation *rel, size_t row) : rel_(rel), row_(row),
                        row_start_(rel->rows_first[row]), row_end_(rel->rows_last[row]) {}

                public:

                    RowIterator begin() { return RowIterator(rel_, rel_->pool[row_start_].right); }
                    RowIterator end() { return RowIterator(rel_, row_end_); }      

                    RowIterator erase(RowIterator it) {
                        size_t current_cell = it.current_cell_;
                        if(current_cell == row_end_) { return it; }
                        size_t next = rel_->pool[current_cell].right;
                        rel_->unlink(current_cell);
                        return RowIterator(rel_, next);
                    }

                    size_t row(void) const noexcept { return row_; }
            };

            class ColIterator {
                private:
                
                    friend class SparseRelation;
                    SparseRelation *rel_{};
                    size_t current_cell_{};

                public:

                    ColIterator(SparseRelation* rel, size_t current_cell) : rel_(rel), current_cell_(current_cell) {}

                    using iterator_category = std::forward_iterator_tag;
                    using value_type = size_t;
                    using difference_type = std::ptrdiff_t;

                    size_t operator*() const { return rel_->pool[current_cell_].row; }

                    ColIterator& operator++() { current_cell_ = rel_->pool[current_cell_].down; return *this; }

                    bool operator==(const ColIterator& other) const { return current_cell_ == other.current_cell_; }
                    bool operator!=(const ColIterator& other) const { return current_cell_ != other.current_cell_; }

                    size_t cell() const noexcept { return current_cell_; }
            };

            class ColRange {
                private:

                    friend class SparseRelation;
                    SparseRelation *rel_{};
                    size_t col_{0};
                    size_t col_start_{0};
                    size_t col_end_{0};

                    ColRange(SparseRelation *rel, size_t col) : rel_(rel), col_(col),
                        col_start_(rel->cols_first[col]),  col_end_(rel->cols_last[col]) {}

                public:

                    ColIterator begin() { return ColIterator(rel_, rel_->pool[col_start_].down); }
                    ColIterator end() { return ColIterator(rel_, col_end_); }      

                    ColIterator erase(ColIterator it) {
                        size_t current_cell = it.current_cell_;
                        if(current_cell == col_end_) { return it; }
                        size_t next = rel_->pool[current_cell].down;
                        rel_->unlink(current_cell);
                        return ColIterator(rel_, next);
                    }

                    size_t col(void) const noexcept { return col_; }
            };

            RowRange row(size_t row) { return RowRange(this, row); }
            ColRange col(size_t col) { return ColRange(this, col); }

            bool get(size_t row_idx, size_t col_idx) const {
                for(size_t current_col : const_cast<SparseRelation*>(this)->row(row_idx)) {
                    if(current_col == col_idx) { return true; }
                }
                return false;
            }

            bool erase(size_t row_idx, size_t col_idx);

            void set(size_t row_idx, size_t col_idx);

            size_t degree(void) const noexcept { return degree_; }

            void print_sr(void);

            void print_borders(void);

        private:

            void unlink(size_t cell_id) {
                auto& element = pool[cell_id];
                pool[element.left].right = element.right;
                pool[element.right].left = element.left;
                pool[element.up].down = element.down;
                pool[element.down].up = element.up;
            }

    };


    class FlattenNfa{
        private:

            struct StateRange {
                const size_t *begin_;
                const size_t *end_;
                inline const size_t *begin() const noexcept { return begin_; }
                inline const size_t *end() const noexcept { return end_; }
                inline bool empty() const noexcept { return begin_ == end_; }
            };

            std::vector<std::vector<size_t>> targets_{};
            std::vector<std::vector<size_t>> offsets_{};
            std::vector<std::vector<size_t>> rev_targets_{};
            std::vector<std::vector<size_t>> rev_offsets_{};
            std::vector<std::vector<size_t>> out_symbols_{};
            std::vector<hash::NfaHasher> hashers_{};
            std::vector<std::vector<hash::StateTransitionsHasher>> state_hashers_{};
            utils::SparseSet<size_t> valid_symbols_{};
            size_t num_of_states_{};
            size_t num_of_symbols_{};

        public:
            
            FlattenNfa(void) = default;
            FlattenNfa(const Nfa& aut, size_t num_of_symbols);
            const utils::SparseSet<size_t>& get_valid_symbols(void) const { return valid_symbols_; }

            bool has_successor(size_t state, size_t symbol) const noexcept { return (offsets_[symbol][state + 1] - offsets_[symbol][state]); }

            size_t get_bucket_count(void) {
                size_t border = (num_of_symbols_ + 3) / 4;
                size_t buckets = 1;
                while(buckets < border) { buckets <<= 1; }
                return buckets;
            }

            inline uint32_t get_bucket_idx(uint64_t hash1, uint64_t hash2, size_t num_of_buckets) noexcept {
                return uint32_t(hash::splitmix64(hash1 + hash::ADDITIVE_CONST * hash2) % num_of_buckets);
            }

            size_t num_of_states() const { return num_of_states_; }

            inline StateRange successors(size_t state, size_t symbol = 0) const noexcept {
                const size_t first = offsets_[symbol][state];
                const size_t last = offsets_[symbol][state + 1];
                return StateRange{targets_[symbol].data() + first, targets_[symbol].data() + last};
            }

            inline StateRange predecessors(size_t state, size_t symbol = 0) const noexcept {
                const size_t first = rev_offsets_[symbol][state];
                const size_t last = rev_offsets_[symbol][state + 1];
                return StateRange{rev_targets_[symbol].data() + first, rev_targets_[symbol].data() + last};
            }

            const std::vector<size_t>& out_symbols(size_t state) const noexcept { return out_symbols_[state]; }

            void print_nfa(void) {
                std::cout << num_of_symbols_ << " ";
                std::cout << "AUTOMATON: " << std::endl;
                for(size_t state = 0; state < num_of_states_; ++state) {
                    std::cout << state << ": " << std::endl;
                    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
                        std::cout << "   " << symbol << " -> ";
                        for(auto target : successors(state, symbol)) {
                        std::cout << target << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
                for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
                    std::cout << "HASHES FOR THE SYMBOL " << symbol << ": " 
                        << hashers_[symbol].get_hashes().hash1 << " " << hashers_[symbol].get_hashes().hash2 << std::endl;
                }
            }

            void print_rev_nfa(void) {
                std::cout << num_of_symbols_ << " ";
                std::cout << "AUTOMATON: " << std::endl;
                for(size_t state = 0; state < num_of_states_; ++state) {
                    std::cout << state << ": " << std::endl;
                    for(size_t symbol = 0; symbol < num_of_symbols_; ++symbol) {
                        std::cout << "   " << symbol << " -> ";
                        for(auto target : predecessors(state, symbol)) {
                        std::cout << target << " ";
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }

            bool are_interchangeable(size_t symbol1, size_t symbol2);
            bool is_covered_by(size_t symbol1, size_t symbol2);
            bool is_reversed_correctly(void);
    };

    
    class MergedTransitions {
        private:

            struct Block {
                size_t first_{};
                size_t last_{};
            };

            struct StateRange {
                const size_t *begin_;
                const size_t *end_;
                inline const size_t *begin() const noexcept { return begin_; }
                inline const size_t *end() const noexcept { return end_; }
                inline bool empty() const noexcept { return begin_ == end_; }
            };

            size_t symbol_{0};
            std::vector<Block> blocks_{};
            std::vector<Block> nodes_{};
            std::vector<size_t> data_{};
            size_t num_of_states_{};

        public:

            MergedTransitions(void) = default;

            MergedTransitions(FlattenNfa& aut, utils::Partition& part, size_t symbol);

            void split_block(FlattenNfa& aut, utils::Partition& part, size_t old_block_idx);

            StateRange get_block(size_t idx) const noexcept {
                const auto& block = blocks_[idx];
                return { data_.data() + block.first_, data_.data() + block.last_};
            }

            StateRange get_node(size_t idx) const noexcept {
                const auto& node = nodes_[idx];
                return { data_.data() + node.first_, data_.data() + node.last_};
            }

            void print_mt(void) {
                std::cout << "SYMBOL: " << symbol_ << std::endl;
                std::cout << "DATA:" << std::endl;
                for(size_t i = 0; i < data_.size(); ++i) {
                    std::cout << data_[i] << " ";
                }
                std::cout << std::endl;
                for(size_t i = 0; i < blocks_.size(); ++i) {
                    std::cout << "BLOCK " << i << ": ";
                    for(size_t idx = blocks_[i].first_; idx < blocks_[i].last_; ++idx) {
                        std::cout << data_[idx] << " ";
                    }
                    std::cout << "(" << blocks_[i].first_ << " " << blocks_[i].last_ << ")" << std::endl;
                }
                for(size_t i = 0; i < nodes_.size(); ++i) {
                    std::cout << "NODE " << i << ": ";
                    for(size_t idx = nodes_[i].first_; idx < nodes_[i].last_; ++idx) {
                        std::cout << data_[idx] << " ";
                    }
                    std::cout << "(" << nodes_[i].first_ << " " << nodes_[i].last_ << ")" << std::endl;
                }
            }
    };

}



namespace mata::nfa::simulation {

class Simulation {
    
    using PRP = utils::PartitionRelationPair;

    using MergedTransitions = detail::MergedTransitions;
    using SparseRelation = detail::SparseRelation;
    using SmartCounters = detail::SmartCounters;
    using FlattenNfa = detail::FlattenNfa;

    using Partition = utils::Partition;
    using Relation = utils::Relation;
    using PartitionRelationPair = utils::PartitionRelationPair;

    using State = utils::State;
    
    private:
        Nfa init_nfa_;
        FlattenNfa nfa_{};
        
        std::vector<MergedTransitions> mt_{};
        
        Partition part_{};
        Relation rel_{1, 1};
        SparseRelation sparse_rel_{};
        
        Partition init_part_{};
        Relation init_rel_{1, 1};
        SparseRelation init_sparse_rel_{};
        
        SmartCounters smart_counters_{};
        
        std::vector<State> repr_{};

        std::vector<std::vector<size_t>> not_rel_1{};
        std::vector<std::vector<size_t>> not_rel_2{};
        
        std::vector<std::vector<size_t>> in{};
        std::vector<std::vector<size_t>> out{};
        
        std::vector<std::vector<std::vector<size_t>>> pred_repr{};

        std::vector<size_t> pred_tmp{};
        
        GenerationSparseSet<size_t> refiner_1{};
        GenerationSparseSet<size_t> refiner_2{};
        
        utils::SparseSet<size_t> used_symbols{};
        
        std::vector<size_t> split_counts_{};
        
        GenerationSparseSet<size_t> touched_1{};
        GenerationSparseSet<size_t> touched_2{};
        
        GenerationSparseSet<size_t> noticed{};

        size_t num_of_symbols_ = 0;
        size_t num_of_states_ = 0;

        size_t inim_num_of_symbols{0};

        void state_structures_reset(const Nfa &nfa);
        void split_by_final_states(void);
        void sim_init(void);
        void split_update_data(size_t idx_c, size_t idx_d);
        bool split_first(void);
        size_t split_first_correction(size_t former_num_of_blocks);
        bool split_second(void);
        void refine(void);
        void split_blocks(void);
        bool at_least_one_split(size_t idx_b, size_t symbol);
        void sim_update_data(void);

        bool is_splitter_first(void);
        bool is_splitter_second(void);

        PartitionRelationPair relation_from_desymbolized_nfa(void);
        
        void refine_notice(size_t idx_a, size_t symbol);
        void refine_remove(size_t idx_a, size_t symbol);

        size_t get_num_of_symbols(void) const { return init_nfa_.delta.get_max_symbol() + 1; };
        
        void split_update_data_prepare(size_t idx_c, size_t idx_d);
        
        size_t skipped_splits{0};
        size_t performed_splits{0};
        size_t skipped_while_touched{0};

    public:
    
        Simulation(const Nfa& nfa) 
            : init_nfa_(nfa) {

            num_of_symbols_ = get_num_of_symbols();
            init_part_ = Partition(init_nfa_.num_of_states());
            init_rel_ = Relation(init_nfa_.num_of_states(), 1);
            init_rel_.set(0, 0, true);
            init_sparse_rel_ = SparseRelation();
            init_sparse_rel_.set(0, 0);

        }
        
        Simulation(const Nfa& nfa, Partition part, Relation rel) 
            : init_nfa_(nfa), init_part_(part), init_rel_(rel) {

            num_of_symbols_ = get_num_of_symbols();
            assert(nfa.num_of_states() == init_part_.num_of_states());
            assert(nfa.num_of_states() == init_rel_.capacity());
            assert(part_.num_of_blocks() == rel_.size());
            assert(rel.is_reflexive());
            assert(rel.is_antisymetric());
            assert(rel.is_transitive());

        }
            
        PartitionRelationPair compute_simulation();
        PartitionRelationPair compute_in_desymbolized_simulation(void);
             
}; // Simulation

} 

#endif // MATA_NFA_SIMULATION_HH
