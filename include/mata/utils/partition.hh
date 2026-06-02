/** @file partition.hh
 *  @brief Definition of a partition
 *
 *  In this context, we consider a carrier set S which contains all
 *  natural numbers from 0 to |S|-1 and nothing else. These numbers are called
 *  states.
 *  Then, partition over S is a set of blocks such that:
 *  - each block contains only states
 *  - each state is represented in exactly one block
 *     - blocks are disjoint
 *     - there is no state which is not represented in any block
 *  - no block is empty
 *
 *  This file provides implementation of a partition P which allows us to:
 *  - find the block which contains a given state in O(1)
 *  - find a representative state of the given block in O(1)
 *  - test whether two states share the same block in O(1)
 *  - test whether all states in a vector A share the same block in O(|A|)
 *  - iterate through the block B in O(|B|)
 *  - iterate through the node N in O(|N|)
 *  - split the whole partition such that each block
 *    is split in two pieces or remains unchanged in O(|S|)
 *  - remember all ancestors of current blocks and access them if necessary
 *    so we can manipulate multiple generations of a partition (before and
 *    after it has been split)
 *
 *  @author Tomáš Kocourek
 */

#ifndef LIBMATA_PARTITION_HH_
#define LIBMATA_PARTITION_HH_

#include <iostream>
#include <optional>
#include <vector>
#include <cassert>
#include "mata/utils/generation-sparse-set.hh"
#include "mata/utils/sparse-multiset.hh"

namespace mata::utils {

/************************************************************************
*
*
*
*                              PARTITION
*                   
*
*
*************************************************************************/

using State = unsigned long;
using StateBlock = std::vector<State>;
using StateBlocks = std::vector<StateBlock>;

/*
*   @struct SplitPair
*   @brief contains information about block which has been split
*
*   The structure SplitPair is created as soon as a block of the partition
*   is split. It contains:
*   - index of the new block which keeps the identity of the split block
*     (first_block_idx)
*   - index of the new block which has been created (second_block_idx)
*   - index of the node which had represented the former block before
*     (node_idx)
*   Using first_block_idx and second_block_idx, we can manipulate with the
*   current generation of the partition. Using node_idx, we are able to work
*   with the older generation of the partition.
*/
struct SplitPair {

    size_t first_block_idx;
    size_t second_block_idx;
    size_t node_idx;

    /**
    *   Initialization of the SplitPair   
    *
    *   @brief constructs a split pair
    *   @param[in] first first part of the split block
    *   @param[in] second second part of the split block
    *   @param[in] node node corresponding to the split block
    */
    SplitPair(size_t first, size_t second, size_t node) : 
        first_block_idx(first), second_block_idx(second), node_idx(node) {}
};

/**
 * @class Partition
 * @brief Partition of a set of states
 *
 * This data structure provides a partition of a set of states S. In this
 * context, the term 'state' refers to any natural number from the
 * interval <0, |S|-1>.
 *
 * This representation defines:
 * - states - element from a consecutive interval of natural numbers
 * - blocks - objects which represent the current generation of the partition.
 *            Each block refers to several states which belong to the block.
 *            The block could be possibly split.
 * - nodes  - objects which represent blocks either from the current generation
 *            of the partition or the block from the previous generations of
 *            the partition (block which had been split). Once a node is
 *            created, it is never changed. When a block is split, two new nodes
 *            are created.
 * - block items
 *          - objects which serve as intermediate data structure between states
 *            and blocks. Each block item contains indices of both corresponding
 *            state and block. Block items are sorted in such way that one could
 *            iterate through each block B or each node N in the O(B) or O(N)
 *            time, respectively.  
 * 
 * Detailed explanation:
 *
 * STATES:
 * States are represented by indices of the vector 'states_' with the
 * unchangeable size |S|. Each element of the vector represents the state
 * by its order in the 'states_' vector. The vector itself contains indices
 * of corresponding block items. Index of the block item corresponding
 * to the state 's' can be get in constant time using 'states_[s]'. 
 *
 * BLOCK ITEMS:
 * Block items are stored in the 'block_items_' vector of the unchangeable size
 * |S|. The block item can by accessed by its index using 
 * block_items[block_item_idx]. Each block item contains its index (which is
 * always the same as its order in the 'block_items_' vector but
 * it is stored directly in the object to simplify manipulations with
 * the partition). It also contains the index of its corresponding state
 * which means that states and block items are bijectively mapped. 
 * In addition, each BlockItem includes an index of the corresponding partition
 * block.
 * The ordering of 'block_items_' vector satisfies the condition that the states
 * of the same block (or node) should always form a contiguous subvector so one
 * could iterate through states in each block (or node) efficiently using
 * 'block_items' vector.
 *
 * BLOCKS:
 * The blocks themselves are represented by the vector 'blocks_' with the
 * size |P|, where P is a partition of states. Each block can be accessed by
 * its index 0 <= i < |P|. The block can by accessed by its index using 
 * blocks_[block_idx]. The block contains an index of itself and an index of its 
 * corresponding node. The total number of blocks can be changed as soon as
 * one block is split. However, the maximal number of blocks is equal to |S|
 * (the case when each block contains only one state). When a block 'B' is
 * split in two pieces 'B1' and 'B2', we create a brand new block 'B2' 
 * and modify the former block 'B' such that it will correspond to its
 * subblock 'B1'. The former block 'B' is thus not represented
 * in the 'blocks_' vector anymore since 'B1' takes over its identity.
 *
 * NODES:
 * Each node represents a current block or a block which has been split before.
 * The node can by accessed by its index using nodes_[node_idx]. If the given
 * node represents an existing block, such block contains an index of that node.
 * In context of nodes which represent former blocks, no block contains their
 * indices. The total number of nodes can be changed as soon as
 * one block is split. In such situation, two new nodes (which represent both
 * new blocks) are created and the former node remains unchanged. Therefore, the
 * maximal number of nodes is equal to 2 * |S| - 1 since once a node is created,
 * it is never changed. Each node contains its own index and two indices of
 * block items (namely 'first' and 'last') which could be used to access
 * block items corresponding to the first and last block items which form a   
 * contiguous subvector of 'block_items_' vector included in such node.
 * When a block is split, the corresponding block items are swapped 
 * in situ such that the indices 'first' and 'last' still surround the
 * corresponding node and both new nodes also point to the contiguous subvector
 * of its block items.
 * Each node N also contains at most three aditional indices:
 * - Its parent X, which is a direct ancestor of N (a node whose
 *   splitting caused existence of N). If there is no ancestor of N,
 *   this index is empty.
 * - Its left child B1, which is a direct descendant of N (a node which
 *   has been created when N was split and the corresponding new block
 *   inherited the identity of the former block)
 * - Its right child B2, which is a direct descendant of N (a node which
 *   has been created when N was split and the corresponding new block
 *   did not inherit the identity of the former block)
 * These indices form a binary tree which allows us to easily access former
 * blocks. If there is a child of a node N, there are always exactly two children
 * (both left and right). Othwerise, N has no child.
 *
 * STATE INFO:
 * Additional data structure which holds redundant information which can be then
 * accessed directly. Primarily, it directly connects state indices with
 * corresponding block indices and with current node index. Otherwise, it will
 * be neccessary to perform several dereferences and vector acessing to connect
 * a state with its block or node.
 * 
 * EXAMPLE:
 * In the example below, we represent a partition {{0, 1, 2, 3, 4}}
 * which had been split to the partition {{0, 2}, {1, 3, 4}}.
 * Thus, we have two blocks: 0 ({0, 2}) and 1 ({1, 3, 4}) which form the current
 * generation of the partition.
 * We also have three nodes: 0 ({0, 1, 2, 3, 4}), 1 ({0, 2}) and 2 ({1, 3, 4})
 * which represent both current generation of the partition and also block which
 * does not exist anymore since it had been split.
 * The block 0 corresponds to the node 1 and the block 1 corresponds to the 
 * node 2.
 * The node 1 contains indices 0 (first) and 1 (last) which means that
 * the blockItems 0 and 1 surround a contiguous subvector of elements in the 
 * node 1 (or in the block 0).
 * Likewise, the node 2 contains indices 2 (first) and 4 (last) which means that
 * the blockItems 2 and 4 surround a contiguous subvector of elements in the
 * node 2 (or in the block 1).
 * Moreover, we also represent the former block which does not exist anymore
 * by the node 0 which contains indices 0 (first) and 4 (last) which means that
 * the blockItems 0 and 4 surround a contiguous subvector of elements in the 
 * node 0. Thus, we know that there had been a block {0, 1, 2, 3, 4} before 
 * it has been split to obtain blocks {0, 2} and {1, 3, 4}.
 * In the picture below, indices of vectors are depicted outside of the vectors.
 * 
 *
 *             0       1       2       3       4
 *             ------- ------- ------- ------- -------
 *            |   0   |   2   |   1   |   4   |   3   |    states_
 *             ------- ------- ------- ------- -------
 *                ↑          ↑ ↑             ↑ ↑
 *                |          \ /             \ /
 *                |           X               X
 *                |          / \             / \
 *             0  ↓    1     ↓ ↓    2   3    ↓ ↓     4
 *             ------- ------- ------- ------- -------
 *            |   0   |   2   |   1   |   4   |   3   |    block_items_
 *       --→→→|-------|-------|-------|-------|--------←←←------------
 *       |    |   0   |   0   |   1   |   1   |   1   |               |
 *       |     ------- ------- ------- ------- -------                |
 *       |        |       |       |       |       |                   |
 *       |     0  ↓       ↓    1  ↓       ↓       ↓                   |
 *       |     ----------------------------------------               |
 *       |    |       1       |            2           |   blocks_    |
 *       |     ----------------------------------------               |
 *       |                |       |                                   |
 *       |     0       1  ↓    2  ↓                                   |
 *       |     ------- ------- -------                                |
 *       -----|   0   |   0   |   2   |   nodes_                      |
 *            |-------|-------|-------|                               |
 *            |   4   |   1   |   4   |                               |
 *             ------- ------- -------                                |
 *                |                                                   |
 *                ----------------------------------------------------
 * 
 *  Using this class, we can:
 *  - find the block which contains a given state in O(1)
 *  - find a representative state of the given block in O(1)
 *  - test whether two states share the same block in O(1)
 *  - test whether all states in a vector A share the same block in O(|A|)
 *  - iterate through the block B in O(|B|)
 *  - iterate through the node N in O(|N|)
 *  - split the whole partition such that each block
 *    is split in two pieces or remains unchanged in O(|S|)
 *  - remember all ancestors of current blocks and access them if necessary
 *    so we can manipulate multiple generations of a partition (before and
 *    after it has been split)
 *
 *  @invariant for each 0 <= i < |S|, block_items_[states_[i]].state() == i
 *  @invariant for each 0 <= i < |S|, states_[block_items_[i].state()] == i
 *  @invariant for each block index b, nodes_[blocks_[b].node_idx_].first().block().idx() == b
 *  @invariant for each subvector [b_i_0, b_i_1] of block_items_ vector such that
 *             block_items_[b_i_0].block_idx != block_items_[b_i_1].block_idx (= B),
 *             it holds that block_items[b_i_1].state_ is representative state of B
 *  @invariant each node is either a leaf node of it has exactly two descendants
 *  @invariant num_of_blocks() <= num_of_states
 *  @invariant num_of_blocks() <= num_of_nodes <= 2 * num_of_states() - 1
 *  @invariant if a state s is contained in a node N, then s is contained in all ancestors of N
 *  @invariant for each 0 <= i < |S|, states_info_[i].block_idx_ == block_items_[states_[i]].block_idx_
 *  @invariant for each 0 <= i < |S|, states_info_[i].node_idx_ == blocks_[block_items_[states_[i]].block_idx_].node_idx_
 *
 */
class Partition {

    public: 

        class StateInfo;
        class BlockItem;
        class Block;
        class Node;
     
    private:
    
        

        using StatesInfo = std::vector<StateInfo>;
        using States = std::vector<size_t>;
        using BlockItems = std::vector<BlockItem>;
        using Blocks = std::vector<Block>;
        using Nodes = std::vector<Node>;   
                
        /**< vector of information about each state */
        StatesInfo states_info_{};

        /**< vector of states refering to the block items */
        States states_{};
        
        /**< vector of block items refering to the states and blocks */
        BlockItems block_items_{};
        
        /**< vector of blocks refering to the nodes  */
        Blocks blocks_{};
        
        /**< vector of nodes refering to the first and last block item 
             of the node */
        Nodes nodes_{};
        
        /**< marked states used to perform splitting. After a split step,
             each block B could possibly be split in two pieces if there is
             both marked and non-marked states in B */
        GenerationSparseSet<State> marked_{};

        /**< if there is a marked state, its corresponding block is considered
             to be touched
             the multiset contains number of marked states in touched blocks
             
             for each touched block b, touched.occurences(b) correspond to that
             number
        */
        SparseMultiset<State> touched_{};

    public:

        /**
        * @class StateInfo
        * @brief Information about processed state
        *
        * Additional data structure which holds redundant information which can be then
        * accessed directly. Primarily, it directly connects state indices with
        * corresponding block indices and with current node index. Otherwise, it will
        * be neccessary to perform several dereferences and vector acessing to connect
        * a state with its block or node.
        *
        **/
        class StateInfo {
            private:
            
                /**< index of itself */
                size_t idx_;
                
                /**< index of the corresponding block_item */
                size_t block_item_idx_;
                
                /**< index of the corresponding block */
                size_t block_idx_;
                
                /**< index of the corresponding node */
                size_t node_idx_;

                /**< reference to the partition which works 
                     with this state info object */
                const Partition& partition_;
                
                // StateInfo class need to access private members
                // of Partition class 
                friend class Partition;

            public:

                /**
                * @brief constructs a state info element
                *
                * @param idx corresponding state
                * @param block_item_idx a block item index corresponding to the state
                * @param block_idx a block index corresponding to the state
                * @param node_idx a node index corresponding to the state
                * @param partition corresponding partition
                *
                * @par Complexity
                * O(1)
                */ 
                StateInfo(size_t idx, size_t block_item_idx, size_t block_idx,
                    size_t node_idx, const Partition& partition) :
                    
                    idx_(idx), block_item_idx_(block_item_idx), block_idx_(block_idx),
                    node_idx_(node_idx), partition_(partition) {}
                
                /**
                * @brief returns identity of the state info element
                *
                * @return state
                *
                * @par Complexity
                * O(1)
                */ 
                size_t idx(void) const { return idx_; }
                
                /**
                * @brief returns an index of a block item corresponding to the state
                *
                * @return block item index
                *
                * @par Complexity
                * O(1)
                */ 
                size_t block_item_idx(void) const { return block_item_idx_; }

                /**
                * @brief returns an index of a block which includes the state
                *
                * @return block index
                *
                * @par Complexity
                * O(1)
                */ 
                size_t block_idx(void) const { return block_idx_; }
                
                /**
                * @brief returns an index of a leaf node which includes the state
                *
                * @return node index
                *
                * @par Complexity
                * O(1)
                */
                size_t node_idx(void) const { return node_idx_; }


        };

        /**
        * @class BlockItem
        * @brief Intermediate between states and blocks
        *
        * Intermediate between states and blocks. Block items are ordered
        * such that the states belonging to the same block form a consecutive
        * subvector of the vector of block items.
        * States and block items are thus bijectively mapped.
        * Block items allow us to iterate through all states contained in one
        * block efficiently, whereas state vector allow us to access each state directly.
        *
        **/
        class BlockItem {
            private:
            
                /**< index of itself */
                size_t idx_;
                
                /**< corresponding state */
                State state_;
                
                /**< index of the corresponding block */
                size_t block_idx_;
                
                /**< reference to the partition which works 
                     with this block item */
                const Partition& partition_;
                
                // BlockItem class need to access private members
                // of Partition class 
                friend class Partition;

            public:

                /**
                * @brief constructs a block item element
                *
                * @param idx identity of the block item element
                * @param state a corresponding state
                * @param block_idx a block index corresponding to the block item
                * @param partition corresponding partition
                *
                * @par Complexity
                * O(1)
                */ 
                BlockItem(size_t idx, size_t state, size_t block_idx,
                    const Partition& partition) :
                    
                    idx_(idx), state_(state), block_idx_(block_idx),
                    partition_(partition) {}
                
                /**
                * @brief returns identity of the block item element
                *
                * @return identity of the block item
                *
                * @par Complexity
                * O(1)
                */ 
                size_t idx(void) const { return idx_; }

                /**
                * @brief returns the corresponding state
                *
                * @return state
                *
                * @par Complexity
                * O(1)
                */ 
                size_t state(void) const { return state_; }
                
                /**
                * @brief returns a constant reference to a block such that the
                * block item is contained in that block
                *
                * @return reference to a block
                *
                * @par Complexity
                * O(1)
                */ 
                const Block& block(void) const {
                    return partition_.blocks_[block_idx_];
                }

                /**
                * @brief returns a constant reference to a leaf node such that the
                * block item is contained in that node
                *
                * @return reference to a leaf node
                *
                * @par Complexity
                * O(1)
                */ 
                const Node& node(void) const { return block().node(); }
        };

        /**
        * @class Block
        * @brief Contains infromation about block from the current generation
        * of the partition.
        *
        * Each block represents a set of states (block of states of the current partition). 
        * It directly points to a single corresponding node.
        *
        **/        
        class Block {
            private:
            
                /**< index of itself */
                size_t idx_;
                
                /**< index of the corresponding node */
                size_t node_idx_;
                
                /**< reference to the partition which works 
                     with this block */
                const Partition& partition_;

                // Blocks need to access private members of Partition class                 
                friend class Partition;
                
            public:
            
                /**
                * @brief constructs a block element
                *
                * @param idx identity of the block element
                * @param node_idx a node index corresponding to the block item
                * @param partition corresponding partition
                *
                * @par Complexity
                * O(1)
                */ 
                Block(size_t idx, size_t node_idx, const Partition& partition) :
                    idx_(idx), node_idx_(node_idx),
                    partition_(partition) {}
                
                /**
                * @brief returns identity of the block element
                *
                * @return identity of the block
                *
                * @par Complexity
                * O(1)
                */ 
                size_t idx(void) const { return idx_; }
                
                /**
                * @brief returns a constant reference to a leaf node such that the
                * block is contained in that node
                *
                * @return reference to a leaf node
                *
                * @par Complexity
                * O(1)
                */ 
                const Node& node(void) const {
                    return partition_.nodes_[node_idx_];
                }  

                /**
                * @brief returns a constant reference to a very first block item
                * which is contained in this block
                *
                * @return reference to a first block item element of the block
                *
                * @par Complexity
                * O(1)
                */
                const BlockItem& first(void) const { return node().first(); }

                /**
                * @brief returns a constant reference to a very last block item
                * which is contained in this block
                *
                * @return reference to a first block item element of the block
                *
                * @par Complexity
                * O(1)
                */
                const BlockItem& last(void) const { return node().last(); }        
                
                /**< iterators which allow us to iterate through all states in the block */
                using const_iterator = typename BlockItems::const_iterator;
                    

                const_iterator begin() const {
                    const_iterator it = partition_.block_items_.begin();
                    std::advance(it, static_cast<long>(node().first().idx()));
                    return it;
                }
                
                const_iterator end() const { 
                    const_iterator it = partition_.block_items_.begin();
                    std::advance(it, static_cast<long>(node().last().idx()) +1);
                    return it;
                }
                
                /**
                * @brief number of elements in this block
                *
                * @return size of the block
                *
                * @par Complexity
                * O(1)
                */
                [[nodiscard]] size_t size(void) const noexcept { 
                    return last().idx() - first().idx() + 1;
                }
                
        };

        /**
        * @class Node
        * @brief Contains infromation about block from the current generation
        * of the partition or from the previous generation of the partition.
        *
        * Each node represents a set of states (block of states of the current partition
        * or the former partition which had been split). 
        * It directly points to both first and last elements of the node stored in a vector
        * of block items.
        *
        **/         
        class Node {
            private:
            
                /**< index of itself */
                size_t idx_;
                
                /**< index of the first block item in the node */
                size_t first_;
                
                /**< index of the last block item in the node */
                size_t last_;

                /**< direct ancestor of the node */
                std::optional<size_t> parent_{std::nullopt};

                /**< direct descendants of the node */
                std::optional<size_t> left_child_{std::nullopt};
                std::optional<size_t> right_child_{std::nullopt};
                
                /**< reference to the partition which works 
                     with this block */                
                const Partition& partition_;
                
                // Nodes need to access private members of Partition class
                friend class Partition;

                // setting children node as soon as the corresponding block is split
                void set_children(size_t left_child, size_t right_child) {
                    left_child_ = left_child;
                    right_child_ = right_child;
                }
                
            public:
            
                /**
                * @brief constructs a node element
                *
                * @param idx identity of the block element
                * @param first very first block item in the node
                * @param last very last block item in the node
                * @param parent optional index of the node's parent
                * @param partition corresponding partition
                *
                * @par Complexity
                * O(1)
                */
                Node(size_t idx, size_t first, size_t last, std::optional<size_t> parent,
                    const Partition& partition) :
                    
                    idx_(idx), first_(first), last_(last), parent_(parent), partition_(partition)
                    {}
            
                /**
                * @brief returns identity of the node element
                *
                * @return identity of the node
                *
                * @par Complexity
                * O(1)
                */ 
                size_t idx(void) const { return idx_; }
                
                /**
                * @brief returns a constant reference to a very first block item
                * which is contained in this node
                *
                * @return reference to a first block item element of the node
                *
                * @par Complexity
                * O(1)
                */
                const BlockItem& first(void) const {
                    return partition_.block_items_[first_];
                }

                /**
                * @brief returns a constant reference to a very last block item
                * which is contained in this node
                *
                * @return reference to a first block item element of the node
                *
                * @par Complexity
                * O(1)
                */
                const BlockItem& last(void) const { 
                    return partition_.block_items_[last_];
                }
                
                /**< iterators which allow us to iterate through all states in the node */
                using const_iterator = 
                    typename BlockItems::const_iterator;
                    
                const_iterator begin() const {  
                    const_iterator it = partition_.block_items_.begin();
                    std::advance(it, static_cast<long>(first().idx()));
                    return it;
                }
                
                const_iterator end() const { 
                    const_iterator it = partition_.block_items_.begin();
                    std::advance(it, static_cast<long>(last().idx()) + 1);
                    return it;
                }
                
                /**
                * @brief number of elements in this node
                *
                * @return size of the node
                *
                * @par Complexity
                * O(1)
                */
                [[nodiscard]] size_t size(void) const { 
                    return last().idx() - first().idx() + 1;
                }

                /**
                * @brief returns a direct ancestor of this node in the node tree
                *
                * If the ancestor does not exist, nullopt will be used.
                *
                * @return parent of the node
                *
                * @par Complexity
                * O(1)
                */
                std::optional<size_t> parent(void) const { return parent_; }

                /**
                * @brief returns a first direct descendant of this node in the node tree
                *
                * If the descendant does not exist, nullopt will be used.
                *
                * @return child of the node
                *
                * @par Complexity
                * O(1)
                */
                std::optional<size_t> left_child(void) const { return left_child_; }

                /**
                * @brief returns a second direct descendant of this node in the node tree
                *
                * If the descendant does not exist, nullopt will be used.
                *
                * @return child of the node
                *
                * @par Complexity
                * O(1)
                */
                std::optional<size_t> right_child(void) const { return right_child_; }
                
                /**
                * @brief checks whether the given block is contained in this node
                *
                * @return true iff block_idx corresponds to the block contained in this node
                *
                * @par Complexity
                * O(1)
                */
                bool contains_block(size_t block_idx) const {
                    const Block& block = partition_.blocks_[block_idx];
                    return first_ <= block.first().idx() && last_ >= block.last().idx(); 
                }
        };

        /**
        * @class BlocksInNode
        * @brief Additional class which allows us to iterate through blocks contained
        * in one node by traversing node tree
        *
        * The node tree of the partition is traversed using depth first search (simulated by a stack).
        * As soon as a leaf node is reached, it is returned.
        *
        **/
        class BlocksInNode {
            private:
                const Partition& partition_;
                std::vector<size_t> stack_{};
                size_t current_block_ = 0;
                size_t node_idx_ = 0;
                bool is_done_ = false;

            public:

                class iterator {
                    public:
                        iterator(BlocksInNode& range, bool end = false) : range_(range), end_(end) {}

                        size_t operator*() const { return range_.current_block_; }

                        iterator& operator++() {
                            range_.advance();
                            end_ = range_.is_done_;
                            return *this;
                        }

                        friend bool operator!=(const iterator& it_1, const iterator& it_2) {
                            return it_1.end_ != it_2.end_;
                        }

                    private:
                        BlocksInNode& range_;
                        bool end_;
                };

                BlocksInNode(const Partition& p, size_t node_idx) : partition_(p), node_idx_(node_idx) {
                    stack_.push_back(node_idx);
                    advance();
                }

                /**
                * @brief resets the object so it can be iterated through again
                *
                * @par Complexity
                * O(1)
                */
                void reset(void) {
                    stack_.clear();
                    stack_.push_back(node_idx_);
                    is_done_ = false;
                    current_block_ = 0;
                }

                iterator begin() { return iterator(*this, is_done_); }
                iterator end() { return iterator(*this, true); }

            private:

                /**
                * @brief performs DFS traversal steps until a next leaf node is reached
                * or the tree is completely traversed
                *
                * @par Complexity
                * O(|N|), where N is the set of nodes (~ O(|S|) with S as a number of states)
                */
                void advance() {
                    while(!stack_.empty()) {
                        size_t node_idx = stack_.back();
                        stack_.pop_back();

                        const Node& node = partition_.nodes_[node_idx];

                        if(!node.left_child()) {
                            current_block_ = node.first().block().idx();
                            return;
                        }

                        stack_.push_back(*node.right_child());
                        stack_.push_back(*node.left_child());
                    }
                    is_done_ = true;
                }
        };
        
        // constructors
        Partition() = default;
        Partition(size_t num_of_states, 
                  const StateBlocks& partition = StateBlocks());
        Partition(const Partition& other);
        Partition(Partition&& other);
        Partition& operator=(Partition&& other);

        /**
        * @brief returns number of states
        *
        * @return number of states
        *
        * @par Complexity
        * O(1)
        */
        size_t num_of_states(void) const { return states_info_.size(); }

        /**
        * @brief returns number of block items
        *
        * @return number of block items
        *
        * @par Complexity
        * O(1)
        */
        size_t num_of_block_items(void) const { return block_items_.size(); }
        
        /**
        * @brief returns number of blocks
        *
        * @return number of blocks
        *
        * @par Complexity
        * O(1)
        */
        size_t num_of_blocks(void) const { return blocks_.size(); }
        
        /**
        * @brief returns number of nodes
        *
        * @return number of nodes
        *
        * @par Complexity
        * O(1)
        */
        size_t num_of_nodes(void) const { return nodes_.size(); }
        
        /**
        * @brief returns number of currently marked states
        *
        * @return number of marked states
        *
        * @par Complexity
        * O(1)
        */
        size_t num_of_marked_states(void) const { return marked_.size(); }
        
        /**
        * @brief splits a partition according to a given sparse set of markings
        *
        * Each block B will be either:
        * - remain unchanged iff all states of B are either marked or unmarked
        * - split in two pieces B1 (all marked states of B) and B2 (all unmarked states of B)
        *
        * @param marked sparse set of marked states
        *
        * @return vector of tuples (first new block index, second new block index)
        *
        * @par Complexity
        * O(|A|), where A is a set of all states
        */
        std::vector<SplitPair> split_blocks(GenerationSparseSet<State>& marked);

        /**
        * @brief splits a partition according to a given sparse set of markings
        * and a block index B
        *
        * The block B will be either:
        * - remain unchanged iff all states of B are either marked or unmarked
        * - split in two pieces B1 (all marked states of B) and B2 (all unmarked states of B)
        *
        * @param marked sparse set of marked states
        * @param block_idx index of a block to be split
        *
        * @return vector of tuples (first new block index, second new block index)
        *
        * @par Complexity
        * O(|A|), where A is a set of all states
        */
        std::vector<SplitPair> split_blocks(
            const GenerationSparseSet<size_t>& marked, size_t block_idx);

        /**
        * @brief splits a partition according to marked states and a given block index B.
        * The states are marked in the inner representation of the partition.
        *
        * The block B will be either:
        * - remain unchanged iff all states of B are either marked or unmarked
        * - split in two pieces B1 (all marked states of B) and B2 (all unmarked states of B)
        *
        * @return vector of tuples (first new block index, second new block index)
        *
        * @par Complexity
        * O(|A|), where A is a set of all states
        */
        std::vector<SplitPair> split_blocks_inner(size_t block_idx);

        /**
        * @brief maps a given state S to the corresponding block B such that
        * S is contained in B 
        *
        * @param state state to be mapped to a block
        *
        * @return index of a block
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_block_idx_of_state(State state) const {
            assert(state < num_of_states() && "Nonexisting state used");
            return states_info_[state].block_idx_;
        }

        /**
        * @brief maps a given state S to the corresponding block item BI such that
        * S is contained in BI
        *
        * @param state state to be mapped to a block item
        *
        * @return index of a block item
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_block_item_idx_of_state(State state) const {
            assert(state < num_of_states() && "Nonexisting state used");
            return states_info_[state].block_item_idx_;
        }

        /**
        * @brief maps a given state S to the corresponding node N such that
        * S is contained in N and N corresponds to an existing block B
        *
        * @param state state to be mapped to a node
        *
        * @return index of a node
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_node_idx_of_state(State state) const {
            assert(state < num_of_states() && "Nonexisting state used");
            return states_info_[state].node_idx_;
        }

        /**
        * @brief maps a given state S to the corresponding state E such that
        * E is a representative state of the block B such that S is contained in B
        *
        * @param state state to be mapped to a representative state
        *
        * @return index of a represetative state
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] State get_block_repr_of_state(State state) const {
            assert(state < num_of_states() && "Nonexisting state used");
            return nodes_[states_info_[state].node_idx_].first().state();
        }

        /**
        * @brief maps a given block index B to the corresponding state E such that
        * E is a representative state of the block B
        *
        * @param block_idx block index to be mapped to a block
        *
        * @return index of a block
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] State get_block_repr_of_block_idx(size_t block_idx) const {
            assert(block_idx < num_of_blocks() && "Nonexisting block index used");
            return blocks_[block_idx].first().state();
        }

        /**
        * @brief maps a given block index B to the corresponding node N such that
        * N is a node of the current generation
        *
        * @param block_idx block index to be mapped to a node
        *
        * @return index of a node
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_node_idx_of_block_idx(size_t block_idx) const {
            assert(block_idx < num_of_blocks() && "Nonexisting block index used");
            return blocks_[block_idx].node_idx_;
        }

        /**
        * @brief returns a size of the given block
        *
        * @param block_idx block index
        *
        * @return block size
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_size_of_block_idx(size_t block_idx) const {
            assert(block_idx < num_of_blocks() && "Nonexisting block index used");
            return blocks_[block_idx].size();
        }

        /**
        * @brief maps a given node index N to the corresponding state E such that
        * E is a representative state of N
        *
        * @param node_idx node index to be mapped to a block representative
        *
        * @return index of a representative state
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] State get_block_repr_of_node_idx(size_t node_idx) const {
            assert(node_idx < num_of_nodes() && "Nonexisting node index used");
            return nodes_[node_idx].first().state();
        }

        /**
        * @brief maps a given node index N to the corresponding block index B such that
        * B and N share their representative states
        *
        * @param node_idx node index to be mapped to a block index
        *
        * @return index of a block
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_block_idx_of_node_idx(size_t node_idx) const {
            assert(node_idx < num_of_nodes() && "Nonexisting node index used");
            return nodes_[node_idx].first().block().idx();
        }

        /**
        * @brief returns a size of the given node
        *
        * @param node_idx node index
        *
        * @return node size
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t get_size_of_node_idx(size_t node_idx) const {
            assert(node_idx < num_of_nodes() && "Nonexisting node index used");
            return nodes_[node_idx].size();
        }

        /**
        * @brief decides whether the given block B is fully contained in the given node N
        *
        * @param node_idx node index N
        * @param block_idx block index B
        *
        * @return true iff B is contained in N
        *
        * @par Complexity
        * O(1)
        */
        bool node_contains_block(size_t node_idx, size_t block_idx) const {
            assert(block_idx < num_of_blocks() && "Nonexisting block used");
            assert(node_idx < num_of_nodes() && "Nonexisting node used");
            size_t repr = blocks_[block_idx].first().idx_;
            return repr >= nodes_[node_idx].first_ && repr <= nodes_[node_idx].last_;
        }

        /**
        * @brief decides whether the given block B is fully contained in the given node N
        *
        * @param node_idx node index
        *
        * @return true iff B is contained in N
        *
        * @par Complexity
        * O(1)
        */
        BlocksInNode get_blocks_in_node(size_t node_idx) const {
            assert(node_idx < num_of_nodes() && "Nonexisting node used");
            return BlocksInNode(*this, node_idx);
        }

        /**
        * @brief decides whether the given node N does not have any descendant
        *
        * @param node_idx node index
        *
        * @return true iff N does not have any descendants
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool is_leaf_node(size_t node_idx) const {
            assert(node_idx < num_of_nodes() && "Nonexisting node used");
            return !nodes_[node_idx].left_child();
        }

        /**
        * @brief decides whether two states are contained in the same block
        *
        * @param first first state to be checked
        * @param second second state to be checked
        *
        * @return true iff first are second are contained in the same block
        *
        * @par Complexity
        * O(1)
        */
        bool in_same_block(State first, State second) const;

        /**
        * @brief decides whether a vector of states contains only states which
        * are contained in the same block
        *
        * @param states vector of states
        *
        * @return true iff all states from the vector states are contained in the same set
        *
        * @par Complexity
        * O(|states|)
        */
        bool in_same_block(const std::vector<State>& states) const;
        
        /**
        * @brief creates a vector of states such that they they are all in the same block as A
        *
        * @param state given state A
        *
        * @return vector of all states such they share the same block as A
        *
        * @par Complexity
        * O(|A|), where A is set of all states
        */
        std::vector<State> states_in_same_block(State state) const;

        // accessing state info, block items, blocks, nodes through indices
        const StateInfo& get_state_info(size_t state_info_idx) const;
        const BlockItem& get_block_item(size_t block_item_idx) const;
        const Block& get_block(size_t block_idx) const;
        const Node& get_node(size_t node_idx) const;
    
        /**
        * @brief converts a partition to a vector of vectors of states
        *
        * @return vector of vector of states
        *
        * @par Complexity
        * O(|A|), where A is set of all states
        */
        StateBlocks partition(void) const;
        
        /**
        * @brief marks a given state in an inner representation
        *
        * @par Complexity
        * O(1)
        */
        void mark(State state);

        /**
        * @brief deletes all markings
        *
        * @par Complexity
        * O(1)
        */
        void clear_markings(void);

        /**
        * @brief decides whether a given state is marked
        *
        * @param state given state
        *
        * @return true iff the state is marked
        *
        * @par Complexity
        * O(1)
        */
        bool is_marked(State state) const;

        /**
        * @brief decides whether a given block contains at least one marked state
        *
        * @param block_idx block index to be checked
        *
        * @return true iff the block is touched
        *
        * @par Complexity
        * O(1)
        */
        bool is_touched(size_t block_idx) const;

        /**
        * @brief decides whether a given block contains only marked states
        *
        * @param block_idx block index to be checked
        *
        * @return true iff the block if fully marked (all states in the block are marked)
        *
        * @par Complexity
        * O(1)
        */
        bool is_fully_marked(size_t block_idx) const;
        
        /**
        * @brief additional function which prints marked blocks
        *
        * @par Complexity
        * O(|A|), where A is a set of all states
        */
        void print_marked_blocks(void);
        
        /**
        * @brief returns a sparse set of touched blocks
        *
        * @return sparse set of touched blocks
        *
        * @par Complexity
        * O(1)
        */
        const SparseMultiset<State> &touched_blocks(void) const { return touched_; }

        // operators
        Partition& operator=(const Partition& other);        
        friend std::ostream& operator<<(std::ostream& os, 
                                        const Partition& p);
        
        const BlockItem& operator[](State state) const;
        bool operator==(const Partition& other) const;
            
}; // Partition

}

#endif // LIBMATA_PARTITION_HH_
