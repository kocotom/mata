#ifndef LIBMATA_SPARSE_MULTISET_HH
#define LIBMATA_SPARSE_MULTISET_HH

#include <cassert>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <limits>
#include <cstdint>
#include <initializer_list>

/**
* @class SparseMultiset
* @brief Sparse multiset of unsigned integral values from a fixed-size domain
*
* A sparse multiset stores nonnegative integral values from the domain [0, domain_size()).
* The data structure allows us to insert a chosen number multiple times. It will
* then remember the number of occurrences. The data structure provides functions
* for constant-time insertion, membership test, number of occurrences test and erasure. 
* Clearing the whole sparse multiset behaves like a constant-time function
* in most cases, O(domain_size()) in extremely rare situations. Iteration through
* contained values is performed in O(size()).
*
* A sparse multiset consists of four vectors which combine information about
* stored values and number of their occurrences.
*
* 1. sparse_:
*    Vector containing domain_size() elements. Each element corresponds to a generation
*    index. If sparse_[i] == generation_ (current generation), then a value `i` is in fact
*    contained in the sparse multiset (at least one occurrence). 
*    Otherwise, the element `i` has 0 occurrences in the sparse multiset.
*    Initially, each generation index is set to 0. The initial generation equals 1, therefore
*    the multiset is initially empty.
*    The current generation index is changed as soon as the clear() function is called. Thus,
*    there is no need to manually erase each sparse_[i] value when the clearing of the whole multiset
*    is performed since the former value corresponds to an old, outdated generation.
*    As soon as the maximal possible generation is reached, each element in sparse_[i] is manually
*    set to 0 and the generation index is manually set to 1.
*    This approach ensures that the membership test works in O(1) and the clearing behaves like
*    constant-time function in most cases.
*    Invariant: sparse_[i] == generation_ iff `i` is contained in the sparse multiset.
* 2. dense_:
*    Vector containing size() elements. Each element corresponds to a value contained in the sparse
*    multiset. There is no guaranteed order of the elements in the dense_ vector. However,
*    it is possible to iterate through all elements contained in the sparse multiset in O(size()) rather
*    than O(domain_size()).
* 3. dense_idx_:
*    Vector containing domain_size() elements. It maps elements from the domain to their positions
*    in the dense_ vector to be able to efficiently erase them in O(1).
*    Invariant: If `i` is contained in the sparse multiset, then dense_[dense_idx_[i]] == i.
* 4. counts_:
*    Vector containing size() elements. The value stored in counts_[i] corresponds to a number
*    of occurrences of the element stored in dense_[i].
*
* The main purposes of this data structure (in comparison with other implementations of sparse set) are:
* 1. The ability to clear the set in constant time (typically)
*    The clear() operation is crucial and it will be called very often. It will be highly inefficient
*    to always iterate through the whole sparse_ vector to invalidate it. In contrast, the incrementation
*    of a current generation is extremely fast.
* 2. The ability to perform a membership test as fast as possible
*    Classical implementations of sparse sets allow this in O(1). However, it typically
*    includes several constant operations and vector accesses to evaluate all membership conditions.
*    Since the membership test will be performed very often, the generation approach is more efficient.
*
* @warning For the sake of efficiency, this implementation does not control domain boundaries in release
* build. The caller is always considered to be responsible for passing allowed values.
*
* @tparam T Unsigned integral type of values stored in the sparse multiset
*
* @invariant Each value stored in the sparse multiset is smaller than domain_size_.
*/
template <typename T>
class SparseMultiset {
    using generation_t = std::uint32_t;
    
    static constexpr generation_t MAX_GENERATION = std::numeric_limits<generation_t>::max();

    static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>,
                  "SparseMultiset requires an unsigned integral type");

    private:
    
        std::vector<generation_t> sparse_{};
        std::vector<T> dense_{};
        std::vector<size_t> dense_idx_{};
        std::vector<size_t> counts_{};

        size_t domain_size_{1};
        generation_t generation_{1};

    public:

        /**
        * @brief Construct an empty sparse multiset according to a given domain size
        *
        * This function creates an empty sparse multiset. The domain size has to be
        * at least 1, otherwise it is set to 1. Sparse vectors are resized to the domain size,
        * whereas the dense vectors reserve enough capacity.
        *
        * @param domain_size initial domain size
        * @par Complexity
        * O(domain_size)
        */        
        explicit SparseMultiset(size_t domain_size = 1) {
            if(!domain_size) { domain_size = 1; }
            sparse_.resize(domain_size, 0);
            dense_idx_.resize(domain_size, 0);
            dense_.reserve(domain_size);
            counts_.reserve(domain_size);
            domain_size_ = domain_size;
        }
        
        /**
        * @brief Construct a sparse multiset according to given values and domain size
        *
        * @param values initial values
        * @param domain_size initial domain size
        * @par Complexity
        * O(domain_size) + O(|values|)
        */ 
        SparseMultiset(std::initializer_list<T> values, size_t domain_size = 1) 
        : SparseMultiset(domain_size) {
            for(T val : values) { insert(val); }
        }

        SparseMultiset(const SparseMultiset&) = default;
        SparseMultiset& operator=(const SparseMultiset&) = default;
        SparseMultiset(SparseMultiset&&) noexcept = default;
        SparseMultiset& operator=(SparseMultiset&&) noexcept = default;

        /**
        * @brief Returns a domain size (upper bound of the allowed values)
        *
        * @return domain size
        *
        * @par Complexity
        * O(1)
        */ 
        [[nodiscard]] size_t domain_size(void) const noexcept { return domain_size_; }
        
        /**
        * @brief Returns the number of distinct values contained in the sparse multiset
        *
        * This corresponds to the size of the dense_ vector. The function returns
        * number of distinct elements of the domain such that they are stored in the multiset
        * at least once.
        *
        * @return number of elements stored in the multiset
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t size(void) const noexcept { return dense_.size(); }
        
        /**
        * @brief Inserts a given element to the sparse multiset
        *
        * Inserts a given element to the sparse multiset. If the element has already been
        * inserted to the multiset, the number of occurrences increases.
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain to be inserted
        *
        * @return number of current occurrences of the inserted value
        *
        * @par Complexity
        * O(1)
        */ 
        size_t insert(T val) {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            if(sparse_[val] == generation_) { return ++counts_[dense_idx_[val]]; }
            sparse_[val] = generation_;
            dense_idx_[val] = dense_.size();
            dense_.push_back(val);
            counts_.push_back(1);
            return 1;
        }

        /**
        * @brief Change a number of occurrences of a given value in the sparse multiset
        *
        * If a given element is not stored in the sparse multiset, it will be inserted.
        * Then the number of occurrences is updated. If the new count is set to 0,
        * the element will be erased.
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        * @param count new number of occurrences
        * @par Complexity
        * O(1)
        */ 
        void set_count(T val, size_t count) {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            if(!count) { 
                erase(val); 
                return;
            }
            if(!occurrences(val)) { insert(val); }
            counts_[dense_idx_[val]] = count;
        }
        
        /**
        * @brief Erases all elements from the sparse multiset
        *
        * This function erases all elements from the sparse multiset. It increases the current
        * generation index which invalidates all generation marks in the sparse_ vector. If the maximal
        * generation index is reached, the sparse_ vector will be cleared manually in linear time.
        *
        * @par Complexity
        * O(domain_size()) in an extremely rare situation when the maximal generation is reached,
        * O(1) otherwise.
        */ 
        void clear(void) {
            dense_.clear();
            counts_.clear();
            if (generation_ == MAX_GENERATION) [[unlikely]] {
                std::fill(sparse_.begin(), sparse_.end(), 0);
                generation_ = 1;
            } else {
                ++generation_;
            }
        }
        
        /**
        * @brief Widens the domain size
        *
        * The function widens the domain size according to a parameter `size`
        * if it is greater than the current domain size. Otherwise, it has no effect.
        *
        * @param size new size of the domain
        *
        * @par Complexity
        * O(size) if the domain is extended
        * O(1) otherwise
        */ 
        void extend_domain(size_t size) {
            if(size <= domain_size_) { return; }
            domain_size_ = size;
            sparse_.resize(size, 0);
            dense_idx_.resize(size, 0);
        }
        
        /**
        * @brief Returns a boolean flag which says whether the multiset is empty
        *
        * @return boolean flag corresponding to the emptiness of the multiset
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool empty(void) const noexcept { return dense_.empty(); }

        /**
        * @brief Returns a number of occurrences of the given element
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        *
        * @return number of occurrences of the given element
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t occurrences(T val) const noexcept {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            if(sparse_[val] != generation_) { return 0; }
            return counts_[dense_idx_[val]];
        }
        
        /**
        * @brief Returns a number of occurrences of the given element.
        * This operator is equivalent to contains(val).
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        *
        * @return number of occurrences of the given element
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t operator[](T val) const noexcept {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            return occurrences(val);
        }

        /**
        * @brief Returns a boolean flag which says whether the given value is
        * contained in the multiset
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        *
        * @return boolean flag which indicates presence of the element
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool contains(T val) const noexcept {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            return sparse_[val] == generation_;
        }

        /**
        * @brief Erases a given value from the multiset
        *
        * Erases all occurrences of a value from the multiset. If the given element is not
        * contained in the multiset, the function won't have any effect.
        *
        * @pre val < domain_size()
        *
        * @param val element to be erased
        *
        * @par Complexity
        * O(1)
        */
        void erase(T val) noexcept {
            assert(val < domain_size_ && "Value outside of the SparseMultiset domain.");
            if (sparse_[val] != generation_) { return; }
            sparse_[val] = 0;
            size_t idx = dense_idx_[val];
            size_t last = dense_.size() - 1;
            T last_val = dense_[last];
            dense_[idx] = last_val;
            counts_[idx] = counts_[last];
            dense_idx_[last_val] = idx; 
            dense_.pop_back();
            counts_.pop_back();
        }

        auto begin() noexcept { return dense_.begin(); }
        auto end() noexcept { return dense_.end(); }

        auto begin() const noexcept { return dense_.begin(); }
        auto end() const noexcept { return dense_.end(); }

        auto cbegin() const noexcept { return dense_.cbegin(); }
        auto cend() const noexcept { return dense_.cend(); }
};



#endif //LIBMATA_SPARSE_MULTISET_HH
