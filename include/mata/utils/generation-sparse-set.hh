#ifndef LIBMATA_GENERATION_SPARSE_SET_HH
#define LIBMATA_GENERATION_SPARSE_SET_HH

#include <cassert>
#include <algorithm>
#include <vector>
#include <type_traits>
#include <limits>
#include <cstdint>
#include <initializer_list>
#include <cstddef>

/**
* @class GenerationSparseSet
* @brief Generation sparse set (shortly denoted as GSS) of unsigned integral 
* values from a fixed-size domain
*
* A GSS stores nonnegative integral values from the domain [0, domain_size()).
* The data structure provides functions for constant-time insertion and membership test.
* Clearing the whole GSS takes constant time in most cases, O(domain_size())
* in extremely rare situations. Iteration through contained values is performed in O(size()).
* Erasing a single element takes O(size()) time.
*
* A GSS consists of two vectors which hold information about stored values.
*
* 1. sparse_:
*    Vector containing domain_size() elements. Each element corresponds to a generation
*    index. If sparse_[i] == generation_ (current generation), then a value `i` is in fact
*    contained in the GSS. Otherwise, the element `i` is not stored in the GSS.
*    Initially, each generation index is set to 0. The initial generation equals 1, therefore
*    the GSS is initially empty.
*    The current generation index is changed as soon as the clear() function is called. Thus,
*    there is no need to manually erase each sparse_[i] value when the clearing of the whole GSS
*    is performed since the former value corresponds to an old, outdated generation.
*    When the maximal possible generation is reached, each element in sparse_[i] is manually
*    set to 0 and the generation index is manually set to 1.
*    This approach ensures that the membership test works in O(1) and the clearing behaves like
*    constant-time function in most cases.
*    Invariant: sparse_[i] == generation_ iff `i` is contained in the GSS.
* 2. dense_:
*    Vector containing size() elements. Each element corresponds to a value contained in the GSS. 
*    There is no guaranteed order of the elements in the dense_ vector. However,
*    it is possible to iterate through all elements contained in the GSS in O(size()) rather
*    than O(domain_size()).
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
* 3. No need to erase individual elements
*    This data structure is primarily constructed for such purposes that there will be no need to erase
*    single elements from the GSS. In contrast, elements are expected to be cleared all at once. Thus,
*    there is no need to tightly bind indices from sparse_ and dense_ vector and the membership test and
*    insertion can be performed faster. However, an inefficient erase function is provided for completeness.
*
* @warning For the sake of efficiency, this implementation does not control domain boundaries in release
* build. The caller is always considered to be responsible for passing allowed values.
*
* @warning For the sake of efficiency, the erase() function should be used as little as possible
* since it requires linear-time w.r.t. number of stored elements. The GSS should be
* primarily used for purposes such that elements are cleared all at once. 
*
* @tparam T Unsigned integral type of values stored in the GSS
*
* @invariant Each value stored in the GSS is smaller than domain_size_.
*/
template <typename T>
class GenerationSparseSet {
    using generation_t = std::uint32_t;

    static constexpr generation_t MAX_GENERATION = std::numeric_limits<generation_t>::max();
    
    static_assert(std::is_integral_v<T> && std::is_unsigned_v<T>,
                  "GenerationSparseSet requires an unsigned integral type");

    private:
    
        std::vector<generation_t> sparse_{};
        std::vector<T> dense_{};

        size_t domain_size_{1};
        generation_t generation_{1};

    public:

        /**
        * @brief Construct an empty GSS according to a given domain size
        *
        * This function creates an empty GSS. The domain size has to be
        * at least 1, otherwise it is set to 1. Sparse vector is resized to the domain size,
        * whereas the dense vector reserves enough capacity.
        *
        * @param domain_size initial domain size
        * @par Complexity
        * O(domain_size)
        */    
        explicit GenerationSparseSet(size_t domain_size = 1) {
            if(!domain_size) { domain_size = 1; }
            sparse_.resize(domain_size, 0);
            dense_.reserve(domain_size);
            domain_size_ = domain_size;
        }
        
        /**
        * @brief Construct a GSS according to given values and domain size
        *
        * @param values initial values
        * @param domain_size initial domain size
        * @par Complexity
        * O(domain_size) + O(|values|)
        */ 
        GenerationSparseSet(std::initializer_list<T> values, size_t domain_size = 1) 
        : GenerationSparseSet(domain_size) {
            for(T val : values) { insert(val); }
        }

        GenerationSparseSet(const GenerationSparseSet&) = default;
        GenerationSparseSet& operator=(const GenerationSparseSet&) = default;
        GenerationSparseSet(GenerationSparseSet&&) noexcept = default;
        GenerationSparseSet& operator=(GenerationSparseSet&&) noexcept = default;

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
        * @brief Returns the number of values contained in the GSS
        *
        * This corresponds to the size of the dense_ vector. The function returns
        * number of elements of the domain such that they are stored in the GSS.
        *
        * @return number of elements stored in the GSS
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] size_t size(void) const noexcept { return dense_.size(); }
        
        /**
        * @brief Inserts a given element to the GSS
        *
        * Inserts a given element to the GSS. If the element has already been
        * inserted to the GSS, the function has no effect.
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain to be inserted
        *
        * @return true iff an actual insertion has been performed
        *
        * @par Complexity
        * O(1) since the dense_ vector always has enough capacity
        */ 
        bool insert(T val) {
            assert(val < domain_size_ && "Value outside of the GenerationSparseSet domain.");
            if(sparse_[val] == generation_) { return false; }
            sparse_[val] = generation_;
            dense_.push_back(val);
            return true;
        }
        
        /**
        * @brief Erases all elements from the GSS
        *
        * This function erases all elements from the GSS. It increases the current
        * generation index which invalidates all generation marks in the sparse_ vector. If the maximal
        * generation index is reached, the sparse_ vector will be cleared manually in linear time.
        *
        * @par Complexity
        * O(domain_size()) in an extremely rare situation when the maximal generation is reached,
        * O(1) otherwise.
        */ 
        void clear(void) {
            dense_.clear();
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
            dense_.reserve(size);
        }
        
        /**
        * @brief Returns a boolean flag which says whether the GSS is empty
        *
        * @return boolean flag corresponding to the emptiness of the GSS
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool empty(void) const noexcept { return dense_.empty(); }

        /**
        * @brief Returns a boolean flag which says whether the given element
        * is contained in the GSS
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        *
        * @return true iff val is stored in the GSS
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool contains(T val) const noexcept {
            assert(val < domain_size_ && "Value outside of the GenerationSparseSet domain.");
            return sparse_[val] == generation_;
        }

        /**
        * @brief Returns a boolean flag which says whether the given element
        * is contained in the GSS. This operator is equivalent to contains(val).
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain
        *
        * @return true iff val is stored in the GSS
        *
        * @par Complexity
        * O(1)
        */
        [[nodiscard]] bool operator[](T val) const noexcept {
            assert(val < domain_size_ && "Value outside of the GenerationSparseSet domain.");
            return contains(val);
        }

        /**
        * @brief Erases a single element from the GSS using linear 
        * iteration through the dense_ vector
        *
        * @pre val < domain_size()
        *
        * @param val value from the domain to be erased
        *
        * @warning For the sake of efficiency, this function should be used as little as possible
        * since it requires linear-time w.r.t. number of stored elements. The GSS should be
        * primarily used for purposes such that elements are cleared all at once. 
        *
        * @par Complexity
        * O(size())
        */
        void erase(T val) noexcept {
            assert(val < domain_size_ && "Value outside of the GenerationSparseSet domain.");
            if (sparse_[val] != generation_) { return; }
            sparse_[val] = 0;
            for (size_t i = 0; i < dense_.size(); ++i) {
                if (dense_[i] == val) {
                    dense_[i] = dense_.back();
                    dense_.pop_back();
                    break;
                }
            }
        }

        auto begin() noexcept { return dense_.begin(); }
        auto end() noexcept { return dense_.end(); }

        auto begin() const noexcept { return dense_.begin(); }
        auto end() const noexcept { return dense_.end(); }

        auto cbegin() const noexcept { return dense_.cbegin(); }
        auto cend() const noexcept { return dense_.cend(); }
};



#endif //LIBMATA_GENERATION_SPARSE_SET_HH
