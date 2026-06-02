#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "mata/utils/sparse-multiset.hh"

TEST_CASE("mata::utils::SparseMultiset") {
    SECTION("Initialization") {
        SparseMultiset<size_t> sm;
        CHECK(sm.domain_size() == 1);
        CHECK(sm.size() == 0);
        
        SparseMultiset<size_t> sm2 = SparseMultiset<size_t>(10);
        CHECK(sm2.domain_size() == 10);
        CHECK(sm2.size() == 0);

        SparseMultiset<size_t> sm3 = SparseMultiset<size_t>({0, 1, 8}, 10);
        CHECK(sm3.domain_size() == 10);
        CHECK(sm3.size() == 3);
    }

    SECTION("Insertion") {
        SparseMultiset<size_t> sm = SparseMultiset<size_t>(20);
        sm.insert(0);
        sm.insert(1);
        sm.insert(2);
        sm.insert(18);
        CHECK(sm.domain_size() == 20);
        CHECK(sm.size() == 4);
        CHECK(sm.occurrences(0) == 1);
        CHECK(sm.occurrences(1) == 1);
        CHECK(sm.occurrences(2) == 1);
        CHECK(sm.occurrences(16) == 0);
        CHECK(sm.occurrences(18) == 1);
        sm.insert(18);
        sm.insert(18);
        sm.insert(18);
        sm.insert(18);
        CHECK(sm.size() == 4);
        CHECK(sm.occurrences(0) == 1);
        CHECK(sm.occurrences(1) == 1);
        CHECK(sm.occurrences(2) == 1);
        CHECK(sm.occurrences(16) == 0);
        CHECK(sm.occurrences(18) == 5);
        CHECK(sm.insert(6) == 1);
        CHECK(sm.insert(6) == 2);
        CHECK(sm.insert(6) == 3);
        CHECK(sm.insert(6) == 4);
        CHECK(sm.insert(6) == 5);
        CHECK(sm.occurrences(6) == 5);
        size_t domain = sm.domain_size();
        for(size_t i = 0; i < domain; ++i) {
            CHECK(sm[i] == sm.occurrences(i));
        }
        sm.clear();
        CHECK(sm.size() == 0);
        for(size_t i = 0; i < domain; ++i) {
            CHECK(!sm[i]);
        }
    }

}
