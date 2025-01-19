#include <numeric>
#include <vector>

#define TEST_ENV

#include "initializer.h"
#include "mutation.h"
#include "random.h"
#include "swap_functions.h"
#include "utilities.h"

constexpr std::size_t SIZE = 10;

template <std::size_t SIZE = SIZE> bool test_swap_functions_swap_ranges()
{
    std::array<int, SIZE> test_arr;
    bool is_correct = true;
    print_vector<int> print;
    for (std::size_t i = 0; i < SIZE; i++)
    {
        std::iota(test_arr.begin(), test_arr.end(), 0);

        PBC_swap::swap_ranges<SIZE>(test_arr.begin(), i, i + 4, 3, 1);
        print(test_arr);
        // Checking if the first is untouched
        is_correct = is_correct && (test_arr[0] == 0);

        assert(is_correct && "Modified first element, something wrong");
    }
    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_swap_functions_rotate()
{
    std::array<int, SIZE> test_arr;
    bool is_correct = true;
    print_vector<int> print;
    for (std::size_t i = 0; i < SIZE; i++)
    {
        std::iota(test_arr.begin(), test_arr.end(), 0);

        PBC_swap::rotate<SIZE>(test_arr.begin(), i, i + 3, i + 5, 1);
        print(test_arr);
        // Checking if the first is untouched
        is_correct = is_correct && (test_arr[0] == 0);

        assert(is_correct && "Modified first element, something wrong");
    }
    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_swap_functions_reverse()
{
    std::array<int, SIZE> test_arr;
    bool is_correct = true;
    print_vector<int> print;
    for (std::size_t i = 0; i < SIZE; i++)
    {
        std::iota(test_arr.begin(), test_arr.end(), 0);

        PBC_swap::reverse<SIZE>(test_arr.begin(), i, i + SIZE - 1, 1);
        print(test_arr);
        // Checking if the first is untouched
        is_correct = is_correct && (test_arr[0] == 0);

        assert(is_correct && "Modified first element, something wrong");
    }
    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_PBC()
{
    print_vector<size_t> print_index;
    std::vector<size_t> indexes(30);
    std::iota(indexes.begin(), indexes.end(), 0);

    std::cout << "Indexes before PBC : \n";

    print_index(indexes, 2);

    std::for_each(indexes.begin(), indexes.end(), [&](size_t &el) { el = PBC(SIZE, el); });

    std::cout << "Indexes after PBC : \n";

    print_index(indexes, 2);

    std::cout << '\n';

    bool is_correct = true;

    for (auto &idx : indexes)
    {
        is_correct = is_correct && idx > 0 && idx < (SIZE + 1);
    }
    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_pair_permutation()
{
    bool is_correct = true;
    Random rng;
    initializer(rng);

    Individual<SIZE> cities;

    std::cout << "Positions Before Pair Permutation :\n";
    cities.print_DNA();

    cities.pair_permutation(rng);

    std::cout << "Positions After Pair Permutation :\n";

    cities.print_DNA();

    std::cout << '\n' << "Doing more shuffling, checking for duplicates and if touching first\n\n";

    for (int i = 0; i < 100; i++)
    {
        cities.pair_permutation(rng);
    }

    // Checking if the first is untouched
    is_correct = is_correct && (cities[0] == 0);

    assert(is_correct && "Modified first element, something wrong");

    // Checking if there are multiple elements
    std::sort(cities.begin(), cities.end());
    is_correct = is_correct && (cities.end() == std::unique(cities.begin(), cities.end()));

    assert(is_correct && "There are duplicates, something wrong");

    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_shift_block()
{
    bool is_correct = true;
    Random rng;
    initializer(rng);

    /// @todo : Aggiungere test per shifting
    /// @todo : Terminare implementazine shift block

    Individual<SIZE> cities;

    std::cout << "Positions Before Shift Block :\n";
    cities.print_DNA();

    cities.shift_block(rng);

    std::cout << "Positions After Shift Block :\n";

    cities.print_DNA();

    std::cout << '\n' << "Doing more shuffling, checking for duplicates and if touching first\n";

    for (int i = 0; i < 100; i++)
    {
        // std::cout << i << '\n';
        cities.shift_block(rng);
    }

    // Checking if the first is untouched
    is_correct = is_correct && (cities[0] == 0);

    assert(is_correct && "Modified first element, something wrong");

    // Checking if there are multiple elements
    std::sort(cities.begin(), cities.end());
    is_correct = is_correct && (cities.end() == std::unique(cities.begin(), cities.end()));

    assert(is_correct && "There are duplicates, something wrong");

    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_permutate_contiguos()
{

    bool is_correct = true;
    Random rng;
    initializer(rng);

    // Initializing cities vector
    Individual<SIZE> cities;

    std::cout << "Positions Before Permutate Contiguos :\n";
    cities.print_DNA();

    cities.permutate_contiguos(rng);

    std::cout << "Positions After Permutate Contiguos :\n";

    cities.print_DNA();

    std::cout << '\n' << "Doing more shuffling, checking for duplicates and if touching first\n";

    for (int i = 0; i < 100; i++)
    {
        cities.permutate_contiguos(rng);
    }

    // Checking if the first is untouched
    is_correct = is_correct && (cities[0] == 0);

    assert(is_correct && "Modified first element, something wrong");

    // Checking if there are multiple elements
    std::sort(cities.begin(), cities.end());
    is_correct = is_correct && (cities.end() == std::unique(cities.begin(), cities.end()));

    assert(is_correct && "There are duplicates, something wrong");

    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_inversion()
{
    bool is_correct = true;
    Random rng;
    initializer(rng);

    // Initializing cities vector
    Individual<SIZE> cities;

    std::cout << "Positions Before Inversion :\n";
    cities.print_DNA();

    cities.inversion(rng);

    std::cout << "Positions After Inversion :\n";

    cities.print_DNA();

    std::cout << '\n' << "Doing more inversions, checking for duplicates and if touching first\n";

    for (int i = 0; i < 100; i++)
    {
        cities.inversion(rng);
    }

    // Checking if the first is untouched
    is_correct = is_correct && (cities[0] == 0);

    assert(is_correct && "Modified first element, something wrong");

    // Checking if there are multiple elements
    std::sort(cities.begin(), cities.end());
    is_correct = is_correct && (cities.end() == std::unique(cities.begin(), cities.end()));

    assert(is_correct && "There are duplicates, something wrong");

    return is_correct;
}

template <std::size_t SIZE = SIZE> bool test_crossover()
{
    Random rng;
    initializer(rng);

    std::cout << "Starting CrossOver Test\n";
    Individual<SIZE> father = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    Individual<SIZE> mother = {0, 9, 8, 7, 6, 5, 4, 3, 2, 1};
    Individual<SIZE> son, daughter;

    bool is_correct = true;

    std::cout << "Father is :\n";
    father.print_DNA();

    std::cout << "Mother is :\n";
    mother.print_DNA();

    father.crossover(mother, daughter, son, rng);

    std::cout << "Son is :\n";
    son.print_DNA();

    std::cout << "Daughter is :\n";
    daughter.print_DNA();

    // Checking if the first is untouched
    std::cout << "Checking Father\n";
    is_correct = is_correct && father.check_health();
    if (!is_correct)
    {
        std::cout << "Father is BAD\n";
        father.print_DNA();
    }

    std::cout << "Checking Mother\n";
    is_correct = is_correct && mother.check_health();
    if (!is_correct)
    {
        std::cout << "Mother is BAD\n";
        mother.print_DNA();
    }

    std::cout << "Checking Son\n";
    is_correct = is_correct && son.check_health();
    if (!is_correct)
    {
        std::cout << "Son is BAD\n";
        son.print_DNA();
    }
    std::cout << "Checking Daughter\n";
    is_correct = is_correct && daughter.check_health();
    if (!is_correct)
    {
        std::cout << "Daughter is BAD\n";
        daughter.print_DNA();
    }
    return is_correct;
}

int main()
{
#ifdef NDEBUG
    std::cout << "Active NDEBUG, test results may not be complete\n";
#endif
    if (test_swap_functions_swap_ranges())
    {
        std::cout << "SWAP RANGES TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "SWAP RANGES TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_swap_functions_rotate())
    {
        std::cout << "ROTATE TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "ROTATE TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_swap_functions_reverse())
    {
        std::cout << "REVERSE TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "REVERSE TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_PBC())
    {
        std::cout << "PBC TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "PBC TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_pair_permutation())
    {
        std::cout << "Pair Permutation TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "Pair Permutation TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_permutate_contiguos())
    {
        std::cout << "Permutate Contiguos TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "Permutate Contiguos TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_shift_block())
    {
        std::cout << "Shift Block TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "Shift Block TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_inversion())
    {
        std::cout << "Inversion TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "Inversion TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    if (test_crossover())
    {
        std::cout << "Crossover TEST : PASSED\n\n";
    }
    else
    {
        std::cout << "Crossover TEST : FAILED : ERROR\n\n";
        exit(-1);
    }

    return 0;
}