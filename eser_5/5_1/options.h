#ifndef __OPTIONS__
#define __OPTIONS__

#include <fstream>
#include <functional>
#include <json.hpp>
#include <string>
#include <unordered_map>

#include "random.h"
#include "vector3D.h"
#include "wave_function.h"

typedef std::function<double(Vector3D)> wave_pdf;
typedef std::function<Vector3D(Random *, double)> generator_func;

enum class WAVE_FUNCTION : unsigned int
{
    GROUND_STATE = 0,
    FIRST_EXCITED,
    SECOND_EXCITED,
    LAST
};

enum class NEW_POS_GENERATOR : unsigned int
{
    UNIFORM = 0,
    NORMAL,
    LAST
};

std::string wave_name(WAVE_FUNCTION wave_enum);
std::string generator_name(NEW_POS_GENERATOR gen_enum);
WAVE_FUNCTION wave_enum(const std::string &wave_name);
NEW_POS_GENERATOR generator_enum(const std::string &gen_name);
bool compare_value(std::ifstream &in);

struct Options
{

    std::size_t n_block{1};               // Nr of Blocks
    std::size_t n_block_step{1};          // Nr of steps per block
    std::size_t n_particle_step{10};      // Nr of steps done by the particle between each measure
    std::size_t n_thermalization_step{0}; // Nr of steps done by the particle for
                                          // thermalization
    double delta{1.};                     // possible range of movement
    Vector3D start_pos{1., 0., 0.};       // Starting Position
    int first_initializer_row{0};         // Random Number generator
    int second_initializer_row{0};        // primes to use

    bool output_file{false};
    bool output_video{false};
    bool output_info{false};
    bool output_mean_radius{false};
    bool reset{false};

    WAVE_FUNCTION func{WAVE_FUNCTION::GROUND_STATE};                 // Wave Function to Sample
    NEW_POS_GENERATOR new_pos_generator{NEW_POS_GENERATOR::UNIFORM}; // Sampling Generator
    std::unordered_map<WAVE_FUNCTION, wave_pdf> convert = {
        {  WAVE_FUNCTION::GROUND_STATE,   ground_state},
        { WAVE_FUNCTION::FIRST_EXCITED,  first_excited},
        {WAVE_FUNCTION::SECOND_EXCITED, second_excited}
    };

    std::unordered_map<NEW_POS_GENERATOR, generator_func> generator = {
        {NEW_POS_GENERATOR::UNIFORM,  Vector3D::generate_unif},
        { NEW_POS_GENERATOR::NORMAL, Vector3D::generate_gauss}
    };
};

void file_parser(std::ifstream &in, Options &opt);
#endif // __OPTIONS__
