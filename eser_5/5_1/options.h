#ifndef __OPTIONS__
#define __OPTIONS__

#include <fstream>
#include <functional>
#include <string>
#include <unordered_map>

#include "random.h"
#include "vector3D.h"
#include "wave_function.h"

typedef std::function<double(Vector3D)> wave_pdf;
typedef std::function<Vector3D(Random*, double)> generator_func;

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
bool compare_value(std::ifstream& in);

struct Options
{

public:
    std::size_t n_block{1};
    std::size_t n_block_step{1};
    std::size_t n_particle_step{10};
    std::size_t n_thermalization_step{0};
    double delta{1.};
    Vector3D start_pos{1., 0., 0.};
    int first_initializer_row{0}, second_initializer_row{0};
    bool output_file{false}, output_video{false}, output_info{false}, output_mean_radius{false}, reset{false};
    WAVE_FUNCTION func{WAVE_FUNCTION::GROUND_STATE};
    NEW_POS_GENERATOR new_pos_generator{NEW_POS_GENERATOR::UNIFORM};
    std::unordered_map<WAVE_FUNCTION, wave_pdf> convert = {
        {WAVE_FUNCTION::GROUND_STATE, ground_state},
        {WAVE_FUNCTION::FIRST_EXCITED, first_excited},
        {WAVE_FUNCTION::SECOND_EXCITED, second_excited}};

    std::unordered_map<NEW_POS_GENERATOR, generator_func> generator = {
        {NEW_POS_GENERATOR::UNIFORM, Vector3D::generate_unif},
        {NEW_POS_GENERATOR::NORMAL, Vector3D::generate_gauss}
    };
};

void file_parser(std::ifstream& in, Options& opt);
#endif // __OPTIONS__
