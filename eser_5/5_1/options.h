#ifndef __OPTIONS__
#define __OPTIONS__

#include <fstream>
#include <functional>
#include <string>
#include <unordered_map>

#include "vector3D.h"
#include "wave_function.h"

typedef std::function<double(Vector3D)> wave_pdf;

enum class WAVE_FUNCTION : unsigned int
{
    GROUND_STATE = 0,
    FIRST_EXCITED,
    SECOND_EXCITED,
    LAST
};

std::string wave_name(WAVE_FUNCTION wave_enum);

struct Options
{

public:
    std::size_t n_block{1};
    std::size_t n_step{10};
    double delta{1.};
    Vector3D start_pos{1., 0., 0.};
    int first_initializer_row{0}, second_initializer_row{0};
    bool output_file{false}, output_video{false};
    WAVE_FUNCTION func{WAVE_FUNCTION::GROUND_STATE};
    std::unordered_map<WAVE_FUNCTION, wave_pdf> convert = {
        {WAVE_FUNCTION::GROUND_STATE, ground_state},
        {WAVE_FUNCTION::FIRST_EXCITED, first_excited},
        {WAVE_FUNCTION::SECOND_EXCITED, second_excited}};


};

void file_parser(std::ifstream &in, Options &opt);

#endif // __OPTIONS__
