#include "options.h"

void file_parser(std::ifstream &in, Options &opt)
{

    std::string appo;
    while (in >> appo)
    {
        if (appo == "BLOCKS")
        {
            in >> opt.n_block;
        }
        else if (appo == "STEPS")
        {
            in >> opt.n_step;
        }
        else if (appo == "DELTA")
        {
            in >> opt.delta;
        }
        else if (appo == "FIRST_ROW_INITIALIZER")
        {
            in >> opt.first_initializer_row;
        }
        else if (appo == "SECOND_ROW_INITIALIZER")
        {
            in >> opt.second_initializer_row;
        }
        else if (appo == "STARTING_POSITION")
        {
            in >> opt.start_pos.x >> opt.start_pos.y >> opt.start_pos.z;
        }
    }
}
