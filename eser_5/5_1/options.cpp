#include "options.h"

std::string wave_name(WAVE_FUNCTION wave_enum)
{
    switch (wave_enum)
    {
    case WAVE_FUNCTION::GROUND_STATE:
        return std::string("GROUND_STATE");
    case WAVE_FUNCTION::FIRST_EXCITED:
        return std::string("FIRST_EXCITED");
    case WAVE_FUNCTION::SECOND_EXCITED:
        return std::string("SECOND_EXCITED");
    default:
        std::cerr << "Could not understand which function was referenced, only values that can be accepted are:\n"
                  << "\t- GROUND_STATE\n"
                  << "\t- FIRST_EXCITED\n"
                  << "Aborting\n";
        exit(-1);
    }
}

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
        else if (appo == "OUTPUT_FILE")
        {
            int value = 0;
            in >> value;
            opt.output_file = (value > 0) ? true : false;
        }
        else if (appo == "OUTPUT_VIDEO")
        {
            int value = 0;
            in >> value;
            opt.output_video = (value > 0) ? true : false;
        }
        else if (appo == "WAVE_FUNCTION")
        {
            std::string value;
            in >> value;
            std::transform(value.begin(), value.end(), value.begin(), toupper);

            for (WAVE_FUNCTION wave_test = WAVE_FUNCTION::GROUND_STATE;
                 wave_test != WAVE_FUNCTION::LAST;
                 wave_test = static_cast<WAVE_FUNCTION>(static_cast<int>(wave_test) + 1))
            {
                if (value == wave_name(wave_test))
                    opt.func = wave_test;
            }
            std::cout << "Wave Function analyzed is : " << wave_name(opt.func) << "\n";
        }
    }
}
