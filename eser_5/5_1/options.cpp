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
                  << "\t- SECOND_EXCITED\n"
                  << "Aborting\n";
        exit(-1);
    }
}

std::string generator_name(NEW_POS_GENERATOR wave_enum)
{
    switch (wave_enum)
    {
    case NEW_POS_GENERATOR::UNIFORM:
        return std::string("UNIFORM");
    case NEW_POS_GENERATOR::NORMAL:
        return std::string("NORMAL");
    default:
        std::cerr << "Could not understand which function was referenced, only values that can be accepted are:\n"
                  << "\t- UNIFORM\n"
                  << "\t- NORMAL\n"
                  << "Aborting\n";
        exit(-1);
    }
}

bool compare_value(std::ifstream& in)
{
    int appo = 0;
    in >> appo;
    return (appo > 0);
}

void file_parser(std::ifstream& in, Options& opt)
{

    std::string appo;
    while (in >> appo)
    {
        if (appo == "BLOCKS")
        {
            in >> opt.n_block;
        }
        else if (appo == "BLOCK_STEPS")
        {
            in >> opt.n_block_step;
        }
        else if (appo == "PARTICLE_STEPS")
        {
            in >> opt.n_particle_step;
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
        else if (appo == "THERMALIZATION_STEPS")
        {
            in >> opt.n_thermalization_step;
        }
        else if (appo == "OUTPUT_FILE")
        {
            opt.output_file = compare_value(in);
        }
        else if (appo == "OUTPUT_VIDEO")
        {
            opt.output_video = compare_value(in);
        }
        else if (appo == "RESET")
        {
            opt.reset = compare_value(in);
        }
        else if (appo == "INFO_FILE")
        {
            opt.output_info = compare_value(in);
        }
        else if (appo == "MEAN_RADIUS_FILE")
        {
            opt.output_mean_radius = compare_value(in);
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
        else if (appo == "GENERATOR")
        {
            std::string value;
            in >> value;
            std::transform(value.begin(), value.end(), value.begin(), toupper);

            for (NEW_POS_GENERATOR gen_test = NEW_POS_GENERATOR::UNIFORM;
                 gen_test != NEW_POS_GENERATOR::LAST;
                 gen_test = static_cast<NEW_POS_GENERATOR>(static_cast<int>(gen_test) + 1))
            {
                if (value == generator_name(gen_test))
                    opt.new_pos_generator = gen_test;
            }
            std::cout << "Generator used is : " << generator_name(opt.new_pos_generator) << "\n";
        }
    }
}
