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
        std::cerr << "Could not understand which function was referenced, only values "
                     "that can be accepted are:\n"
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
        std::cerr << "Could not understand which function was referenced, only values "
                     "that can be accepted are:\n"
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
    nlohmann::json data = nlohmann::json::parse(in, nullptr, true, true);
    in.close();

    // Simulation Step Parameters
    data["blocks"].get_to(opt.n_block);
    data["block_steps"].get_to(opt.n_block_step);
    data["particle_steps"].get_to(opt.n_particle_step);

    // Movement Parameters
    data["delta"].get_to(opt.delta);

    // Random Number Generator Parameters
    data["first_row_initializer"].get_to(opt.first_initializer_row);
    data["second_row_initializer"].get_to(opt.second_initializer_row);

    // Starting Position
    data["starting_position"][0].get_to(opt.start_pos.x);
    data["starting_position"][1].get_to(opt.start_pos.y);
    data["starting_position"][2].get_to(opt.start_pos.z);

    // Thermalization Parameters
    data["thermalization_steps"].get_to(opt.n_thermalization_step);

    // Output Information Flags
    data["output_file"].get_to(opt.output_file);
    data["output_video"].get_to(opt.output_video);
    data["reset"].get_to(opt.reset);
    data["info_file"].get_to(opt.output_info);
    data["mean_radius_file"].get_to(opt.output_mean_radius);

    // Wave Function and Generator types
    std::string value = "";
    data["wave_function"].get_to(value);
    std::transform(value.begin(), value.end(), value.begin(), toupper);
    for (WAVE_FUNCTION wave_test = WAVE_FUNCTION::GROUND_STATE;
         wave_test != WAVE_FUNCTION::LAST;
         wave_test = static_cast<WAVE_FUNCTION>(static_cast<int>(wave_test) + 1))
    {
        if (value == wave_name(wave_test))
            opt.func = wave_test;
    }
    std::cout << "Wave Function analyzed is : " << wave_name(opt.func) << "\n";

    data["generator"].get_to(value);
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
