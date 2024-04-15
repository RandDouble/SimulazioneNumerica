/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "system.h"
#include <iostream>

using namespace std;

/// @brief Function for parsing command line option, it can be extended
/// @param argc
/// @param argv
/// @return
/// @author Stefano Pilosio
/// @date 25th March 2024
bool check_write_to_file(int argc, char *argv[]);

int main(int argc, char *argv[])
{

    int nconf = 1;
    System SYS;
    SYS.initialize();
    SYS.initialize_properties();
    SYS.block_reset(0);
    bool write_to_file = check_write_to_file(argc, argv); // Bool checking id user want to wirte resulting config to file

    for (int i = 0; i < SYS.get_nbl(); i++) // loop over blocks
    {
        for (int j = 0; j < SYS.get_nsteps(); j++) // loop over steps in a block
        {
            SYS.step();
            SYS.measure();
            if (j % 10 == 0) // write every 10 step
            {
                if (write_to_file)
                {
                    SYS.write_XYZ(nconf); // Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                }
                nconf++;
            }
        }
        SYS.averages(i + 1); // Perform averages for blocks
        SYS.block_reset(i + 1);
    }
    SYS.finalize();

    return 0;
}

bool check_write_to_file(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0)
        {
            std::cout << "USAGE : " << argv[0] << " [OPTIONS... ]\n"
                      << "\t--write-config :  To output all the configuration assumed by system\n"
                      << "\t--help, -h     :  To display this help message\n";
            exit(0);
        }
        if (strcmp(argv[i], "--write-config") == 0)
        {
            return true;
        }
        std::cerr << argv[i] << " : Unknown Option\nIgnored\n";
    }
    return false;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
