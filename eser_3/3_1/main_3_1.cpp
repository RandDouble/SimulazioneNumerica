#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "finanza.h"
#include "random.h"
#include "utilities.h"

int main()
{
    Random rng;
    initializer(rng);

    constexpr std::size_t n_block = 100;
    constexpr unsigned int n_interval = 100;
    FinancialOption fnopt(&rng);

    std::ofstream f_out;

    std::vector v_C_final_value(n_block, 0.);
    std::vector v_P_final_value(n_block, 0.);
    std::vector v_C_step_value(n_block, 0.);
    std::vector v_P_step_value(n_block, 0.);

    // Progressive Mean Block Container
    std::vector C_final_mean(n_block, values{0., 0.});
    std::vector P_final_mean(n_block, values{0., 0.});
    std::vector C_step_mean(n_block, values{0., 0.});
    std::vector P_step_mean(n_block, values{0., 0.});

    auto C_final_begin = v_C_final_value.begin();
    auto P_final_begin = v_P_final_value.begin();
    auto C_step_begin = v_C_step_value.begin();
    auto P_step_begin = v_P_step_value.begin();

    for (std::size_t i = 0; i < n_block; ++i)
    {
        // One step iteration
        v_C_final_value[i] = fnopt.call_option_final_mean();
        v_P_final_value[i] = fnopt.put_option_final_mean();

        // Discretized array
        v_C_step_value[i] = fnopt.call_option_step_mean(n_interval);
        v_P_step_value[i] = fnopt.put_option_step_mean(n_interval);

        // Progressive Mean
        C_final_mean[i] = {calc_mean(C_final_begin, C_final_begin + i), calc_std(C_final_begin, C_final_begin + i)};
        P_final_mean[i] = {calc_mean(P_final_begin, P_final_begin + i), calc_std(P_final_begin, P_final_begin + i)};
        C_step_mean[i] = {calc_mean(C_step_begin, C_step_begin + i), calc_std(C_step_begin, C_step_begin + i)};
        P_step_mean[i] = {calc_mean(P_step_begin, P_step_begin + i), calc_std(P_step_begin, P_step_begin + i)};
    }

    // Output to file
    f_out.open("cost_option.csv");
    f_out << "Call_final,Call_final_err,Put_final,Put_final_err,Call_step,Call_step_err,Put_step,Put_step_err\n";

    for (std::size_t i = 0; i < n_block; i++)
    {
        f_out << C_final_mean[i] << ',' << P_final_mean[i] << ',' << C_step_mean[i] << ',' << P_step_mean[i] << "\n";
    }

    f_out.close();

    return 0;
}
