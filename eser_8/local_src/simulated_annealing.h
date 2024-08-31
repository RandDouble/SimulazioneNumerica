#include <iomanip>
#include <sstream>
#include <vector>

#include <utilities.h>

#include "metropolis.h"
#include "random.h"
#include "wave_function.h"

#ifndef __SIMULATED_ANNEALING__
#define __SIMULATED_ANNEALING__

class SimulatedAnnealing
{
private:
    // const double starting_temperature;
    // const double ending_temperature;
    const double integration_delta; // Montecarlo Movement

    const size_t n_temp_step;
    const size_t n_tempering_step;

    double m_actual_beta;

    PsiParam delta_move;
    const double m_cooling_rate;
    const PsiParam m_reduction_rate;

    size_t monte_carlo_step = 100;
    size_t monte_carlo_block = 100;

    MetropolisMemory<double> sampler;
    WaveFunction wave;

    std::stringstream param_history;

private:
    double boltzmann_weight(const double energy); // Calculate boltzmann weight
    void update_temperature();
    PsiParam propose_new_param();
    double scaling_movement_law() const;

public:
    SimulatedAnnealing(const double starting_temperature,
                       const double ending_temperature,
                       const size_t n_temp_step,
                       const size_t n_tempering_step,
                       const double integration_delta,
                       const PsiParam& delta_move,
                       const PsiParam& parameter_reduction_rate,
                       const PsiParam& wave_param) :
        integration_delta(integration_delta),
        n_temp_step(n_temp_step),
        n_tempering_step(n_tempering_step),
        m_actual_beta{(1. / starting_temperature)},
        delta_move(delta_move),
        // cooling_rate((starting_temperature - ending_temperature) / n_temp_step),
        m_cooling_rate(
            std::pow(ending_temperature / starting_temperature, 1. / n_temp_step)),
        m_reduction_rate(parameter_reduction_rate),
        sampler(0.),
        wave(wave_param)
    {
        assert(!std::isinf(starting_temperature) && "Starting temperature is infinite");
        assert(!std::isinf(ending_temperature) && "Ending temperature is infinite");
        assert(!std::isinf(m_actual_beta) && "Actual beta is infinite");
        assert(!std::isnan(m_actual_beta) && "Actual beta is NaN");
        assert(!std::isinf(m_cooling_rate) && "Cooling rate is infinite");
        assert(!std::isnan(m_cooling_rate) && "Cooling rate is NaN");
    }

    void set_metropolis_step(const size_t n_step) { monte_carlo_step = n_step; }
    void set_metropolis_block(const size_t n_block) { monte_carlo_block = n_block; }
    double get_monte_carlo_acceptance() const { return sampler.get_acceptance(); }
    PsiParam get_params() const { return wave.param(); }
    double get_temperature() const { return 1. / m_actual_beta; }
    double get_beta() const { return m_actual_beta; }

    void initialize_rng(size_t rows_to_skip = 0ull);

    void temperature_cycle();
    size_t tempering_cycle();

    values evaluate_energy(const double x_start);

    void write_results(const std::string& filename);
    void save_seed();
};

#endif // __SIMULATED_ANNEALING__
