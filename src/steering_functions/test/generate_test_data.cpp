/**
 * Generate test data for Python-C++ consistency verification.
 *
 * Outputs a JSON file containing distances, controls and path endpoints
 * computed by all steering function state spaces for a fixed set of
 * random start/goal pairs.
 */
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "steering_functions/cc_dubins_state_space/cc00_dubins_state_space.hpp"
#include "steering_functions/cc_dubins_state_space/cc0pm_dubins_state_space.hpp"
#include "steering_functions/cc_dubins_state_space/cc_dubins_state_space.hpp"
#include "steering_functions/cc_dubins_state_space/ccpm0_dubins_state_space.hpp"
#include "steering_functions/cc_dubins_state_space/ccpmpm_dubins_state_space.hpp"
#include "steering_functions/cc_reeds_shepp_state_space/cc00_reeds_shepp_state_space.hpp"
#include "steering_functions/dubins_state_space/dubins_state_space.hpp"
#include "steering_functions/hc_reeds_shepp_state_space/hc00_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_reeds_shepp_state_space/hc0pm_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_reeds_shepp_state_space/hc_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_reeds_shepp_state_space/hcpm0_reeds_shepp_state_space.hpp"
#include "steering_functions/hc_reeds_shepp_state_space/hcpmpm_reeds_shepp_state_space.hpp"
#include "steering_functions/reeds_shepp_state_space/reeds_shepp_state_space.hpp"
#include "steering_functions/steering_functions.hpp"

using namespace std;
using namespace steering;

static constexpr double KAPPA          = 1.0;
static constexpr double SIGMA          = 1.0;
static constexpr double DISCRETIZATION = 0.1;
static constexpr int    NUM_SAMPLES    = 50;
static constexpr int    SEED           = 42;

/* deterministic random in [lower, upper] */
static double rnd(double lower, double upper)
{
    return rand() * (upper - lower) / RAND_MAX + lower;
}

static State random_state()
{
    State s;
    s.x     = rnd(-10.0, 10.0);
    s.y     = rnd(-10.0, 10.0);
    s.theta = rnd(-M_PI, M_PI);
    s.kappa = (rand() % 2) ? rnd(-KAPPA, KAPPA) : 0.0;
    s.d     = 0.0;
    return s;
}

/* Write a state as JSON object. */
static void write_state(ofstream& f, const State& s)
{
    f << "{\"x\":" << s.x
      << ",\"y\":" << s.y
      << ",\"theta\":" << s.theta
      << ",\"kappa\":" << s.kappa
      << ",\"d\":" << s.d << "}";
}

/* Write a control as JSON object. */
static void write_control(ofstream& f, const Control& c)
{
    f << "{\"delta_s\":" << c.delta_s
      << ",\"kappa\":" << c.kappa
      << ",\"sigma\":" << c.sigma << "}";
}

/* Write results for one state space. Returns whether any data was written. */
template<typename SS>
void write_ss_results(ofstream& f, const string& name, SS& ss,
                      const State& start, const State& goal, bool last)
{
    f << "        \"" << name << "\": {\n";

    double dist = ss.get_distance(start, goal);
    f << "          \"distance\": " << dist << ",\n";

    auto controls = ss.get_controls(start, goal);
    f << "          \"controls\": [";
    for (size_t ci = 0; ci < controls.size(); ++ci)
    {
        if (ci > 0) f << ",";
        write_control(f, controls[ci]);
    }
    f << "],\n";

    auto path = ss.get_path(start, goal);
    f << "          \"path_length\": " << path.size() << ",\n";
    f << "          \"path_end\": ";
    if (!path.empty())
        write_state(f, path.back());
    else
        f << "null";
    f << "\n";

    f << "        }";
    if (!last) f << ",";
    f << "\n";
}

int main(int argc, char** argv)
{
    string outpath = "test_data.json";
    if (argc > 1)
        outpath = argv[1];

    srand(SEED);

    // -------- state spaces --------
    CC_Dubins_State_Space          cc_dub_fwd(KAPPA, SIGMA, DISCRETIZATION, true);
    CC_Dubins_State_Space          cc_dub_bwd(KAPPA, SIGMA, DISCRETIZATION, false);
    CC00_Dubins_State_Space        cc00_dub_fwd(KAPPA, SIGMA, DISCRETIZATION, true);
    CC00_Dubins_State_Space        cc00_dub_bwd(KAPPA, SIGMA, DISCRETIZATION, false);
    CC0pm_Dubins_State_Space       cc0pm_dub_fwd(KAPPA, SIGMA, DISCRETIZATION, true);
    CC0pm_Dubins_State_Space       cc0pm_dub_bwd(KAPPA, SIGMA, DISCRETIZATION, false);
    CCpm0_Dubins_State_Space       ccpm0_dub_fwd(KAPPA, SIGMA, DISCRETIZATION, true);
    CCpm0_Dubins_State_Space       ccpm0_dub_bwd(KAPPA, SIGMA, DISCRETIZATION, false);
    CCpmpm_Dubins_State_Space      ccpmpm_dub_fwd(KAPPA, SIGMA, DISCRETIZATION, true);
    CCpmpm_Dubins_State_Space      ccpmpm_dub_bwd(KAPPA, SIGMA, DISCRETIZATION, false);
    Dubins_State_Space             dub_fwd(KAPPA, DISCRETIZATION, true);
    Dubins_State_Space             dub_bwd(KAPPA, DISCRETIZATION, false);
    Reeds_Shepp_State_Space        rs_ss(KAPPA, DISCRETIZATION);
    CC00_Reeds_Shepp_State_Space   cc00_rs(KAPPA, SIGMA, DISCRETIZATION);
    HC_Reeds_Shepp_State_Space     hc_rs(KAPPA, SIGMA, DISCRETIZATION);
    HC00_Reeds_Shepp_State_Space   hc00_rs(KAPPA, SIGMA, DISCRETIZATION);
    HC0pm_Reeds_Shepp_State_Space  hc0pm_rs(KAPPA, SIGMA, DISCRETIZATION);
    HCpm0_Reeds_Shepp_State_Space  hcpm0_rs(KAPPA, SIGMA, DISCRETIZATION);
    HCpmpm_Reeds_Shepp_State_Space hcpmpm_rs(KAPPA, SIGMA, DISCRETIZATION);

    // -------- generate samples --------
    vector<pair<State, State>> samples;
    samples.reserve(NUM_SAMPLES);
    for (int i = 0; i < NUM_SAMPLES; ++i)
    {
        State start = random_state();
        State goal;
        if (i == 0)
            goal = start;
        else if (i == 1)
        {
            goal       = start;
            goal.x    += 0.01;
            goal.y    += 0.01;
            goal.theta += 0.01;
        }
        else
            goal = random_state();
        samples.emplace_back(start, goal);
    }

    // -------- write JSON --------
    ofstream f(outpath);
    f << std::setprecision(17);
    f << "{\n";

    f << "  \"kappa\": " << KAPPA << ",\n";
    f << "  \"sigma\": " << SIGMA << ",\n";
    f << "  \"discretization\": " << DISCRETIZATION << ",\n";
    f << "  \"seed\": " << SEED << ",\n";
    f << "  \"num_samples\": " << NUM_SAMPLES << ",\n";

    f << "  \"samples\": [\n";
    for (size_t si = 0; si < samples.size(); ++si)
    {
        const auto& [start, goal] = samples[si];
        f << "    {\n";
        f << "      \"start\": "; write_state(f, start); f << ",\n";
        f << "      \"goal\": ";  write_state(f, goal);  f << ",\n";

        f << "      \"results\": {\n";
        write_ss_results(f, "dubins_forwards",        dub_fwd,       start, goal, false);
        write_ss_results(f, "dubins_backwards",       dub_bwd,       start, goal, false);
        write_ss_results(f, "reeds_shepp",            rs_ss,         start, goal, false);
        write_ss_results(f, "cc00_dubins_forwards",   cc00_dub_fwd,  start, goal, false);
        write_ss_results(f, "cc00_dubins_backwards",  cc00_dub_bwd,  start, goal, false);
        write_ss_results(f, "cc0pm_dubins_forwards",  cc0pm_dub_fwd, start, goal, false);
        write_ss_results(f, "cc0pm_dubins_backwards", cc0pm_dub_bwd, start, goal, false);
        write_ss_results(f, "ccpm0_dubins_forwards",  ccpm0_dub_fwd, start, goal, false);
        write_ss_results(f, "ccpm0_dubins_backwards", ccpm0_dub_bwd, start, goal, false);
        write_ss_results(f, "ccpmpm_dubins_forwards", ccpmpm_dub_fwd,start, goal, false);
        write_ss_results(f, "ccpmpm_dubins_backwards",ccpmpm_dub_bwd,start, goal, false);
        write_ss_results(f, "cc_dubins_forwards",     cc_dub_fwd,    start, goal, false);
        write_ss_results(f, "cc_dubins_backwards",    cc_dub_bwd,    start, goal, false);
        write_ss_results(f, "cc00_reeds_shepp",       cc00_rs,       start, goal, false);
        write_ss_results(f, "hc00_reeds_shepp",       hc00_rs,       start, goal, false);
        write_ss_results(f, "hc0pm_reeds_shepp",      hc0pm_rs,      start, goal, false);
        write_ss_results(f, "hcpm0_reeds_shepp",      hcpm0_rs,      start, goal, false);
        write_ss_results(f, "hcpmpm_reeds_shepp",     hcpmpm_rs,     start, goal, false);
        write_ss_results(f, "hc_reeds_shepp",         hc_rs,         start, goal, true);
        f << "      }\n";

        f << "    }";
        if (si + 1 < samples.size()) f << ",";
        f << "\n";
    }
    f << "  ]\n";
    f << "}\n";

    f.close();
    cerr << "Wrote " << outpath << " (" << samples.size() << " samples x 19 state spaces)\n";
    return 0;
}
