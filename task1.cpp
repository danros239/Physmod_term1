#include <fstream>
#include "json.hpp"
#include "function.hpp"

using json = nlohmann::json;

template<typename T, template <typename F_T> class Func, template <typename S_T> class State>
void run_test(Abstract_Solver<T, Func> s, State<T> &p, Func<T> &foo, 
    std::string filename,
    T dt,
    T time_limit,
    int algo)
{
    std::ofstream out(filename);
    out << "t x v" << std::endl;
    for(int i=0; i<1000; i++)
    {
        switch (s.type)
        {
        case algos::RK4:
            s.rk_solve(p, foo, dt);
            break;

        case algos::HEUN:
            s.heun_solve(p, foo, dt);
            break;

        case algos::EULER:
            s.euler_solve(p, foo, dt);
            break;

        default:
            break;
        
        }
        out << " " << p.t << " " << p.x[0] << " " << p.x[1] << std::endl;
    }
}

int main(int argc, char* argv[])
{
    std::ifstream file(argv[1]);
    json data = json::parse(file);
    file.close();

    std::cout << data.dump() << std::endl;
    std::string algorithm = data["algorithm"];

    // MPendulum<float> force(data["omega"]);
    MP_Friction<float> force(data["omega"], data["friction"]);
    MP_Periodic_Force<float> force_p(data["omega_0"], data["omega"], data["f"], data["friction"]);

    std::vector<float> state = {data["x0"], data["y0"]};
    NState<float> basic_pendulum(2, state, 0);


    Abstract_Solver<float, MP_Periodic_Force> as;
    if(data["algorithm"] == "rk4")
        as.type = algos::RK4;
    if(data["algorithm"] == "heun")
        as.type = algos::HEUN;
    if(data["algorithm"] == "euler")
        as.type = algos::EULER;

    std::string alg_type = data["algorithm"], func_type = data["name"];

    std::string filename = "data/" + alg_type + "_" + func_type + ".csv";

    run_test<float, MP_Periodic_Force, NState> (as, basic_pendulum, force_p, filename, data["step"], 1000, as.type);
    
}