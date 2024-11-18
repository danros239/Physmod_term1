#include <fstream>
#include "json.hpp"
#include "function.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    std::ifstream file(argv[1]);
    json data = json::parse(file);
    std::cout << data.dump() << std::endl;
    
    Abstract_Solver<float, MPendulum> as;

    MPendulum<float> force(1.0f);
    SingleState<float> basic_pendulum(data["x0"], data["y0"]);

    as.rk_solve(basic_pendulum, force, data["step"]);

    basic_pendulum.print();
    file.close();
}