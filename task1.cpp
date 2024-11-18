#include <fstream>
#include "json.hpp"
#include "function.hpp"

using json = nlohmann::json;

int main(int argc, char* argv[])
{
    std::ifstream file(argv[1]);
    json data = json::parse(file);
    file.close();

    std::cout << data.dump() << std::endl;
    std::string algorithm = data["algorithm"];
    
    Abstract_Solver<float, MPendulum> as;

    MPendulum<float> force(1.0f);
    SingleState<float> basic_pendulum(data["x0"], data["y0"]);

    std::ofstream output("euler.csv");
    output << "t x v" << std::endl;
    for(int i=0; i<1000; i++)
    {
        as.euler_solve(basic_pendulum, force, data["step"]);
        basic_pendulum.print();
        output << " " << basic_pendulum.t << " " << basic_pendulum.x << " " << basic_pendulum.v << std::endl;
    }
    output.close();


}