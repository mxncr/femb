#pragma once

#include <string>

namespace GEO {
    class Mesh;
}

namespace femb {

    bool fem_simulation(GEO::Mesh& M,
            const std::string& dirichlet_region,
            const std::string& dirichlet_value,
            const std::string& neumann_region,
            const std::string& neumann_value,
            const std::string& sourceterm_value,
            const std::string& diffusion_value,
            const std::string& solution_name
            );

}
