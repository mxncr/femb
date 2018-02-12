#pragma once

/* to disable exprtk warnings */
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wall"
#include "OGF/femb/third_party/exprtk.hpp"
#pragma GCC diagnostic pop

#include <string>

namespace femb {
    class ScalarFunction {
        public:
            ScalarFunction(const std::string& expr);
            double operator()(const double point[3]);

        protected:
            exprtk::expression<double>      expr_;
            exprtk::parser<double>          parser_;
            exprtk::symbol_table<double>    symbol_table_;
            double                          pt_[3];
    };
}

