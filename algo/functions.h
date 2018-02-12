#pragma once
#include "OGF/femb/third_party/exprtk.hpp"
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

