#include "OGF/femb/algo/functions.h"

#include <iostream>

namespace femb {
    ScalarFunction::ScalarFunction(const std::string& expr){
        symbol_table_.add_variable("x", pt_[0]);
        symbol_table_.add_variable("y", pt_[1]);
        symbol_table_.add_variable("z", pt_[2]);
        symbol_table_.add_constants();
        expr_.register_symbol_table(symbol_table_);

        if(!parser_.compile(expr,expr_) ) {
            std::cout << "[exprtk parsing] Wrong expression: \n "
                << parser_.error() << std::endl;;
            return;
        }
    }

    double ScalarFunction::operator()(const double point[3]){
        pt_[0] = point[0];
        pt_[1] = point[1];
        pt_[2] = point[2];
        return expr_.value();
    }
}
