#ifndef INCSATGC_EXTENDSOLVERS_H
#define INCSATGC_EXTENDSOLVERS_H


#include "core/Solver.h"
#include "cadical.hpp"
#ifdef CRYPTOMINISAT
    #include "cryptominisat.h"
#endif

//idea of this file is to adapt certain functions to use other solvers,
//so that the cardinality encodings can work with them
//for that, mainly override the addClause_ function from NSPACE::Solver
//and for cryptominisat also newVar() to add new vars when done so in encoding

namespace CaDiCaLAdaptor {
    class Solver : public NSPACE::Solver{
    public:
        //go from Lit to literal that cadical can use
        static int LitToVar(NSPACE::Lit &lit);
        //overrides addClause_ so that the clauses are added to underlying cadical solver
        bool addClause_(NSPACE::vec<NSPACE::Lit>& ps) override;
        void assume(const NSPACE::vec<NSPACE::Lit> &assumps);
        //no need to override newVar as cadical adds new vars on the fly as needed

        //make new function to add many variables at once
        void newVars(int num_vars);

        int num_clauses = 0;
        CaDiCaL::Solver solver;
    };

}

#ifdef CRYPTOMINISAT
namespace CryptominisatAdaptor {
    class Solver : public NSPACE::Solver{
    public:
        static std::vector<CMSat::Lit> LitsToCMLits(const NSPACE::vec<NSPACE::Lit>& ps);
        //overrides addClause_ so that the clauses are added to underlying cadical solver
        bool addClause_(NSPACE::vec<NSPACE::Lit>& ps) override;
        //override to add new var to both base solver and cryptominisat solver
        NSPACE::Var newVar(bool polarity = true, bool dvar = true) override;
        //make new function to add many variables at once
        void newVars(int num_vars);

        int num_clauses = 0;
        CMSat::SATSolver solver;
    };

}
#endif


#endif //INCSATGC_EXTENDSOLVERS_H
