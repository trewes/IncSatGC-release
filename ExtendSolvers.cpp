#include "ExtendSolvers.h"

int CaDiCaLAdaptor::Solver::LitToVar(NSPACE::Lit &lit) {
    //go from Lit to literal that cadical can use
    bool sign = NSPACE::sign(lit); //sign or not
    int var = NSPACE::var(lit); //remove sign and divide by two
    var = var + 1; //cadical vars start at 1 instead of 0
    int new_var = sign ? -var : var;
    assert(new_var != 0);
    return new_var;
}

bool CaDiCaLAdaptor::Solver::addClause_(NSPACE::vec<NSPACE::Lit> &ps) {
    for(int i = 0; i < ps.size(); i++){
        NSPACE::Lit lit = ps[i];
        int new_var = LitToVar(lit);
        solver.add(new_var);
    }
    solver.add(0);
    num_clauses++;
    return true;
}

void CaDiCaLAdaptor::Solver::assume(const NSPACE::vec<NSPACE::Lit> &assumps) {
    //add assumptions from a vec of lits, used for the assumptions we get from the cardinality constraints
    for (int i = 0; i < assumps.size(); i++) {
        NSPACE::Lit lit = assumps[i];
        int new_var = LitToVar(lit);
        assert(new_var != 0);
        solver.assume(new_var);
    }
}

void CaDiCaLAdaptor::Solver::newVars(int num_vars) {
    //reserve space for new vars in cadical solver
    solver.reserve(solver.vars() + num_vars);
    //still add each var individually to base solver
    for (int i = 0; i < num_vars; ++i) {
        NSPACE::Solver::newVar();
    }
}


#ifdef CRYPTOMINISAT
std::vector<CMSat::Lit> CryptominisatAdaptor::Solver::LitsToCMLits(const NSPACE::vec<NSPACE::Lit>& ps) {
    std::vector<CMSat::Lit> CMLits;
    CMLits.reserve(ps.size());
    for (int i = 0; i < ps.size(); ++i) {
        NSPACE::Lit lit = ps[i];
        CMLits.emplace_back(NSPACE::var(lit), NSPACE::sign(lit));
    }
    return CMLits;
}

bool CryptominisatAdaptor::Solver::addClause_(NSPACE::vec<NSPACE::Lit> &ps) {
    num_clauses++;
    return solver.add_clause( LitsToCMLits(ps) );
}

NSPACE::Var CryptominisatAdaptor::Solver::newVar(bool polarity, bool dvar) {
    //add new var to both underlying cryptominisat solver
    //and base class glucose solver, again to ensure correct variable indices in cardinality encoding
    solver.new_var();
    return NSPACE::Solver::newVar(polarity, dvar);
}

void CryptominisatAdaptor::Solver::newVars(int num_vars) {
    //add vars all at once in cryptominisat solver
    solver.new_vars(num_vars);
    //still add each var individually to base solver
    for (int i = 0; i < num_vars; ++i) {
        NSPACE::Solver::newVar();
    }
}
#endif


