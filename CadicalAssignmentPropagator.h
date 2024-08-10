#ifndef INCSATGC_CADICALPROPAGATORDIRECT_H
#define INCSATGC_CADICALPROPAGATORDIRECT_H


#include <iostream>
#include <vector>
#include <queue>

#include <boost/dynamic_bitset.hpp>
#include "cadical.hpp"
#include "Graph.h"
#include "GraphMatrix.h"
#include "Statistics.h"

using ColorMap = std::vector<int>; //maps to each vertex the color it has been assigned, starts at color 0, -1 if uncolored
constexpr int NoColor = -1; //mark vertex unassigned in ColorMap
constexpr int NoVertex = -1; //mark vertex non-existing various data structures
using Bitset = boost::dynamic_bitset<>;

//declare IncSatGC
class IncSatGC;

class CadicalAssignmentPropagator : public CaDiCaL::ExternalPropagator{

public:
    bool is_lazy = false; // lazy propagator only checks complete assignments
    bool are_reasons_forgettable = true; // Reason external clauses can be deleted
    bool is_tainting = true; // The external clauses must trigger restore (unless frozen)

    ~CadicalAssignmentPropagator () override ;

    // Notify the propagator about assignments to observed variables.
    // The notification is not necessarily eager. It usually happens before
    // the call of propagator callbacks and when a driving clause is leading
    // to an assignment.
    //
    // void notify_assignment (int lit, bool is_fixed) override ;
    void notify_assignment (const std::vector<int>& lits) override;
    void notify_new_decision_level () override ;
    void notify_backtrack (size_t new_level) override ;

    // Check by the external propagator the found complete solution (after
    // solution reconstruction). If it returns false, the propagator must
    // provide an external clause during the next callback.
    //
    bool cb_check_found_model (const std::vector<int> &model) override ;

    // Ask the external propagator for the next decision literal. If it
    // returns 0, the solver makes its own choice.
    //
    int cb_decide () override ;

    // Ask the external propagator if there is an external propagation to make
    // under the current assignment. It returns either a literal to be
    // propagated or 0, indicating that there is no external propagation under
    // the current assignment.
    //
    int cb_propagate () override ;

    // Ask the external propagator for the reason clause of a previous
    // external propagation step (done by cb_propagate). The clause must be
    // added literal-by-literal closed with a 0. Further, the clause must
    // contain the propagated literal.
    //
    int cb_add_reason_clause_lit (int propagated_lit) override ;

    // The following two functions are used to add external clauses to the
    // solver during the CDCL loop. The external clause is added
    // literal-by-literal and learned by the solver as an irredundant
    // (original) input clause. The clause can be arbitrary, but if it is
    // root-satisfied or tautology, the solver will ignore it without learning
    // it. Root-falsified literals are eagerly removed from the clause.
    // Falsified clauses trigger conflict analysis, propagating clauses
    // trigger propagation. In case chrono is 0, the solver backtracks to
    // propagate the new literal on the right decision level, otherwise it
    // potentially will be an out-of-order assignment on the current level.
    // Unit clauses always (unless root-satisfied, see above) trigger
    // backtracking (independently from the value of the chrono option and
    // independently from being falsified or satisfied or unassigned) to level
    // 0. Empty clause (or root falsified clause, see above) makes the problem
    // unsat and stops the search immediately. A literal 0 must close the
    // clause.
    //
    // The external propagator indicates that there is a clause to add.
    // The parameter of the function allows the user to indicate that how
    // 'forgettable' is the external clause. Further, there is a member-parameter
    // called 'is_tainting' that indicates if the external clauses are irredundant
    // (almost always they are, see below the rare exceptions).
    // There are overall four possible scenarios defined by these two Boolean flags:
    // - Non-Forgettable & Tainting Clause [default value].
    //    - Policy to Forget: The clause will not be deleted, unless it is
    //      directly implied by some of the other irredundant clauses (e.g. root
    //      satisfied).
    //    - Tainting: The clause will be considered as irredundant and thus the
    //      literals of it will be tainted and might trigger restore steps
    //      (except frozen literals).
    //    - Example: It is always safe and correct to use this clause type.
    // - Forgettable & Tainting Clause.
    //    - Policy to Forget: The clause will be considered for deletion during
    //      clause database reduction rounds.
    //    - Tainting: The clause will be considered as irredundant and thus the
    //      literals of it will be tainted and might trigger restore steps (except
    //      frozen literals).
    //    - Example usage: Theory lemmas in SMT solvers. Anything that is
    //      irredundant from the SAT perspective (hence must be considered for
    //      tainting), but allowed to be deleted because the external propagator
    //      (e.g. theory solver) can always re-derive it.
    // - Forgettable & Not-Tainting Clause.
    //    - Policy to Forget: The clause will be considered for deletion during
    //      clause database reduction rounds.
    //    - Tainting: The clause is assumed to be derivable in propositional
    //      logic from the current set of irredundant and redundant clauses (i.e.
    //      (\phi \land \rho) \implies C, as in the side condition of Learn- in
    //      incremental inprocessing), therefore the literals of it are not
    //      tainted and do not trigger restore steps.
    //    - Example usage: Use it ONLY to share learned clauses between parallel
    //      SAT solvers. It is NOT correct to use this clause type unless the
    //      clause is indeed directly derivable by the SAT solver alone.
    // - Non-Forgettable & Non-Tainting Clause:
    //    - Policy to Forget: The clause will not be deleted, unless it is
    //      directly implied by some of the other irredundant clauses (e.g. root
    //      satisfied).
    //    - Redundancy: The clause is assumed to be derivable in propositional
    //      logic from the current set of irredundant and redundant clauses (i.e.
    //      (\phi \land \rho) \implies C, as in the side condition of Learn- in
    //      incremental inprocessing), therefore the literals of it are not
    //      tainted and do not trigger restore steps.
    //    - Example usage: Same as previous case, the only difference is that the
    //      clauses will be kept around (can be useful in cases where their
    //      derivation was expensive).
    //
    //
    // Reason clauses of external propagation steps are assumed to be forgettable,
    // and tainting by default. Use parameter 'reason_forgettable' to change it.
    //
    // DISCLAIMER: Any tainting clause is correct to be added ONLY if it is
    // guaranteed that the SAT solver uses only RUP-based clause addition
    // techniques (so no BVA) and the literals are either tainted or frozen.
    // See our Incremental Inprocessing SAT'19 paper for more details on it.
    //
    // The naming of the clause types will might change later (to coordinate
    // with IPASIR-2).
    //
    bool cb_has_external_clause (bool& is_forgettable) override;

    // The actual function called to add the external clause.
    //
    int cb_add_external_clause_lit () override ;


    /*
     * New section for new data members and functions
     */

    //need reference to IncSatGC instance to access data and call its functions, in particular the find conflict ones
    IncSatGC* INSTANCE; //use raw pointer to not take ownership of the superclass
    //constructor that constructs reference to INSTANCE
    explicit CadicalAssignmentPropagator(IncSatGC& reference, int full_num_colors_);

    Statistics& stats;
    const Options& options;

    Graph::NeighborList neighbors;
    //basic problem data
    int num_vertices;
    int full_num_colors;
    int num_colors;
    void update_num_colors(int num_colors_);
    //helper functions to go from vertex, color to variable index and vice versa
    [[nodiscard]] int xvi_to_lit(int v, int i) const;
    [[nodiscard]] std::pair<int, int> lit_to_xvi(int lit) const;

    int current_level = 0;
    std::deque<int> propagations;
    std::deque< std::vector<int> > clauses;
    //once a literal has been propagated, remember the reason clause in literal_to_clause, either positive or negative. mostly needed for clique based pruning
    std::vector<std::vector<int>> literal_to_clause_pos;
    std::vector<std::vector<int>> literal_to_clause_neg;
    //store the decisions and propagations done, i.e. the trail
    std::vector< std::vector<int> > current_trail;
    //store the assignment of variables so far, done with a signed char* of vals and negatively indexed as in cadical
    size_t vsize = 0; // actually allocated variable data size
    int max_var = 0;  // internal maximum variable index
    signed char *vals = nullptr; // assignment [-max_var,max_var]
    void enlarge_vals(int new_max_var);
    [[nodiscard]] signed char val(int lit) const;
    void set_val(int lit, signed char val);
    static int sign (int lit) { return (lit > 0) - (lit < 0); } //branchless sign computation
    int num_uncolored;

    //functions to do the assignment updates and propagation
    void assign_color(int lit);
    void unassign_color(int lit);

    void propagation_helper(int lit, int plit);
    void propagate_assignment(int lit);

    //functions to produce the next decision literal
    [[nodiscard]] int next_decision() const;
    [[nodiscard]] int num_uncolored_neigbors(int v) const;
    [[nodiscard]] int dsatur_decision() const;
    [[nodiscard]] int pass_same(int v1, int v2) const; //helper for pass vertex choice
    [[nodiscard]] int pass_vertex_decision(const std::vector<int> &candidates) const;
    [[nodiscard]] int largest_degree_decision() const;
    [[nodiscard]] int find_lowest_color(int selected_vertex) const;

    //data structures to keep track of aturation degree for branching decision
    //initialise all vertices to be uncolored and have zero saturation
    ColorMap vertex_color;
    //record how many colors are unavailable for a vertex, either by assignment of neighbours or by assignment of negative literal
    std::vector<int> forbidden_colors;

    [[nodiscard]] bool check_consisteny() const;

    // new section for integrating zykov approach functionality into assignment porpagator
    // i.e., an mgraph that is kept up to date according to the color asignments and functions to explain prunings
    MGraph mgraph;
    bool mgraph_is_used;
    MGraph original_graph; //needed for checking if edge already existed in original graph
    //vector to track representing vertex of a color class
    std::vector<int> color_rep;
    //vector in which to store any external clauses we want to add
    std::vector<std::vector<int>> external_clauses;
    //functions and fields to compute and store maximal cliques, often not of maximum size
    std::vector<Bitset> maximal_cliques;
    int max_clique_size;
    void compute_cliques();
    [[nodiscard]] bool need_to_recompute_cliques() const;
    //some helper functions to determine whether cliques need to be recomputed and which bounds we compute
    long long num_assignments_last_clique_computation = 0; //need to recompute cliques if assignments were done
    [[nodiscard]] bool compute_clique_clauses() const;
    [[nodiscard]] bool compute_mycielsky_clauses() const;
    bool first_call_after_backtrack; //only call expensive bounds after a backtrack

    //helper function to put subgraph explanation into a vector
    [[nodiscard]] std::vector<int> subgraph_explanation_clause(const Bitset &nodeset) const;
    //function that checks cliques for being of size > num_colors, and adds their explanation
    void check_for_clique_clauses();
    //given a clique, produce the explanation clause and add it to external_clauses
    void add_clique_explanation_clause(const Bitset &clique);
    //function that tries to compute mycielsky extension based off of cliques
    void check_for_mycielsky_clauses();
    //given a mycielsky bound, produce the explanation clause and add it to external_clauses
    void add_mycielsky_explanation_clause(const MGraph::SubGraph &subgraph);

    //functions for clique-based pruning as pre-processing of the subproblem graph
    // and to compute dominated vertices as a "good" decision to be made
    //store vertices whose neighbourhood changed, only those are relevant for neighbourhood search
    std::vector<int> touched_vertices;
    //only a decision, return the lit that was decided on
    int find_dominated_vertex_decisions();
    // void find_vertex_fusions();
    // void find_edge_additions();
    // void find_clique_based_pruning(); //need to fix clique based pruning

};



#endif //INCSATGC_CADICALPROPAGATORDIRECT_H
