#include "CadicalAssignmentPropagator.h"

#include "IncSatGC.h" //make type information available to use here for pointer of IncSatGC object

# define  DEBUG_PRINT(X) //do {X} while(0)  //debug print that can be disabled easily

# define  PROP_TIMING(X) //do {X} while(0)  //call to timing functions that can easily be disabled

CadicalAssignmentPropagator::~CadicalAssignmentPropagator() = default;


void CadicalAssignmentPropagator::notify_assignment(const std::vector<int>& lits) {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorAssignment););
    stats.prop_num_assignments += static_cast<int>(lits.size());
    DEBUG_PRINT(std::cout << "Was notified of assignments " << lits << " at level " << current_level << "\n";);
    for(int lit : lits) {
        DEBUG_PRINT(
            int v = std::get<0>(lit_to_xvi(lit));
            int i = std::get<1>(lit_to_xvi(lit));
            std::cout << "Was notified of assignment of literal " << lit << " = x_" << v << "," << i
            << " at level " << current_level << "\n";
        );

        if(val(lit) != 0){
            DEBUG_PRINT(std::cout << "propagated lit " << lit << " already has an assignment " << val(lit) << "\n";);
            //might get notified about the same assignment twice
            assert(val(lit) == (lit > 0 ? 1 : -1));
            continue;
        }
        if(lit > 0 and vertex_color[std::get<0>(lit_to_xvi(lit))] != NoColor) {
            //already assigned a color, should cause a conflict and not be communicated to assignments and graph
            assert(std::find(propagations.begin(), propagations.end(), -lit) != propagations.end()
                or std::find(propagations.begin(), propagations.end(), -xvi_to_lit(std::get<0>(lit_to_xvi(lit)), vertex_color[std::get<0>(lit_to_xvi(lit))])) != propagations.end());
            DEBUG_PRINT(std::cout << "skipped " << lit << " as vertex was already assigned a color\n";);
            continue;
        }

        set_val(lit, sign (lit));
        current_trail.back().push_back(lit);


        //assignment of a color to a vertex,
        //update saturation and vertex colors, but only if this is the first assignemnt, i.e. vertex is still uncolored
        assign_color(lit);

        DEBUG_PRINT(std::cout << "Forbidden colors " << forbidden_colors << "\n";);
        DEBUG_PRINT(std::cout << "Vertex colors " << vertex_color << "\n";);
        DEBUG_PRINT(std::cout << "Color reps " << color_rep << "\n";);

        //handle x_v,i = true and x_v,i = false on a propagation level
        propagate_assignment(lit);
    }
    assert(check_consisteny());
    PROP_TIMING(stats.end_phase(Statistics::PropagatorAssignment););
}

void CadicalAssignmentPropagator::notify_new_decision_level() {
    assert(current_level + 1 == current_trail.size());
    current_level++;
    DEBUG_PRINT(std::cout << "Was notified of new decision level " << new_level << "\n";);
    //assert there are no more propagations
    assert(propagations.empty());
    assert(clauses.empty());

    current_trail.emplace_back();
    assert(check_consisteny());
    if(mgraph_is_used) {
        mgraph.notify_new_level();
        assert(mgraph.added_edges_limiter.size() == mgraph.current_level);
    }
    stats.prop_max_level = std::max(stats.prop_max_level, current_level);
    touched_vertices.clear();
    //pre-emptive resize of stats vector for larger levels
    if(current_level >= stats.prop_node_depth_history.size()){
        int new_size = current_level * 2;
        stats.prop_node_depth_history.resize(new_size, 0);
        if(options.use_clique_explanation_clauses){
            stats.prop_clique_pruning_level.resize(new_size,0);
        }
        if(options.use_mycielsky_explanation_clauses){
            stats.prop_myc_pruning_level.resize(new_size,0);
        }
        if(options.enable_positive_pruning){
            stats.prop_positive_pruning_level.resize(new_size,0);
        }
        if(options.enable_negative_pruning){
            stats.prop_negative_pruning_level.resize(new_size,0);
        }
    }
    stats.prop_node_depth_history[current_level]++;
}

void CadicalAssignmentPropagator::notify_backtrack(size_t new_level) {
    DEBUG_PRINT(std::cout << "Was notified of backtrack: " << new_level << "\n";);
    PROP_TIMING(stats.start_phase(Statistics::PropagatorBacktrack););
    stats.prop_num_backtracks++;
    int length = current_level - new_level;
    if(length >= stats.prop_backtrack_size.size()) {
        stats.prop_backtrack_size.resize(2 * length, 0);
    }
    stats.prop_backtrack_size[length]++;
    if(options.enable_detailed_backtracking_stats) { //collect current level, level we backtrack to, and reason
        stats.prop_detailed_backtrack_list.insert(stats.prop_detailed_backtrack_list.end(), {current_level, new_level, stats.backtrack_reson});
        stats.backtrack_reson = 0;
    }

    assert(current_level > new_level);
    while (current_level > new_level){
        for(auto l_it = current_trail.back().rbegin(); l_it != current_trail.back().rend(); ++l_it){
            int lit = *l_it;
            set_val(lit,0);
            unassign_color(lit);
        }
        current_trail.pop_back();
        current_level--;
        assert(current_level + 1 == current_trail.size());
    }

    if(mgraph_is_used) {
        mgraph.notify_backtrack_level(static_cast<int>(new_level));
        assert(mgraph.added_edges_limiter.size() == mgraph.current_level);
        DEBUG_PRINT(mgraph.print(););
    }

    //reset all propagation info
    propagations.clear();
    clauses.clear();
    first_call_after_backtrack = true;
    touched_vertices.clear();

    assert(check_consisteny());
    PROP_TIMING(stats.end_phase(Statistics::PropagatorBacktrack););
}

bool CadicalAssignmentPropagator::cb_check_found_model(const std::vector<int> &model) {
    DEBUG_PRINT(std::cout << "EP checked found model\n";);
    ColorMap colors(num_vertices, NoColor);
    for (int v = 0; v < num_vertices; v++) {
        for (int c = 0; c < num_colors; c++) {
            int index = full_num_colors * v + c;
            if(model[index] > 0){ //vertex i gets color j
                colors[v] = c;
            }
        }
    }
    bool is_valid_model = INSTANCE->is_valid_coloring(colors);
    assert(is_valid_model);
    return is_valid_model;
}

int CadicalAssignmentPropagator::cb_decide() {
    DEBUG_PRINT(std::cout << "EP was asked for next decision literal\n";);
    PROP_TIMING(stats.start_phase(Statistics::PropagatorDecide););
    int decision_literal = 0;
    if(options.use_dominated_vertex_decisions) {
        //returns the literal to contract two dominated vertices or 0 if none was found
        decision_literal = find_dominated_vertex_decisions();
    }
    if(decision_literal == 0) {
        decision_literal = next_decision();
    }
    stats.prop_num_decisions++;
    PROP_TIMING(stats.end_phase(Statistics::PropagatorDecide););
    return decision_literal;
}

int CadicalAssignmentPropagator::cb_propagate() {
    DEBUG_PRINT(std::cout << "EP was asked if there is an external propagation\n";);
    assert(clauses.size() == propagations.size());
    if(propagations.empty() and (options.enable_positive_pruning or options.enable_negative_pruning)) {
        //check if there are prunings to be propagated
        // find_clique_based_pruning();
    }
    while (not propagations.empty()){
        //get next propagation and store corresponding reason clause in literal_to_clause
        int literal = propagations.front();
        int abslit = std::abs(literal);
        propagations.pop_front();
        std::vector<int> reason_clause = clauses.front();
        clauses.pop_front(); // delete last clause
        if(val(literal) == sign(literal)) {
            continue; //literal was already assigned to value of propagation, take next propagation
        }
        //make sure the literal is part of the given clause
        assert(std::find(reason_clause.begin(), reason_clause.end(), literal) != reason_clause.end());
        if(literal > 0) {
            literal_to_clause_pos[abslit] = reason_clause;
        }
        else {
            literal_to_clause_neg[abslit] = reason_clause;
        }

        DEBUG_PRINT(
            int v = std::get<0>(lit_to_xvi(literal));
            int i = std::get<1>(lit_to_xvi(literal));
            std::cout << "cb_propagate " << literal << " (x_" << v << "," << i << ")\n";
        );
        stats.prop_num_propagations++;
        return literal;
    }
    assert(propagations.empty() and clauses.empty());
    return 0;
}

int CadicalAssignmentPropagator::cb_add_reason_clause_lit(int propagated_lit) {
    DEBUG_PRINT(std::cout << "EP was asked for the next lit reason of external propagation " << propagated_lit << "\n";);
    //given a previously propagated literal, produce the reason literal by literal now
    int abslit = std::abs(propagated_lit);
    auto &reason = (propagated_lit > 0 ? literal_to_clause_pos[abslit] : literal_to_clause_neg[abslit]);
    if (not reason.empty()){
        auto lit = reason.back();
        reason.pop_back();
        DEBUG_PRINT(std::cout << "cb_add_reason_clause_lit " << lit << "\n";);
        return lit;
    }
    //gave reason clause, increment stat
    stats.prop_num_reason_clauses++;
    return 0;
}

bool CadicalAssignmentPropagator::cb_has_external_clause(bool& is_forgettable) {
    DEBUG_PRINT(std::cout << "EP was asked if there are external clauses\n";);
    //only look for external clauses if there are no more propagations to be made to keep assignemnts and graph in sync
    //also only call after backtracks since assignment propagator goes through too many nodes to compute a clique each time
    if(not propagations.empty() or not first_call_after_backtrack) {
        return false;
    }
    if(compute_clique_clauses()) {
        check_for_clique_clauses();
    }
    if(compute_mycielsky_clauses()) {
        check_for_mycielsky_clauses();
    }
    first_call_after_backtrack = false;

    if(not external_clauses.empty()){
        DEBUG_PRINT(std::cout << "cb_has_external_clause : true, because external_clauses not empty\n";);
        is_forgettable = true; //the explanation clauses can be forgotten, we only really want to use them to backtrack once
        return true;
    }
    return false;
}

int CadicalAssignmentPropagator::cb_add_external_clause_lit() {
    DEBUG_PRINT(std::cout << "EP was asked for the next external clause literal\n";);
    std::vector<int> &clause = external_clauses.back();
    if (not clause.empty()){
        int lit = clause.back();
        clause.pop_back();
        DEBUG_PRINT(std::cout << "cb_add_external_clause_lit " << lit << "\n";);
        return lit;
    }
    // delete last clause
    external_clauses.pop_back();
    stats.prop_num_external_clauses++;
    return 0;
}


/*
 * New section for new data members and functions
 */

CadicalAssignmentPropagator::CadicalAssignmentPropagator(IncSatGC &reference, int full_num_colors_)
    : stats(reference.stats), options(reference.options)
{
    INSTANCE = &reference;
    neighbors = reference.graph.get_neighbor_list();
    num_vertices = INSTANCE->num_vertices;
    full_num_colors = full_num_colors_;
    num_colors = full_num_colors;

    propagations = {};
    clauses = {};
    literal_to_clause_pos = std::vector<std::vector<int>>(full_num_colors * num_vertices + 1);
    literal_to_clause_neg = std::vector<std::vector<int>>(full_num_colors * num_vertices + 1);
    current_trail.emplace_back();
    enlarge_vals(full_num_colors * num_vertices + 1);
    num_uncolored = num_vertices;

    //branching decision data structures
    vertex_color.resize(num_vertices, NoColor);
    forbidden_colors.resize(num_vertices, 0);

    mgraph_is_used = (options.use_clique_explanation_clauses or options.use_dominated_vertex_decisions);
    if(mgraph_is_used) {//only keep reduced graph active if it is used for pruning
        mgraph = MGraph(num_vertices, INSTANCE->graph.ecount(), INSTANCE->graph.elist());
    }
    original_graph = MGraph(num_vertices, INSTANCE->graph.ecount(), INSTANCE->graph.elist());
    color_rep.resize(full_num_colors, NoVertex);
    first_call_after_backtrack = false;

    //pre-allocate space for some of the stat vectors
    stats.prop_node_depth_history.resize(500,0);
    stats.prop_backtrack_size.resize(50,0);
    if(options.enable_detailed_backtracking_stats) {
        stats.prop_detailed_backtrack_list.reserve(100000);
    }
    if(options.use_clique_explanation_clauses){
        stats.prop_clique_pruning_level.resize(500,0);
    }
    if(options.use_mycielsky_explanation_clauses){
        stats.prop_myc_pruning_level.resize(500,0);
    }
    if(options.enable_positive_pruning){
        stats.prop_positive_pruning_level.resize(500,0);
    }
    if(options.enable_negative_pruning){
        stats.prop_negative_pruning_level.resize(500,0);
    }
}

void CadicalAssignmentPropagator::update_num_colors(const int num_colors_) {
    num_colors = num_colors_;
}

int CadicalAssignmentPropagator::xvi_to_lit(const int v, const int i) const {
    return full_num_colors * v + i + 1;
}

std::pair<int, int> CadicalAssignmentPropagator::lit_to_xvi(int lit) const {
    lit = std::abs(lit);
    return {(lit-1) / full_num_colors, (lit-1) % full_num_colors};
}


void CadicalAssignmentPropagator::enlarge_vals(const int new_max_var) {
    size_t new_vsize = 2 * new_max_var;
    signed char *new_vals;
    const size_t bytes = 2u * new_vsize;
    new_vals = new signed char[bytes]; // g++-4.8 does not like ... { 0 };
    memset (new_vals, 0, bytes);
    new_vals += new_vsize;

    if (vals) {
        memcpy (new_vals - max_var, vals - max_var, 2u * max_var + 1u);
        vals -= vsize;
        delete[] vals;
    }
    vals = new_vals;
    vsize = new_vsize;
    max_var = new_max_var;
}

signed char CadicalAssignmentPropagator::val(const int lit) const {
    assert (-max_var <= lit);
    assert (lit);
    assert (lit <= max_var);
    return vals[lit];
}

void CadicalAssignmentPropagator::set_val(const int lit, const signed char val) {
    assert (-1 <= val);
    assert (val <= 1);
    assert (-max_var <= lit);
    assert (lit);
    assert (lit <= max_var);
    vals[lit] = val;
    vals[-lit] = val; //cadical assigns -val and then takes care when retreiving the value, we want the assignment regardless of sign!
}

void CadicalAssignmentPropagator::assign_color(const int lit) {
    auto [v,i] = lit_to_xvi(lit);
    if(lit < 0) {
        //assignment of negative literal means that a color was forbidden
        forbidden_colors[v]++;
        if(not mgraph_is_used) {
            return;
        }
        //if color rep exists, add an edge from vertex to rep of forbidden color
        if(color_rep[i] != NoVertex ) {
            DEBUG_PRINT(std::cout << "Added edge from " << v << " with color !=" << i << " to "
                        << mgraph.vertex_rep[color_rep[i]] << " with color " << i << "\n";);
            mgraph.separate_vertices(mgraph.vertex_rep[v], mgraph.vertex_rep[color_rep[i]]);
            touched_vertices.push_back(mgraph.vertex_rep[v]);
            touched_vertices.push_back(mgraph.vertex_rep[color_rep[i]]);
        }
    }
    else {
        assert(lit > 0);
        assert(vertex_color[v] == NoColor); //since we take care to enforce at-most-one color constraint
        vertex_color[v] = i;
        num_uncolored--;
        if(not mgraph_is_used) {
            return;
        }
        DEBUG_PRINT(std::cout << "Assigned " << v << " to " << i << " and merged with color rep " << color_rep[i] << "\n";);
        //first time color is assigned, separate vertices which were assigned a different color
        if(color_rep[i] == NoVertex) {
            color_rep[i] = v; //v is first vertex with color i
            //check for vertices where this color is forbidden and add an edge
            for (int u : INSTANCE->complement_graph_adjacency[v]) {
                if(mgraph.nodeset[u] and val(xvi_to_lit(u,i)) == -1) {
                    mgraph.separate_vertices(mgraph.vertex_rep[u], mgraph.vertex_rep[v]);
                    touched_vertices.push_back(mgraph.vertex_rep[u]);
                    touched_vertices.push_back(mgraph.vertex_rep[v]);
                }
            }
        }
        //same color assigned, merge to vertices of color class
        else if(not mgraph.has_edge(v,mgraph.vertex_rep[color_rep[i]])){ //edge might be there but a future propagation will produce the conflict
            mgraph.contract_vertices(v, mgraph.vertex_rep[color_rep[i]]); //merge v into vertices already colored with i
            touched_vertices.push_back(mgraph.vertex_rep[v]);
        }
        //separate vertex from other color classes
        for (int c = 0; c < num_colors; c++) {
            if(c != i and color_rep[c] != NoVertex and not mgraph.has_edge(mgraph.vertex_rep[v], mgraph.vertex_rep[color_rep[c]])) {
                DEBUG_PRINT(std::cout << "Added edge from " << v << " with color " << i << " to "
                            << mgraph.vertex_rep[color_rep[c]] << " with color " << c << "\n";);
                mgraph.separate_vertices(mgraph.vertex_rep[v], mgraph.vertex_rep[color_rep[c]]);
                touched_vertices.push_back(mgraph.vertex_rep[v]);
                touched_vertices.push_back(mgraph.vertex_rep[color_rep[c]]);
            }
        }
    }
}

void CadicalAssignmentPropagator::unassign_color(const int lit) {
    auto [v,i] = lit_to_xvi(lit);
    if(lit < 0) {
        //color is allowed again
        forbidden_colors[v]--;
    }
    else{
        //vertex was assigned color by this literal, reset color
        assert(lit > 0);
        assert(vertex_color[v] == i);
        vertex_color[v] = NoColor;
        num_uncolored++;
        if(color_rep[i] == v) {
            color_rep[i] = NoVertex; // v was rep of color i and was uncolored, color i has no rep anymore
        }
    }
}


void CadicalAssignmentPropagator::propagation_helper(const int lit, const int plit) {
    assert(lit > 0);
    if(val(plit) == -1) {
        //lit is already been set by propagations, do not propagate again
        return;
    }
    if(val(plit) == 0) {
        propagations.push_back(-plit);
        //reason: not the propagated lit or not the assigned lit, they are adjacent, or same vertex but different color
        clauses.push_back({-plit, -lit});
    }
    else {
        assert(val(plit) == 1);
        //case that plit is already assigned to 1 leads to conflict and thus backtrack, propagate it first
        propagations.push_front(-plit);
        clauses.push_front({-plit, -lit});
    }
    DEBUG_PRINT(std::cout << "Added literal " << -plit << " x_" << std::get<0>(lit_to_xvi(plit)) << "," << std::get<1>(lit_to_xvi(plit))
        << " to propagation queue with reason " << -plit << " " << -lit << "\n";);
}

void CadicalAssignmentPropagator::propagate_assignment(const int lit) {
    auto [v,i] = lit_to_xvi(lit);
    if(lit > 0) {
        DEBUG_PRINT(std::cout << "propagate literal " << lit << " x_" << v << "," << i << " = 1\n";);
        //propagate color of vertex, i.e. set all x_u,i to false for all neighbours u
        for(int u : neighbors[v]) {
            int plit = xvi_to_lit(u, i);
            propagation_helper(lit, plit);
        }
        // also set all other x_v,c to false for all other colors c != i, this enforces the amo colors constraint
        for(int color = 0; color < num_colors; ++color) {
            if(color == i) {
                continue;
            }
            int plit = xvi_to_lit(v, color);
            propagation_helper(lit, plit);
        }
    }
    else {
        //propagate that vertex is not assigned that color, i.e. do nothing
        DEBUG_PRINT(std::cout << "propagate literal " << lit << " x_" << v << "," << i << " = 0\n";);
    }
}

int CadicalAssignmentPropagator::next_decision() const {
    if(num_uncolored == 0) {
        //all vertices have been covered by a color but there are still undecided variables
        // the rest of the variables could then be treated as don't cares, we just let cadical decide the rest until we get a full model
        // we could also propgate any unset variables to 0 instead of deciding on them, not sure if that makes a difference
        DEBUG_PRINT(std::cout << "There were no more uncolored vertices\n";);
        return 0;
    }
    switch (options.assignment_propagator_decision_strategy) {
        case Options::CadicalAssignment:
            return 0; //let cadical decide
        case Options::Dsatur:
        case Options::DsaturPass:
            return dsatur_decision();
        case Options::LargestUncoloredDegree:
            return largest_degree_decision();
        default:
            throw std::runtime_error("Invalid decision strategy or Zykov propagator.");
    }
}

int CadicalAssignmentPropagator::num_uncolored_neigbors(const int v) const {
    return static_cast<int>(std::count_if(neighbors[v].begin(), neighbors[v].end(),
                                    [&](const int u) {return vertex_color[u] == NoColor;}));
}

int CadicalAssignmentPropagator::dsatur_decision() const {
    //dsatur strategy, find vertex with largest degree saturation and assign it the lowest possible color
    //takes forbidden colors into account, i.e. not just adjacent colors
    int max_saturation = -1;
    std::vector<int> candidates;
    for(int vertex = 0; vertex < num_vertices; vertex++){
        //find uncolored vertex with higher saturation or same saturation and larger degree, update
        if(vertex_color[vertex] == NoColor) {
            if(forbidden_colors[vertex] > max_saturation) {
                max_saturation = forbidden_colors[vertex];
                candidates = {vertex};
            }
            else if(forbidden_colors[vertex] == max_saturation) {
                candidates.push_back(vertex);
            }
        }
    }
    //should only hapen if all vertices have been colored, we check that beforehand
    assert(max_saturation != -1 and not candidates.empty());

    int chosen_vertex;
    if(options.assignment_propagator_decision_strategy == Options::DsaturPass) {
        //use pass tie-breaking strategy
        chosen_vertex = pass_vertex_decision(candidates);
    }
    else {
        //break ties by maximum degree in uncolored subgrapgh
        chosen_vertex = *std::max_element(candidates.begin(), candidates.end(),
            [this](const int v1, const int v2) {return num_uncolored_neigbors(v1) < num_uncolored_neigbors(v2);});
    }

    //find lowest possible color for max saturated vertex, if any exists (do not create new color)
    int lowest_color = find_lowest_color(chosen_vertex);

    int dlit = xvi_to_lit(chosen_vertex, lowest_color);
    assert(vertex_color[chosen_vertex] == NoColor);
    assert(val(dlit) == 0);
    DEBUG_PRINT(std::cout << "Decided on " << dlit << "(x_" << chosen_vertex << "," << lowest_color
                          << ", saturation " << max_saturation << ")\n";);
    return dlit; //return positive assignment since it has a stronger effect than simply forbidding the color
}


int CadicalAssignmentPropagator::pass_same(const int v1, const int v2) const {
    //only defined for adjacent vertices
    assert(original_graph.has_edge(v1, v2));
    //this function counts the number of common available colors for v1 and v2
    int same_count = 0;
    for(int i = 0 ;i < num_colors; i++){
        //check if color available for both vertices
        if(val(xvi_to_lit(v1, i)) == 0 and val(xvi_to_lit(v2, i)) == 0){
            same_count++;
        }
    }
    return same_count;
}

int CadicalAssignmentPropagator::pass_vertex_decision(const std::vector<int> &candidates) const {
    int pass_vertex = candidates.front();
    int pass_value = 0;
    std::vector<int> same_values(candidates.size(), 0); //vector to store overall pass values in
    for(int i = 0; i < candidates.size(); i++){
        int v = candidates[i];
        for(int j = i + 1; j < candidates.size(); j++){
            int other_v = candidates[j];
            //compute pass value if vertices are different, value is 0 if verices are adajcent
            if(v != other_v and original_graph.has_edge(v, other_v)){
                int pass_same_value = pass_same(v,other_v); //pass_same is symmetric, only compute once
                same_values[i] += pass_same_value;
                same_values[j] += pass_same_value;
            }
        }
        if(same_values[i] > pass_value){
            pass_vertex = v;
            pass_value = same_values[i];
        }
    }
    DEBUG_PRINT(std::cout << "pass vertex " << pass_vertex << " value " << pass_value << "\n";);
    return pass_vertex;
}

int CadicalAssignmentPropagator::largest_degree_decision() const {
    //largest degree strategy, find uncolored vertex with largest degree in uncolored subgraph
    //break ties with general vertex degree
    int largest_degree_vertex = NoVertex;
    int max_uncolored = -1;
    for(int vertex = 0; vertex < num_vertices; vertex++){
        //find uncolored vertex with more uncolored neighbors or same saturation and larger degree, update
        if(vertex_color[vertex] == NoColor) {
            int num_uncolored = num_uncolored_neigbors(vertex);
            if(num_uncolored > max_uncolored
                or (num_uncolored == max_uncolored and neighbors[vertex].size() > neighbors[largest_degree_vertex].size())){
                max_uncolored = num_uncolored;
                largest_degree_vertex = vertex;
            }
        }
    }
    //should only hapen if all vertices have been colored, we check that beforehand
    assert(largest_degree_vertex != NoVertex);

    //find lowest possible color for max saturated vertex, if any exists (do not create new color)
    int lowest_color = find_lowest_color(largest_degree_vertex);


    int dlit = xvi_to_lit(largest_degree_vertex, lowest_color);
    assert(vertex_color[largest_degree_vertex] == NoColor);
    assert(val(dlit) == 0);
    DEBUG_PRINT(std::cout << "Decided on " << dlit << "(x_" << largest_degree_vertex << "," << lowest_color
                          << ", saturation " << largest_degree_vertex << ")\n";);
    return dlit;
}

int CadicalAssignmentPropagator::find_lowest_color(int selected_vertex) const {
    int lowest_color = num_colors;
    for (int color = 0; color < num_colors; color++) {
        assert(val(xvi_to_lit(selected_vertex, color)) != 1);
        if(val(xvi_to_lit(selected_vertex, color)) == 0) {
            //color still available
            lowest_color = color;
            break;
        }
    }
    //should not happen because cadical should produce a conflict before this already
    assert(lowest_color != num_colors);
    return lowest_color;
}


bool CadicalAssignmentPropagator::check_consisteny() const {
    //check that assignment and vertex color match
    for (int v = 0; v < num_vertices; v++) {
        for (int i = 0; i < num_colors; ++i) {
            int lit = xvi_to_lit(v, i);
            if(val(lit) <= 0) {
                assert(vertex_color[v] != i);
            }
            if(val(lit) == 1) {
                assert(vertex_color[v] == i); //check equality since we enforce amo color constraint
            }
        }
    }

    //check that num_uncolored count is correct
    assert(std::count(vertex_color.begin(), vertex_color.end(), NoColor) == num_uncolored);

    if(not mgraph_is_used) {
        return true;
    }
    //verify that current mgraph fully corresponds to graph induced by partial coloring
    for (int c = 0; c < num_colors; c++) {
        for (int v = 0; v < num_vertices; v++) {
            if(val(xvi_to_lit(v,c)) == 1) {
                // if(std::find(propagations.begin(), propagations.end(), -xvi_to_lit(v, vertex_color[v])) != propagations.end()
                //     or std::find(propagations.begin(), propagations.end(), -xvi_to_lit(v, c)) != propagations.end()) {
                //     continue;
                // }
                assert(vertex_color[v] == c);
                assert(color_rep[c] != NoVertex);
                //assert that all vertices of same color are contracted in reduced graph
                assert(mgraph.has_edge(v, mgraph.vertex_rep[color_rep[c]]) or mgraph.is_contracted(v, mgraph.vertex_rep[color_rep[c]]));
            }
            else if(val(xvi_to_lit(v,c)) == -1) {
                //assert that if a color was forbidden and its rep exists, there is supposed to be an edge between the two
                assert(color_rep[c] == NoVertex or mgraph.has_edge(v, mgraph.vertex_rep[color_rep[c]]));
            }
        }
    }
    //test that all colored vertices form a clique
    Bitset testclique(num_vertices);
    for (int c = 0; c < num_colors; c++) {
        if(color_rep[c] != NoVertex) {
            testclique.set(mgraph.vertex_rep[color_rep[c]]);
        }
    }
    assert(mgraph.is_clique(testclique));

    return true;
}

void CadicalAssignmentPropagator::compute_cliques() {
    PROP_TIMING(stats.start_phase(Statistics::PropagatorComputeCliques););
    maximal_cliques.clear();
    max_clique_size = -1;
    num_assignments_last_clique_computation = stats.prop_num_assignments;
    //run any clique algorithm to store largest cliques in maximal_cliques
    max_clique_size = mgraph.greedy_cliques(maximal_cliques, options.prop_clique_limit);
    assert(not maximal_cliques.empty());
    stats.prop_num_clique_computations++;
    stats.prop_num_maximal_cliques_computed += static_cast<int>(maximal_cliques.size());
    if(max_clique_size == num_colors) {
        stats.prop_num_tight_cliques_computed += static_cast<int>(maximal_cliques.size());
    }
    DEBUG_PRINT(std::cout << "Found " << maximal_cliques.size() << " cliques of size " << max_clique_size << "\n";);
    PROP_TIMING(stats.end_phase(Statistics::PropagatorComputeCliques););
}

bool CadicalAssignmentPropagator::need_to_recompute_cliques() const {
    //cliques might not yet exist or be valid anymore if new assignments were made on the same level
    return maximal_cliques.empty() or (num_assignments_last_clique_computation != stats.prop_num_assignments);
}

bool CadicalAssignmentPropagator::compute_clique_clauses() const {
    assert(propagations.empty());
    return options.use_clique_explanation_clauses and first_call_after_backtrack;
}

bool CadicalAssignmentPropagator::compute_mycielsky_clauses() const {
    assert(propagations.empty());
    return options.use_mycielsky_explanation_clauses and first_call_after_backtrack and external_clauses.empty();
}

std::vector<int> CadicalAssignmentPropagator::subgraph_explanation_clause(const Bitset &nodeset) const {
    std::vector<int> clause;
    //for each v in clique and u in bag r(v), add literal -x_u,i for color i of bag r(v)
    for (int v = nodeset.find_first(); v != Bitset::npos; v = nodeset.find_next(v)) {
        assert(v == mgraph.vertex_rep[v] and mgraph.nodeset[v]);
        int v_col = vertex_color[v];
        if(v_col == NoColor) {
            //might be that vertex is uncolored but we set a literal to false, forbidding it from a different color
            for (int c = 0; c < num_colors; c++) {
                if(color_rep[c] == NoVertex) { //color has not been assigned yet and thus no edge added
                    continue;
                }
                //edge is there and color variable is set and color rep is in clique
                int c_rep = mgraph.vertex_rep[color_rep[c]];
                if(nodeset[c_rep] and mgraph.has_edge(v, c_rep) and not original_graph.has_edge(v, c_rep)) {
                    //add x_(v,c) to clause
                    assert(val(xvi_to_lit(v,c)) == -1);
                    clause.push_back( xvi_to_lit(v,c));
                    assert(clause.back() != 0);
                }
            }
            continue; //part of clique but not colored yet
        }
        for (int u : mgraph.bag[v]) {
            assert(mgraph.is_contracted(u,v));
            assert(val(xvi_to_lit(u,v_col)) == 1);
            //add -x_(u,v_col) to clause
            clause.push_back( - xvi_to_lit(u,v_col));
            assert(clause.back() != 0);
        }
    }
    return clause;
}

void CadicalAssignmentPropagator::check_for_clique_clauses() {
    if(need_to_recompute_cliques()) {
        compute_cliques();
    }
    assert(not maximal_cliques.empty());
    PROP_TIMING(stats.start_phase(Statistics::PropagatorCliqueClauses););
    //check for cliques of size > num_colors and add clauses that one of the edges of the clique has to be contracted to avoid a new color
    //(because a clique of that size would need more colors than we are currently allowing)
    if(max_clique_size > num_colors) {
        for(const auto & clique : maximal_cliques) {
            //clique is too large and one edge should have been contracted instead, add the reason clause for this
            add_clique_explanation_clause(clique);
            stats.prop_clique_pruning_level[current_level]++;
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorCliqueClauses););
}

void CadicalAssignmentPropagator::add_clique_explanation_clause(const Bitset &clique) {
    DEBUG_PRINT(std::cout << "Trying to explain clique " << clique << " of size " << clique.count() << "\n";);
    assert(mgraph.is_clique(clique));
    external_clauses.push_back(subgraph_explanation_clause(clique));
    stats.prop_num_clique_successes++;
    stats.backtrack_reson = 1;
    DEBUG_PRINT(std::cout << "Gave explanation " << external_clauses.back() << " for clique " << clique << " of size " << clique.count() << "\n";);
}

void CadicalAssignmentPropagator::check_for_mycielsky_clauses() {
    if(need_to_recompute_cliques()) {
        compute_cliques();
    }
    assert(not maximal_cliques.empty());
    PROP_TIMING(stats.start_phase(Statistics::PropagatorMycielskyClauses););
    assert(not maximal_cliques.empty());
    assert(max_clique_size <= num_colors); //otherwise clique pruning would have been done instead
    int gap = num_colors - max_clique_size;
    if(gap < options.mycielsky_threshold) {
        for(const auto & clique : maximal_cliques) {
            if(gap >= stats.mycielsky_calls.size()) {
                stats.mycielsky_calls.resize(gap + 1, 0);
            }
            stats.mycielsky_calls[gap]++;
            //start algorithm with clique subgraph; given as bitset and build graph from that
            MGraph::SubGraph subgraph(clique);
            //try to extend subgraph to next generalised mycielsky, need at least gap + 1 succesful iterations
            int iterations = mgraph.mycielsky_extension_clique(subgraph, gap + 1);
            if(iterations > gap) { //only explain subgraph if the proven bound was good enough to cross threshold
                assert(iterations == gap + 1);
                add_mycielsky_explanation_clause(subgraph);
                if(gap >= stats.mycielsky_sucesses.size()) {
                    stats.mycielsky_sucesses.resize(gap + 1, 0);
                }
            stats.mycielsky_sucesses[gap]++;
            stats.prop_myc_pruning_level[current_level]++;
            stats.backtrack_reson = 2;
            }
        }
    }
    PROP_TIMING(stats.end_phase(Statistics::PropagatorMycielskyClauses););
}

void CadicalAssignmentPropagator::add_mycielsky_explanation_clause(const MGraph::SubGraph &subgraph) {
    DEBUG_PRINT(std::cout << "Trying to explain mycielsky subgraph " << subgraph.nodes << " of size " << subgraph.num_vertices << "\n";);
    external_clauses.push_back(subgraph_explanation_clause(subgraph.nodeset));
    DEBUG_PRINT(std::cout << "Gave explanation " << external_clauses.back() << " for mycielsky subgraph " << subgraph.nodes << " of size " << subgraph.num_vertices << "\n";);
}


int CadicalAssignmentPropagator::find_dominated_vertex_decisions() {
    assert(not INSTANCE->complement_graph_adjacency.empty());
    //touched vertices contains all vertices who had their neighbourhood changed, either through contraction or adding an edge
    for (int i : touched_vertices){
        if(not mgraph.nodeset[i] or vertex_color[i] == NoColor) {
            //vertex might have been contracted into another, or has not been colored yet
            continue;
        }
        for (int j : INSTANCE->complement_graph_adjacency[i]){ //can only dominate each other if not adjacent
            if(not mgraph.nodeset[j] or mgraph.has_edge(i,j) or vertex_color[j] != NoColor) {
                continue;
            }
            //vertex i cannot become dominated by edge additions, so only need to check whether it dominates j
            if(mgraph.gmatrix[j].is_subset_of(mgraph.gmatrix[i])) {
                //j is dominated by i since N(i) contains all of N(j) or vice versa
                assert(vertex_color[i] != NoColor and vertex_color[j] == NoColor);
                //assign vertex j the same color as vertex i
                int lit = xvi_to_lit(j, vertex_color[i]);
                assert(val(lit) == 0);
                DEBUG_PRINT(
                std::cout << "DomVertex decided on " << lit  << " = x_" << j << "," << vertex_color[i] << " since "<< i << " and " << j << " dominate each other"
                          << " on level " << current_level << "\n";
                );
                stats.prop_num_dominated_vertex_decisions++;
                return lit;
            }
        }
    }
    return 0;
}




//
// void CadicalAssignmentPropagator::find_vertex_fusions() {
//     PROP_TIMING(stats.start_phase(Statistics::PropagatorPositivePruning););
//     assert(max_clique_size == num_colors);
//     //test for vertex fusion, i.e., u in C and non-adjacent v in V s.t. v is adjacent to all of C\u, then we can merge u,v
//     for(auto & clique : maximal_cliques) {
//         bool exit = false;
//         for (int u = clique.find_first(); u != Bitset::npos and not exit; u = clique.find_next(u)) {
//             clique.reset(u); //for easier adjacency check with vertices v, intersection has to be clique minus u
//             for (auto v : INSTANCE->complement_graph_adjacency[u]) {
//                 if(not mgraph.nodeset[v] or mgraph.has_edge(u,v)
//                     or ((vertex_color[u] == NoColor) == (vertex_color[v] == NoColor))) {//exactly one colored
//                     continue;
//                 }
//                 if(clique.is_subset_of(mgraph.gmatrix[v])) {
//                     //in this case, we can merge u and v, i.e. assign them the same color
//                     assert(vertex_color[u] != NoColor or vertex_color[v] != NoColor);
//                     int propagated_lit = vertex_color[u] != NoColor ? xvi_to_lit(v, vertex_color[u]) : xvi_to_lit(u, vertex_color[v]);
//                     propagations.push_back(propagated_lit);
//                     //as reason clause, we give all edges of the clique and edges from v to the clique
//                     clique.set(u); //include u in clique explanation!
//                     std::vector<int> clause = subgraph_explanation_clause(clique);
//                     clique.reset(u);
//                     //now consider edges from v to clique adn where they come from
//                     for (int w = clique.find_first(); w != Bitset::npos; w = clique.find_next(w)) {
//                         assert(mgraph.vertex_rep[w] == w);
//                         assert(mgraph.has_edge(v,w));
//                         int w_col = vertex_color[w];
//                         if(w_col != NoColor) {
//                             //vertex w already explained in clique explanation
//                             continue;
//                         }
//                         if(vertex_color[v] == NoColor) {
//                             //both uncolored, edge must have existed in original graph
//                             assert(original_graph.has_edge(v,w));
//                         }
//                         else {
//                             //vertex v is colored, check that x_w,vcol is false to have produced this edge
//                             assert(val(xvi_to_lit(w,vertex_color[v])) == -1);
//                             clause.push_back( xvi_to_lit(w, vertex_color[v]));
//                             assert(clause.back() != 0);
//                         }
//                     }
//
//                     clause.push_back(propagated_lit);
//                     clauses.push_back(clause);
//                     stats.prop_positive_prunings++;
//                     DEBUG_PRINT(std::cout << "Vertex fusion propagates " << propagated_lit << " merge " << u << ","
//                         << v << " with reason " << clauses.back() << " on level " << current_level
//                         << " clique " << clique << " and N(v) " << mgraph.gmatrix[v] << "\n";);
//                     exit = true;//only want to find one conflict per clique, as a balance with speed
//                     break;
//                 }
//             }
//             clique.set(u);
//         }
//     }
//     PROP_TIMING(stats.end_phase(Statistics::PropagatorPositivePruning););
// }
//
// void CadicalAssignmentPropagator::find_edge_additions() {
//     PROP_TIMING(stats.start_phase(Statistics::PropagatorNegativePruning););
//     assert(max_clique_size == num_colors);
//     //test for possible edge addition, i.e., for u,v in V\C, test that for all w, (u,w) or (v,w) in E
//     //only iterate over vertices whose neighbourhood changed
//     for(int u : mgraph.vertices){
//         for (auto v : INSTANCE->complement_graph_adjacency[u]) {
//             if(u > v or not mgraph.nodeset[v] or mgraph.has_edge(u,v)
//                 or ((vertex_color[u] == NoColor) == (vertex_color[v] == NoColor))) { //exactly one colored
//                 continue;
//             }
//             for(auto & clique : maximal_cliques) {
//                 if(clique[u] or clique[v]) {
//                     continue;
//                 }
//                 //got u,v in V\C, now check clique for condition C subset (N(u) U N(v))
//                 if(clique.is_subset_of (mgraph.gmatrix[u] | mgraph.gmatrix[v])) {
//                     //in this case, we can add an edge between u and v, i.e. assign them different colors0
//                     assert((vertex_color[u] != NoColor or vertex_color[v] != NoColor)
//                             and not (vertex_color[u] != NoColor and vertex_color[v] != NoColor));
//                     int propagated_lit = (vertex_color[u] != NoColor) ? -xvi_to_lit(v, vertex_color[u]) : -xvi_to_lit(u, vertex_color[v]);
//                     assert(propagated_lit != 0);
//                     propagations.push_back(propagated_lit);
//                     //as reason clause, we give all edges of the clique and edges from x and y to the clique
//                     //add colored vertex to clique and explain it automatically
//                     clique.set(u);
//                     clique.set(v);
//                     std::vector<int> clause = subgraph_explanation_clause(clique);
//                     clique.set(u);
//                     clique.set(v);
//
//                     clause.push_back(propagated_lit);
//                     clauses.push_back(clause);
//                     stats.prop_negative_prunings++;
//                     DEBUG_PRINT(std::cout << "Edge addition propagates " << propagated_lit << " edge " << u << ","
//                             << v << " with reason " << clauses.back() << " on level " << current_level
//                             << " clique " << clique << " and N(u) " << mgraph.gmatrix[u] << " and N(v) " << mgraph.gmatrix[v] << "\n";);
//                 }
//             }
//         }
//     }
//     PROP_TIMING(stats.end_phase(Statistics::PropagatorNegativePruning););
// }
//
// void CadicalAssignmentPropagator::find_clique_based_pruning() {
//     assert(not INSTANCE->complement_graph_adjacency.empty());
//     if(not first_call_after_backtrack) { //computing clique pruning at every node is too expensive
//         return;
//     }
//    if(need_to_recompute_cliques() and (options.enable_positive_pruning or options.enable_negative_pruning)) {
//     compute_cliques();
//    }
//     //prunings only hold if clique is same size as num_colors we are looking for
//     if(max_clique_size != num_colors) {
//         return;
//     }
//     //check for positive and negative pruning, i.e., vertex fusions and edge additions
//     if(options.enable_positive_pruning) {
//         find_vertex_fusions();
//     }
//     if(options.enable_negative_pruning and propagations.empty()) {
//         find_edge_additions();
//     }
//     //since these are prunings that need to hold (otherwise problem is not k-colorable),
//     //it is okay to propagate all of them at the same time, even if they conflict each other.
//     //with correct handling, this simply leads to a conflict and backtrack
// }
