#include "Options.h"


Options::Options() {
    encoding = AssignmentEncoding;
    checker = SparseTrianglesChecker;
    strategy = TopDown;
    solver = CaDiCaL;
    verbosity = Normal;
    disable_preprocessing = false;
    reduce_graph = false;
    use_clique_in_ordering = false;
    use_mycielsky_lb = false;
    remove_trivial_cj = false;
    assignment_encoding_amo = false;
    write_cnf_only = false;
    assignment_propagator_decision_strategy = CadicalAssignment;
    zykov_propagator_decision_strategy = CadicalZykov;
    disable_cardinality_constraints = false;
    zykov_coloring_algorithm = None;
    prop_clique_limit = std::numeric_limits<int>::max();
    use_clique_explanation_clauses = false;
    use_mycielsky_explanation_clauses = false;
    mycielsky_threshold = 1;
    use_dominated_vertex_decisions = false;
    enable_positive_pruning = false;
    enable_negative_pruning = false;
    original_paper_configuration = false;
    coloringfilepath = "";
    stats_csvfile = "";
    enable_detailed_backtracking_stats = false;
}

Options::Options(int argc, char **argv) : Options() {
    try {
    std::string sep;
    for (int i = 0; i < argc; ++i) {
        full_cmd += sep;
        full_cmd.append(argv[i]);
        sep = " ";
    }
    int tmp_num_colors;
    po::options_description general_options("Main options");
    //don't add default values to parameters here, as they are given in the basic constructor
    general_options.add_options()
            ("help,h", "Print help message")
            //algorithm configuration
            ("encoding,e", po::value<int>((int*)&encoding),
                    option_description("Set encoding", tostr_encoding).c_str())
            ("checker,c", po::value<int>((int*)&checker),
                    option_description("Set checker algorithm for CEGAR encoding", tostr_checker).c_str())
            ("strategy,s", po::value<int>((int*)&strategy),
                    option_description("Set search strategy", tostr_strategy).c_str())
            ("solver,g", po::value<int>((int*)&solver),
                    option_description("Set solver", tostr_solver).c_str())
            ("verbosity,v", po::value<int>((int*)&verbosity),
                    option_description("Set verbosity level", tostr_verbosity).c_str())
            ("num-colors,k", po::value<int>(&tmp_num_colors),
                    "Solve for specific number of colors")
            ("disable-preprocessing", po::bool_switch(&disable_preprocessing),
                    "Disable preprocessing and initial bounds")
            ("reduce-graph,r", po::bool_switch(&reduce_graph),
                    "Enable removing vertices in preprocessing")
            ("use-clique,o", po::bool_switch(&use_clique_in_ordering),
                    "Enable fixing the first vertices of a clique")
            ("mycielsky-bound,m", po::bool_switch(&use_mycielsky_lb),
                    "Compute initial mycielsky lower bound ")
            ("remove-cj", po::bool_switch(&remove_trivial_cj),
                    "Removes some of the literals from the cardinality constraints if possible")
            ("assignment-amo", po::bool_switch(&assignment_encoding_amo),
                    "Uses at-most-one color clauses for assignment encoding")

            ("write-cnf-only", po::bool_switch(&write_cnf_only),
                    "Option to build but only write cnf to file and not solve (works for -s 2)")

            ("assignment-strategy", po::value<int>((int*)&assignment_propagator_decision_strategy),
            option_description("Options for decision strategy in Assignment propagator", tostr_assignment_strategy).c_str())

            ("zykov-strategy", po::value<int>((int*)&zykov_propagator_decision_strategy),
                    option_description("Options for decision strategy in Zykov propagator", tostr_zykov_strategy).c_str())
            ("disable-cardinality", po::bool_switch(&disable_cardinality_constraints),
                    "Options to disable cardinality constraints and only use cliques to assert coloring size")
            ("coloring-algorithm", po::value<int>((int*)&zykov_coloring_algorithm),
                    option_description("Options for coloring algorithm used in Zykov search tree", tostr_color_algorithm).c_str())

            ("clique-limit", po::value<int>(&prop_clique_limit),
                "Sets a limit of how many clique can be searched for at once in the Zykov propagator")
            ("clique-explanations", po::bool_switch(&use_clique_explanation_clauses),
                    "Enables learning clique explanation clauses in external propagator for zykov encoding")
            ("myc-explanations", po::bool_switch(&use_mycielsky_explanation_clauses),
                    "Enables learning mycielsky graph explanation clauses in external propagator for zykov encoding")
            ("myc-threshold", po::value<int>(&mycielsky_threshold),
                    "Calls mycielsky lower bound algorithm if current_k - lb < threshold (default = 1)")
            ("dominated-decisions", po::bool_switch(&use_dominated_vertex_decisions),
                    "Enables contracting dominated vertices before other variable decisions in external propagator for zykov encoding")
            ("positive-pruning", po::bool_switch(&enable_positive_pruning),
                    "enables clique-based positive pruning in Zykov propagator")
            ("negative-pruning", po::bool_switch(&enable_negative_pruning),
            "enables clique-based negative pruning in Zykov propagator")


            ("original-paper-configuration", po::bool_switch(&original_paper_configuration),
                    "Runs the cegar algorithm with the configuration as close as possible as the original Cegar paper")

#ifdef CRYPTOMINISAT
            ("cms-threads", po::value<int>(&cms_threads),
                    "Change number threads used by Cryptominisat, default is 1")
#endif
            //if given, write or append collected statistics to a csv file
            ("logfile,l", po::value(&stats_csvfile),
                "Path to write collected statistics to")
            //optional filepath of where to write a found coloring
            ("coloring,f", po::value(&coloringfilepath),
                "Path to write the coloring to")
            //optional filepath of where to write a found coloring
            ("detailed-btstats", po::bool_switch(&enable_detailed_backtracking_stats),
                "Enables collection of more detailed backtracking stats")
            ;

    //parse positional option which is the filename
    po::options_description positional_options("Positional options");
    positional_options.add_options()
            ("inputfile",po::value<std::string>(&filepath), "input file");

    po::positional_options_description filename_positional;
        filename_positional.add("inputfile", 1);

    po::options_description cmdline_options;
    cmdline_options.add(general_options).add(positional_options);

    //store the set options in the variable map
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
              .options(cmdline_options)
              .positional(filename_positional)
              .run(), vm);
    po::notify(vm);

    //print help if option is set or no inputfile is given
    if (vm.count("help") or (not vm.count("inputfile"))) {
        std::cout << "USAGE: " << argv[0] << " [options] inputfile \n"
                  << " where input is a graph in dimacs format.\n"  ;
        std::cout << general_options << "\n";
        std::exit(0);
    }
    if(not std::filesystem::exists(filepath)){
        throw std::runtime_error("Input file does not exist.");
    }
    filename = std::filesystem::path(filepath).filename().string();

    //take care that options are set properly and don't conflict

    //cannot have preprocessing disabled but use clique in ordering
    if(disable_preprocessing){
        if(use_clique_in_ordering) {
            throw po::error("Cannot have pre-processing disabled and use a clique in ordering.");
        }
        if(reduce_graph) {
            throw po::error("Cannot have pre-processing disabled and use graph reduction.");
        }
    }

    //if option for specific number of colors is set, also check that a number of colors is given
    if (vm.count("strategy") && strategy == Options::SingleK){
        if (not vm.count("num-colors")){
            throw po::error("If strategy is to check only a specific number of colors, a number has to be supplied.");
        }
        specific_num_colors = tmp_num_colors;

    }
    else if(vm.count("num-colors")){
        if (verbosity >= Normal){
            std::cout << "Warning: number of color is given in top-down/bottom-up, "
                         "this means that number will be used as upper/lower bound. "
                         "This should only be used for testing!\n";
        }
        specific_num_colors = tmp_num_colors;
    }
    //cannot combine full maxsat encoding with any search strategy, and solver does not matter
    if(vm.count("encoding") && encoding == Options::FullMaxSAT){
        if (vm.count("strategy")){
            throw po::error("If Full MaxSAT encoding is used, no search strategy can be set.");
        }
        if (vm.count("solver")){
            throw po::error("If Full MaxSAT encoding is used, no solver is used and only the wncf is written.");
        }
        if (vm.count("checker")){
            throw po::error("If Full MaxSAT encoding is used, no solver is used and only the wncf is written.");
        }
    }

    if(vm["original-paper-configuration"].as<bool>()) {
        if(vm.count("encoding") or vm.count("strategy") or vm.count("checker") or vm.count("solver)")
            or vm.count("remove-cj")) {
            std::cout << "Warning: with original paper configuration, other settings are overriden.\n";
        }
        assert(original_paper_configuration);
        disable_preprocessing = true;
        reduce_graph = false;
        use_clique_in_ordering = false;
        use_mycielsky_lb = false;
        remove_trivial_cj = false;
        solver = Glucose;
        encoding = CEGAR;
        checker = PaperChecker;
        strategy = BottomUp;
    }

    //check that options are within range
    if( not(
            (Encoding::AssignmentEncoding <= encoding and encoding <= Encoding::AssignmentPropagator) and
            (CheckAlgorithm::NaiveChecker <= checker and checker <= CheckAlgorithm::PaperChecker) and
            (SearchStrategy::TopDown <= strategy and strategy <= SearchStrategy::SingleK) and
            (Solver::Glucose <= solver and solver <= Solver::Cryptominisat) and
            (Verbosity::Quiet <= verbosity and verbosity <= Verbosity::Debug) and
            (AssignmentPropagatorDecisionStrategy::CadicalAssignment <= assignment_propagator_decision_strategy and
            assignment_propagator_decision_strategy <= AssignmentPropagatorDecisionStrategy::LargestUncoloredDegree) and
            (ZykovPropagatorDecisionStrategy::CadicalZykov <= zykov_propagator_decision_strategy and
                zykov_propagator_decision_strategy <= ZykovPropagatorDecisionStrategy::BagSize) and
            (ZykovPropagatorColoringAlgorithm::None <= zykov_coloring_algorithm and
                zykov_coloring_algorithm <= ZykovPropagatorColoringAlgorithm::IteratedSEQ)
            )){
        throw po::error("Option was out of range.");
    }

    //if using either propagator, check that solver is cadical
    if(encoding == AssignmentPropagator or encoding == ZykovPropagator) {
        if (solver != CaDiCaL) {
            throw std::runtime_error("Propagators are only supported by Cadical");
        }
        if(not use_clique_explanation_clauses) {
            if(use_mycielsky_explanation_clauses) {
                throw std::runtime_error("Mycielsky explanations only work with clique explanations");
            }
            if(enable_positive_pruning or enable_negative_pruning) {
                throw std::runtime_error("Clique based pruning only work with clique explanations");
            }

        }
    }

    if(encoding == ZykovPropagator and disable_cardinality_constraints and not use_clique_explanation_clauses) {
        throw po::error("Cannot disable cardinality encodign without using clique explanations");
    }


    if(encoding == AssignmentPropagator){
        if(use_clique_explanation_clauses or use_mycielsky_explanation_clauses) {
            throw po::error("Clique and Mycielsky pruning are not supported by Assignment propagator");
        }
        if(enable_positive_pruning or enable_negative_pruning) {
            throw po::error("Positive and negative pruning are not supported by Assignment propagator");
        }
    }

#ifdef CRYPTOMINISAT
    if(cms_threads < 0 or mycielsky_threshold < 0 or prop_clique_limit < 0){
        throw po::error("One of the given options can't be negative");
    }
#endif

    if(zykov_coloring_algorithm != None or enable_detailed_backtracking_stats) {
        std::cout << "Warning: options of coloring heuristic during propagator and detailed backtracking statistics are only for testing.\n";
        if(enable_detailed_backtracking_stats and stats_csvfile.empty()) {
            throw po::error("Only enable detailed backtracking stats if also writing stats to csv file.");
        }
    }

    } catch (const po::error &ex) {
        std::cerr << ex.what() << '\n';
        std::exit(1);
    }
}



void Options::print() const {
    std::cout << "c Options : Instance                = " << filename << "\n";
    std::cout << "c Options : Encoding                = " << tostr_encoding[encoding] << "\n";
    std::cout << "c Options : Checker                 = " << (encoding == Encoding::CEGAR ? tostr_checker[checker] : "none" ) << "\n";
    std::cout << "c Options : Strategy                = " << tostr_strategy[strategy] << "\n";
    std::cout << "c Options : Solver                  = " << tostr_solver[solver] << "\n";
    std::cout << "c Options : Verbosity               = " << tostr_verbosity[verbosity] << "\n";
    if(strategy == SearchStrategy::SingleK and specific_num_colors.has_value()){
    std::cout << "c Options : Solving for k           = " << specific_num_colors.value() << "\n";
    }
    std::cout << "c Options : Preprocessing           = " << (not disable_preprocessing ? "True" : "False") << "\n";
    std::cout << "c Options : Graph reduction         = " << (reduce_graph ? "True" : "False") << "\n";
    std::cout << "c Options : Use clique in ordering  = " << (use_clique_in_ordering ? "True" : "False") << "\n";
    std::cout << "c Options : Use mycielsky bound     = " << (use_mycielsky_lb ? "True" : "False") << "\n";
    std::cout << "c Options : Remove trivial cj       = " << (remove_trivial_cj ? "True" : "False") << "\n";
    if(encoding == AssignmentEncoding) {
    std::cout << "c Options : Assignment at-most-one  = " << (assignment_encoding_amo ? "True" : "False") << "\n";
    }
    if(encoding == AssignmentPropagator) {
    std::cout << "c Options : Assignment decision strategy = " << tostr_assignment_strategy[assignment_propagator_decision_strategy] << "\n";}
    else if(encoding == ZykovPropagator) {
    std::cout << "c Options : Zykov decision strategy = " << tostr_zykov_strategy[zykov_propagator_decision_strategy] << "\n";
    }
    if(encoding == AssignmentPropagator or encoding == ZykovPropagator) {
    std::cout << "c Options : Use clique bounding     = " << (use_clique_explanation_clauses ? "True" : "False") << "\n";
    std::cout << "c Options : Use mycielsky bounding  = " << (use_mycielsky_explanation_clauses ? "True" : "False") << "\n";
    std::cout << "c Options : Use dominated decision  = " << (use_dominated_vertex_decisions ? "True" : "False") << "\n";
    std::cout << "c Options : Use positive pruning    = " << (enable_positive_pruning ? "True" : "False") << "\n";
    std::cout << "c Options : Use negative pruning    = " << (enable_negative_pruning ? "True" : "False") << "\n";
    }
    std::cout << "c #################################\n";
}

void Options::print_header() const {
    std::cout << "c Running : " << full_cmd << "\n";
#ifdef __unix__
    char hostname[1024];
    gethostname(hostname, 1024);
    std::cout << "c Machine : " << hostname << ", pid " << getpid() << "\n";
#endif
    time_t currentTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::cout << "c Date    : " << std::put_time(localtime(&currentTime), "%F %T") << "\n";
}



std::string Options::enum_names_to_string(const std::vector<std::string>& enum_strings) {
    std::string s = "[";
    std::string sep;
    for(const std::string& str : enum_strings){
        s += sep + str;
        sep = ", ";
    }
    return s += "]";
}

std::string Options::option_description(std::string desc, const std::vector<std::string> &enum_strings) {
    if(not enum_strings.empty()) {
        desc += " " + enum_names_to_string(enum_strings);
    }
    return desc;
}









