# IncSatGC - Incremental SAT Solving for the Graph Coloring problem

This repository contains the code for my master thesis project 
looking at different SAT Solving approaches to the graph coloring problem.
These include: 
1. Assignment Encoding
2. Partial Order Encoding
3. Approaches using the Zykov Encoding as in [[1]](#1)
   1. As MaxSAT Encoding (only writing encoding to a .wncf file)
   2. Full SAT Encoding
   3. Incremental SAT Encoding 
4. Propagators using the IPASIR-UP interface from [[2]](#2)
   1. A propagator for the Zykov encoding that propagates the transitivity constraints
      and allows to prune the search tree with lower bounds
   2. A propagator for the assignment encoding that additionally builds a reduced graph on which to compute lower bounds to prune the search tree


## Dependencies and Installation

We use 
[Open-WBO](https://github.com/sat-group/open-wbo) 
with [Glucose 4.2](https://github.com/audemard/glucose)
to build the incremental cardinality constraints with the Totalizer Encoding.
Additionally, we implemented support to use recent versions of 
[CaDiCal](https://github.com/arminbiere/cadical) and 
[Cryptominisat 5](https://github.com/msoos/cryptominisat)
as solvers too, though the latter has some quirks[[3]](#3).
For Cadical, the “ipasirup-for-rc2” branch (b3b2c85) has to be used, 
since the user interface of that version was used for the propagators.
[Boost](https://www.boost.org/) and ZLIB are also required.
Further, a binary of [CliSAT](https://github.com/psanse/CliSAT) is used to compute an initial clique.

Then to create and use the binary ``IncSatGC`` it should be as simple as
```
mkdir build && cd build
cmake .. -DBOOST_ROOT=/path/to/boost -DOPENWBO_ROOT_DIR=/path/to/open-wbo -DCADICAL_ROOT_DIR=/path/to/cadical -DCLISAT_BINARY_PATH=/path/to/clisat/binary  
make
```

## Usage

```
./IncSatGC inputfile [options]  
```
The input is a graph in dimacs format 
and available options can be obtained with ``-h, --help``.
To Switch between the encodings given above, use ``-e, --encoding``:
0. Assignment Encoding
1. Full MaxSAT Encoding, only writes wcnf file
2. Full Zykov encoding
3. Incremental Zykov Encoding
4. Partial Order Encoding
5. Zykov Propagator
6. Assignment Propagator

The best options are not set by default, always recommended are ``-rom`` 
and additionally ``--remove-cj`` for Zykov based encodings. 
The bottom-up strategy (set by ``-s 1``) seems to work better for Zykov based encodings are the propagators.
For the propagators ``--clique-explanations`` and ``--mycielsky-explanations`` must be set to take advantage of 
their pruning capabilities. 
Other options should be set according to their explanations from ``-h, --help``.



### Reference

<a id="1">[1]</a> 
Glorian, Gael, et al. "An incremental sat-based approach to the graph colouring problem." 
Principles and Practice of Constraint Programming: 25th International Conference, CP 2019

<a id="1">[2]</a>
Fazekas, Katalin, et al. "IPASIR-UP: user propagators for CDCL." 26th International Conference on Theory and Applications of Satisfiability Testing (SAT 2023). Schloss Dagstuhl-Leibniz-Zentrum für Informatik, 2023.

<a id="1">[3]</a> 
Cryptominisat has name clashes with macros in Glucose, 
namely 'var_Undef', 'l_True', 'l_False' and 'l_Undef'.
One has to make a new Cryptominisat installation where one prepends 
those names with 'CM', i.e. 'l_True'->'CMl_True' 
and enable the support of Cryptominisat with ``-DCRYPTOMINISAT=TRUE``
and ``-DCRYPTOMINISAT_ROOT_DIR=/path/to/CMS``.
Otherwise, it is not linked and compiles without the Cryptominisat dependency.
