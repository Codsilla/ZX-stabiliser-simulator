# ZX-stabiliser-simulator
Compute amplitude of a ZX diagram through stabiliser decomposition

Classical simulator of quantum computing with stabilser decompositions through ZX-calculus.
This is a cleaner implementation (and slightly faster on bigger instances) of [this previous implementation](https://github.com/Codsilla/ZX-simulator/) 

This whole project is a re-implementation of the simulator of [quizx](https://github.com/Quantomatic/quizx), but now including graphs cuts (using KaHyPar
hypergraph partionning libray), a bunch of other optimizations such as $\alpha$ estimations and new decompositions.

Note that you'll need to have a working build of KahyPar with bindings for rust which you can find [here](https://github.com/tuomas56/kahypar-rs). You will then have to update the dependencies to include your path to KaHyPar-rs.
