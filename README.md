A simple C++ genetic algorithm class that's meant to be easily usable and extensible.

Current TODO list:
- [ ] Split the code into multiple files
- [x] Fix the indent
- [ ] Make the genetic algorithm class more customizable, i.e. custom gene types, more general and scalable evolution, ...
- [ ] Recode the canon algorithm entirely (utterly bad code)
- [ ] Reduce copies within the algorithms (optimization)

The example is using *onilib*.
Requires *C++11* to compile (for <random>, range-based loops and probably more I'm missing).
Give _-std=c++11_ as a compiler argument. Recent g++ or clang++ are required in order to do this.
