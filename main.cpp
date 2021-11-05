#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cassert>
#include <random>
#include "base.cpp"
#include "bigint.cpp"
#include "fractions.cpp"
#include "matrix.cpp"
#include "calculate_resistance.cpp"
#include "unionfind.cpp"
#include "permutation_util.cpp"
#include "ribbon_graph.cpp"
#include "ribbongen.cpp"
#include "graph_io.cpp"
#include "graph_util.cpp"
#include "spantree_util.cpp"
#include "polynomial.cpp"
#include "poly_stats.cpp"

/*
References:
K7 <-> Heawood
Shrikhande <-> Dyck
T3 <-> K3-3
T4 <-> Cube
Double square is self dual

            R = 1       R = 2       R = 1/2
K7          2/7         1/3         2/9
Heawood     13/21       26/29       13/34
Shrikhande  5/16
Dyck        31/48
T3          2/9         1/4         2/11
K3-3        5/9         10/13       5/14
Cube        7/12
T4          1/4

chi(k_n) = ceil((n - 3)(n - 4) / 12)
chi(k_m,n) = ceil((m - 2)(n - 2) / 4)
*/

int main() {
    graph G(10);
    for(int i = 0; i < 5; ++i)
        for(int j = 5; j < 10; ++j)
            add_edge(G, i, j);
    ribbon_graph r = bogo_embed(G, 3);
    int g = r.genus();
    std::cout << g << std::endl;
    assert(r.valid());
    assert(r.smooth());
    std::cout << "dually spanning G:\n";
    for(int i = 0; i <= 2 * g; ++i) std::cout << i << ": " << dual_spanset_num(r, i) << std::endl;
    std::cout << "quasi-trees G:\n";
    for(int i = 0; i <= 2 * g; ++i) std::cout << i << ": " << quasi_num(r, i) << std::endl;
    std::cout << "tutte poly of G:\n";
    auto poly1 = tutte<int>(G);
    std::cout << poly1 << std::endl;
    std::cout << "kr poly of G:\n";
    auto poly2 = krushkal_ribbon_poly<int>(r);
    std::cout << poly2 << std::endl;
    ribbon_graph s = r.dual();
    std::cout << "kr poly of G^*:\n";
    auto poly3 = krushkal_ribbon_poly<int>(s);
    std::cout << poly3 << std::endl;
}