#pragma once

// O(n^3) reiknirit sem reiknar nákvæmt viðnám
template<typename T>
T calculate_resistance_exact(weighted_graph<T> &g, int x, int y) {
    // Til einföldunar látum við y vera stærra
    if(x > y) {
        std::swap(x, y);
    }
    // Reiknum út laplacian fylki netsins
    matrix<T> laplacian(g.size(), g.size());
    for(int i = 0; i < g.size(); ++i) {
        for(auto p : g[i]) {
            laplacian(i, p.second) -= T(1) / p.first;
        }
    }
    for(int i = 0; i < g.size(); ++i) {
        T sm = T(0);
        for(int j = 0; j < g.size(); ++j) {
            if(i == j) {
                continue;
            }
            sm -= laplacian(i, j);
        }
        laplacian(i, i) = sm;
    }
    // Viðnám er þá gefið sem hlutfall þessara ákveða
    matrix<T> cof11(g.size() - 1, g.size() - 1);
    matrix<T> cofxy(g.size() - 2, g.size() - 2);
    for(int i = 1; i < g.size(); ++i) {
        for(int j = 1; j < g.size(); ++j) {
            cof11(i - 1, j - 1) = laplacian(i, j);
        }
    }
    for(int i = 0; i < g.size(); ++i) {
        if(i == x || i == y) {
            continue;
        }
        int ind1 = i < x ? i : (i < y ? i - 1 : i - 2);
        for(int j = 0; j < g.size(); ++j) {
            if(j == x || j == y) {
                continue;
            }
            int ind2 = j < x ? j : (j < y ? j - 1 : j - 2);
            cofxy(ind1, ind2) = laplacian(i, j);
        }
    }
    // Viðnám reiknað með matrix tree theorem
    T detxy, det11;
    int rxy, r11;
    cof11.rref(det11, r11);
    cofxy.rref(detxy, rxy);
    return detxy / det11;
}

// Reiknar viðnám milli tveggja punkta
// exact breyta segir til um hvort ítrun eða nákvæm aðferð er notuð
template<typename T>
T calculate_resistance(weighted_graph<T> &g, int x, int y, bool exact=false) {
    if(exact) {
        return calculate_resistance_exact(g, x, y);
    }
    // Spenna er þýtt fall, reiknum það með ítrunaraðferð
    std::vector<T> volt(g.size(), T(1) / T(2)), csm(g.size(), T(0));
    for(int i = 0; i < g.size(); ++i) {
        for(auto p : g[i]) {
            csm[i] +=  T(1) / p.first;
        }
    }
    volt[x] = T(1);
    volt[y] = T(0);
    const T ITER_EPS = T(1) / T(1000000000);
    T effr = T(0);
    while(true) {
        for(int i = 0; i < g.size(); ++i) {
            if(i == x || i == y) {
                continue;
            }
            T sm = T(0);
            for(auto p : g[i]) {
                sm += volt[p.second] / p.first;
            }
            sm /= csm[i];
            volt[i] = sm;
        }
        // Straumur um legg a - b er gefinn með (volt[a] - volt[b]) / resistance[a, b]
        // Látum I vera straum út úr x (sem hefur spennu 1)
        // Þá gefur lögmál Ohms að R = V/I = 1/I því spennan er 1
        T current = T(0);
        for(auto p : g[x]) {
            current += (volt[x] - volt[p.second]) / p.first;
        }
        T oldr = effr;
        effr = 1.0 / current;
        if(abs(effr - oldr) < ITER_EPS) {
            break;
        }
    }
    return effr;
}

// Reiknar viðnám milli tveggja hlutmengja
// Gerir það með því að herpa A og B í hnúta
// Kallar svo bara á fall að ofan
template<typename T>
T calculate_resistance(weighted_graph<T> &g, std::vector<int> &A, std::vector<int> &B, bool exact=false) {
    std::vector<int> new_inds(g.size(), -1);
    for(int i : A) {
        new_inds[i] = 0;
    }
    for(int i : B) {
        new_inds[i] = 1;
    }
    int cur_ind = 2;
    for(int i = 0; i < g.size(); ++i) {
        if(new_inds[i] != -1) {
            continue;
        }
        new_inds[i] = cur_ind;
        cur_ind++;
    }
    weighted_graph h(cur_ind);
    std::map<std::pair<int,int>,T> resistances;
    for(int i = 0; i < g.size(); ++i) {
        for(int j = 0; j < g[i].size(); ++j) {
            int a = new_inds[i], b = new_inds[g[i][j].second];
            if(!resistances.count(std::make_pair(a, b))) {
                resistances[std::make_pair(a, b)] = T(0);
            }
            resistances[std::make_pair(a, b)] += T(1) / g[i][j].first;
        }
    }
    for(auto p : resistances) {
        int a = p.first.first, b = p.first.second;
        h[a].push_back(std::make_pair(T(1) / p.second, b));
        h[b].push_back(std::make_pair(T(1) / p.second, a));
    }
    return calculate_resistance(h, 0, 1, exact);
}
