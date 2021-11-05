#pragma once

template<typename T>
T snoob(T x) {
    T y = x & -x, z = x + y;
    return z | ((x ^ z) >> 2) / y;
}

long long spanning_tree_num(graph& g) {
    if(g.size() == 0) return 0;
    if(g.size() == 1) return 1;
    int x = 0, y = 1;
    matrix<fraction<long long>> laplacian(g.size(), g.size());
    for(int i = 0; i < g.size(); ++i) {
        for(int j : g[i]) {
            laplacian(i, j) -= fraction<long long>(1);
        }
    }
    for(int i = 0; i < g.size(); ++i) {
        fraction<long long> sm(0);
        for(int j = 0; j < g.size(); ++j) {
            if(i == j) {
                continue;
            }
            sm -= laplacian(i, j);
        }
        laplacian(i, i) = sm;
    }
    matrix<fraction<long long>> cof11(g.size() - 1, g.size() - 1);
    for(int i = 1; i < g.size(); ++i) {
        for(int j = 1; j < g.size(); ++j) {
            cof11(i - 1, j - 1) = laplacian(i, j);
        }
    }
    fraction<long long> det11;
    int r11;
    cof11.rref(det11, r11);
    assert(det11.q == 1);
    return det11.p;
}

int spanset_num(graph& g, int ex) {
    int tot = g.size() - 1 + ex, res = 0;
    std::vector<int> e1, e2;
    for(int i = 0; i < g.size(); ++i) {
        for(int j : g[i]) {
            if(i < j) {
                e1.push_back(i);
                e2.push_back(j);
            }
        }
    }
    assert(e1.size() <= 60);
    assert(tot <= e1.size());
    long long msk = 0;
    for(int i = 0; i < tot; ++i) msk |= (1LL << i);
    while(msk < (1LL << e1.size())) {
        unionfind uf(g.size());
        for(int i = 0; i < e1.size(); ++i) {
            if(msk & (1LL << i)) {
                uf.unite(e1[i], e2[i]);
            }
        }
        if(uf.c == 1) res++;
        msk = snoob(msk);
    }
    return res;
}

int spanset_num(ribbon_graph& r, int ex) {
    graph G = r.to_graph();
    return spanset_num(G, ex);
}

int dual_spanset_num(ribbon_graph& r, int ex, int fix = -1) {
    std::vector<int> edg = r.edge_list();
    ribbon_graph s = r.dual();
    int tot = r.v - 1 + ex, res = 0;
    if(ex < 0 || tot > r.edges()) return 0;
    if(fix != -1 && (fix < 0 || fix >= r.edges())) return 0;
    assert(r.edges() <= 60);
    long long msk = 0;
    for(int i = 0; i < tot; ++i) msk |= (1LL << i);
    while(msk < (1LL << r.edges())) {
        if(fix != -1 && !(msk & (1LL << fix))) continue;
        std::vector<int> keep, dualkeep;
        for(int i = 0; i < r.edges(); ++i) {
            if(msk & (1LL << i)) {
                keep.push_back(edg[i]);
            } else {
                dualkeep.push_back(edg[i]);
            }
        }
        ribbon_graph h = r.span_subgraph(keep);
        ribbon_graph hd = s.span_subgraph(dualkeep);
        if(h.components() == 1 && hd.components() == 1) res++;
        msk = snoob(msk);
    }
    return res;
}

int quasi_num(ribbon_graph& r, int ex, int fix = -1) {
    int tot = r.v - 1 + ex, res = 0;
    if(ex < 0 || tot > r.edges()) return 0;
    if(fix != -1 && (fix < 0 || fix >= r.edges())) return 0;
    assert(r.edges() <= 60);
    std::vector<int> edg = r.edge_list();
    long long msk = 0;
    for(int i = 0; i < tot; ++i) msk |= (1LL << i);
    while(msk < (1LL << r.edges())) {
        if(fix != -1 && !(msk & (1LL << fix))) continue;
        std::vector<int> keep;
        for(int i = 0; i < r.edges(); ++i) {
            if(msk & (1LL << i)) {
                keep.push_back(edg[i]);
            }
        }
        ribbon_graph sub = r.span_subgraph(keep);
        if(sub.components() == 1 && sub.boundary_components() == 1) res++;
        msk = snoob(msk);
    }
    return res;
}