template <typename T>
polynom<2, T> tutte(graph &g) {
    std::vector<std::pair<int,int>> edg;
    unionfind ufg(g.size());
    for(int i = 0; i < g.size(); ++i) {
        int loopcnt = 0;
        for(int j : g[i]) {
            if(i < j) edg.push_back(std::make_pair(i, j));
            else if(i == j) loopcnt++;
            ufg.unite(i, j);
        }
        for(int j = 0; j < loopcnt / 2; ++j) edg.push_back(std::make_pair(i, i));
    }
    polynom<2, T> res;
    for(long long msk = 0; msk < (1LL << edg.size()); ++msk) {
        unionfind uf(g.size());
        for(int i = 0; i < edg.size(); ++i) {
            if(msk & (1LL << i)) {
                uf.unite(edg[i].first, edg[i].second);
            }
        }
        int e1 = uf.c - ufg.c, e2 = uf.c + __builtin_popcountll(msk) - g.size();
        res[e1][e2] += T(1);
    }
    return res;
}

template <typename T>
polynom<2, T> weighted_tutte(weighted_graph<T> &g) {
    std::vector<std::pair<int,int>> edg;
    std::vector<T> weights;
    unionfind ufg(g.size());
    for(int i = 0; i < g.size(); ++i) {
        std::vector<std::pair<T, int>> loops;
        for(auto p : g[i]) {
            int j = p.second;
            if(i < j) {
                edg.push_back(std::make_pair(i, j));
                weights.push_back(p.first);
            } else if(i == j) {
                loops.push_back(std::make_pair(p.first, i));
            }
            ufg.unite(i, j);
        }
        std::sort(loops.begin(), loops.end());
        for(int j = 0; j < loops.size(); j += 2) {
            edg.push_back(std::make_pair(loops[j].second, loops[j].second));
            weights.push_back(loops[j].first);
        }
    }
    polynom<2, T> res;
    for(long long msk = 0; msk < (1LL << edg.size()); ++msk) {
        unionfind uf(g.size());
        T w = T(1);
        for(int i = 0; i < edg.size(); ++i) {
            if(msk & (1LL << i)) {
                uf.unite(edg[i].first, edg[i].second);
                w /= weights[i];
            }
        }
        int e1 = uf.c - ufg.c, e2 = uf.c + __builtin_popcountll(msk) - g.size();
        res[e1][e2] += w;
    }
    return res;
}

template <typename T>
polynom<3, T> bollobas_riordan_c(ribbon_graph& r) {
    std::vector<int> edg = r.edge_list();;
    int gk = r.components();
    polynom<3, T> res;
    polynom<1, T> t;
    t[1] = T(1); 
    t[0] = T(-1);
    for(long long msk = 0; msk < (1LL << edg.size()); ++msk) {
        std::vector<int> keep;
        for(int i = 0; i < r.edges(); ++i) {
            if(msk & (1LL << i)) {
                keep.push_back(edg[i]);
            }
        }
        ribbon_graph h = r.span_subgraph(keep);
        res[h.genus()][h.nullity()] += t.pow(h.components() - gk);
    }
    return res;
}

template <typename T>
polynom<4, T> krushkal_ribbon_poly(ribbon_graph& r) {
    std::vector<int> edg = r.edge_list();
    ribbon_graph s = r.dual();
    int rk = r.components(), sk = s.components(), g = r.genus();
    polynom<4, T> res;
    for(long long msk = 0; msk < (1LL << edg.size()); ++msk) {
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
        assert(h.nullity() == hd.components() - sk + g + h.genus() - hd.genus());
        res[h.components() - rk][hd.components() - sk][h.genus()][hd.genus()] += T(1);
    }
    return res;
}

template <typename T>
polynom<4, T> krushkal_weighted(weighted_ribbon_graph<T>& w) {
    std::vector<int> edg = w.r.edge_list();
    std::vector<T> weights;
    for(int i : edg) weights.push_back(w.ohm[i]);
    weighted_ribbon_graph<T> v = w.dual();
    int rk = w.r.components(), sk = v.r.components(), g = v.r.genus();
    polynom<4, T> res;
    for(long long msk = 0; msk < (1LL << edg.size()); ++msk) {
        std::vector<int> keep, dualkeep;
        T x = T(1);
        for(int i = 0; i < w.r.edges(); ++i) {
            if(msk & (1LL << i)) {
                keep.push_back(edg[i]);
                x /= weights[i];
            } else {
                dualkeep.push_back(edg[i]);
            }
        }
        ribbon_graph h = w.r.span_subgraph(keep);
        ribbon_graph hd = v.r.span_subgraph(dualkeep);
        assert(h.nullity() == hd.components() - sk + g + h.genus() - hd.genus());
        res[h.components() - rk][hd.components() - sk][h.genus()][hd.genus()] += x;
    }
    return res;
}