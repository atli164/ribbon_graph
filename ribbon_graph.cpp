#pragma once 

struct ribbon_graph {
    int v, h;
    std::vector<int> lam, tau, rho;
    ribbon_graph() : v(0), h(0), lam(0), tau(0), rho(0) { }
    ribbon_graph(int _v, int _h, std::vector<int> _l, std::vector<int> _t, std::vector<int> _r) 
        : v(_v), h(_h), lam(_l), tau(_t), rho(_r) { }
    ribbon_graph(const ribbon_graph& r) 
        : v(r.v), h(r.h), lam(r.lam), tau(r.tau), rho(r.rho) { }
    bool valid() {
        // All are functions from H
        if(lam.size() != h) return false;
        if(tau.size() != h) return false;
        if(rho.size() != h) return false;
        // Check if rho is permutation
        if(!is_permutation(rho)) return false;
        // Check lam rho = lam
        if(lam != composeperms(lam, rho)) return false;
        // Check tau fixed-point free involution 
        for(int i = 0; i < h; ++i) if(tau[tau[i]] != i) return false;
        return true;
    }
    bool smooth() {
        std::vector<std::vector<int>> lam_inv(v);
        for(int i = 0; i < h; ++i) {
            lam_inv[lam[i]].push_back(i);
        }
        for(int i = 0; i < v; ++i) {
            std::map<int,int> newind;
            for(int j = 0; j < lam_inv[i].size(); ++j) {
                newind[lam_inv[i][j]] = j;
            }
            std::vector<int> rhopart(lam_inv[i].size());
            for(int j = 0; j < lam_inv[i].size(); ++j) {
                if(!newind.count(rho[lam_inv[i][j]])) return false;
                rhopart[j] = newind[rho[lam_inv[i][j]]];
            }
            if(orbit_num(rhopart) > 1) return false;
        }
        return true;
    }
    bool is_loop(int i) {
        return lam[tau[i]] == lam[i];
    }
    bool loopless() {
        for(int i = 0; i < h; ++i) if(is_loop(i)) return false;
        return true;
    }
    int vertices() {
        return v;
    }
    int edges() {
        return h / 2;
    }
    int components() {
        unionfind uf(v);
        for(int i = 0; i < h; ++i) {
            uf.unite(lam[i], lam[tau[i]]);
        }
        return uf.c;
    }
    int boundary_components() {
        unionfind uf(h);
        std::vector<int> deg(v, 0);
        for(int i = 0; i < h; ++i) {
            uf.unite(i, rho[tau[i]]);
            deg[lam[i]]++;
        }
        int isol = 0;
        for(int i = 0; i < v; ++i) if(deg[i] == 0) isol++;
        return isol + uf.c;
    }
    int rank() {
        return vertices() - components();
    }
    int nullity() {
        return edges() - vertices() + components();
    }
    int genus() {
        return (components() + nullity() - boundary_components()) / 2;
    }
    ribbon_graph dual() {
        std::vector<int> rt = composeperms(rho, tau);
        unionfind uf(rt.size());
        for(int i = 0; i < rt.size(); ++i) {
            uf.unite(i, rt[i]);
        }
        int vp = uf.c;
        std::map<int,int> newind;
        for(int i = 0; i < rt.size(); ++i) {
            if(!newind.count(uf.find(i))) {
                int x = newind.size();
                newind[uf.find(i)] = x;
            }
        }
        std::vector<int> newlam(h);
        for(int i = 0; i < h; ++i) {
            newlam[i] = newind[uf.find(i)];
        }
        return ribbon_graph(vp, h, newlam, tau, rt);
    }
    bool is_bridge(int i) {
        return dual().is_loop(i);
    }
    bool bridgeless(int i) {
        return dual().loopless();
    }
    graph to_graph() {
        graph g(v);
        for(int i = 0; i < h; ++i) {
            g[lam[i]].push_back(lam[tau[i]]);
        }
        return g;
    }
    ribbon_graph span_subgraph(std::vector<int> keep) {
        std::vector<int> is_in(h, 0), new_ind(h, -1), old_ind;
        int newh = 0;
        for(int x : keep) {
            if(!is_in[x]) {
                is_in[x] = 1;
                new_ind[x] = newh;
                old_ind.push_back(x);
                newh++;
                is_in[tau[x]] = 1;
                new_ind[tau[x]] = newh;
                old_ind.push_back(tau[x]);
                newh++;
            }
        }
        std::vector<int> newlam(newh), newtau(newh), newrho(newh);
        for(int i = 0; i < newh; ++i) {
            newlam[i] = lam[old_ind[i]];
            newtau[i] = i ^ 1;
            newrho[i] = old_ind[i];
            do {
                newrho[i] = rho[newrho[i]];
            } while(new_ind[newrho[i]] == -1);
            newrho[i] = new_ind[newrho[i]];
        }
        return ribbon_graph(v, newh, newlam, newtau, newrho);
    }
    std::vector<int> edge_list() {
        std::vector<int> edg;
        for(int i = 0; i < h; ++i) {
            if(i < tau[i]) {
                edg.push_back(i);
            }
        }
        return edg;
    }
    /*
    void split_edge(int i) {
        assert(i >= 0 && i < e);
        tau.push_back(e + 1);
        tau.push_back(e);
        int x = lam[i], y = lam[tau[i]];
        lam[tau[i]] = v;
        lam.push_back(v);
        lam.push_back(y);
        rho.push_back(tau[i]);
        rho.push_back(rho[tau[i]]);
        int z = -1;
        for(int j = 0; j < e; ++j) {
            if(rho[j] == tau[i]) {
                z = j;
                break;
            }
        }
        assert(z != -1);
        rho[z] = e + 1;
        rho[tau[i]] = e;
        e += 2;
        v++;
        assert(valid());
        assert(smooth());
    }
    void delete_edge(int x) {
        int i1 = x, i2 = tau[x];
        if(i2 < i1) std::swap(i1, i2);
        for(int i = 0; i < e; ++i) {
            if(rho[i] == i1) rho[i] = rho[i1];
            if(rho[i] == i2) rho[i] = rho[i2];
        }
        tau.erase(tau.begin() + i2);
        tau.erase(tau.begin() + i1);
        lam.erase(lam.begin() + i2);
        lam.erase(lam.begin() + i1);
        rho.erase(rho.begin() + i2);
        rho.erase(rho.begin() + i1);
        e -= 2;
        for(int i = 0; i < e; ++i) {
            if(tau[i] > i2) tau[i] -= 2;
            else if(tau[i] > i1) tau[i]--;
            if(rho[i] > i2) rho[i] -= 2;
            else if(rho[i] > i1) rho[i]--;
        }
        assert(valid());
        assert(smooth());
    }
    void contract_edge(int z) {
        int i1 = z, i2 = tau[z];
        if(i1 > i2) std::swap(i1, i2);
        for(int i = 0; i < e; ++i) {
            if(rho[i] == i1) rho[i] = rho[i2];
            if(rho[i] == i2) rho[i] = rho[i1];
        }
        int x = lam[i1], y = lam[i2];
        if(x > y) std::swap(x, y);
        for(int i = 0; i < e; ++i) {
            if(lam[i] == y) lam[i] = x;
        }
        for(int i = 0; i < e; ++i) {
            if(lam[i] > y) lam[i]--;
        }
        v--;
        tau.erase(tau.begin() + i2);
        tau.erase(tau.begin() + i1);
        lam.erase(lam.begin() + i2);
        lam.erase(lam.begin() + i1);
        rho.erase(rho.begin() + i2);
        rho.erase(rho.begin() + i1);
        e -= 2;
        for(int i = 0; i < e; ++i) {
            if(tau[i] > i2) tau[i] -= 2;
            else if(tau[i] > i1) tau[i]--;
            if(rho[i] > i2) rho[i] -= 2;
            else if(rho[i] > i1) rho[i]--;
        }
        assert(valid());
        assert(smooth());
    }
    int component_num() {
        unionfind uf(v);
        for(int i = 0; i < e; ++i) {
            uf.unite(lam[i], lam[tau[i]]);
        }
        return uf.c;
    }
    std::vector<rotsys> factor_components() {
        unionfind uf(v);
        for(int i = 0; i < e; ++i) {
            uf.unite(lam[i], lam[tau[i]]);
        }
        std::vector<rotsys> res(uf.c);
        std::vector<int> factor_ind(v, -1);
        int cur = 0;
        for(int i = 0; i < v; ++i) {
            if(uf.find(i) == i) {
                factor_ind[i] = cur;
                cur++;
            }
        }
        for(int i = 0; i < v; ++i) {
            factor_ind[i] = factor_ind[uf.find(i)];
        }
        std::vector<int> v2v(v, -1), e2e(e, -1), vcnt(uf.c, 0), ecnt(uf.c, 0);
        for(int i = 0; i < v; ++i) {
            v2v[i] = vcnt[factor_ind[i]];
            vcnt[factor_ind[i]]++;
        }
        for(int i = 0; i < e; ++i) {
            e2e[i] = ecnt[factor_ind[lam[i]]];
            ecnt[factor_ind[lam[i]]]++;
        }
        for(int i = 0; i < uf.c; ++i) {
            res[i].v = vcnt[i];
            res[i].e = ecnt[i];
            res[i].lam = std::vector<int>(e);
            res[i].tau = std::vector<int>(e);
            res[i].rho = std::vector<int>(e);
        }
        for(int i = 0; i < e; ++i) {
            res[factor_ind[lam[i]]].lam[e2e[i]] = v2v[lam[i]];
            res[factor_ind[lam[i]]].tau[e2e[i]] = e2e[tau[i]];
            res[factor_ind[lam[i]]].rho[e2e[i]] = e2e[rho[i]];
        }
        return res;
    }
    void delete_vertex(int x) {
        while(true) {
            int to_del = -1;
            for(int i = 0; i < e; ++i) {
                if(lam[i] == x) to_del = i;
            }
            if(to_del == - 1) break;
            delete_edge(to_del);
        }
        for(int i = 0; i < e; ++i) {
            if(lam[i] > x) lam[i]--;
        }
        v--;
    }
    void split_vertex(int i, int j) {

    }
    */
};

template<typename T>
struct weighted_ribbon_graph {
    ribbon_graph r;
    std::vector<T> ohm;
    weighted_ribbon_graph(ribbon_graph _r) : r(_r), ohm(r.h, T(1)) {}
    weighted_ribbon_graph(ribbon_graph _r, std::vector<T> _ohm) : r(_r), ohm(_ohm) { }
    bool valid() {
        // Check if half edges have same resistance
        for(int i = 0; i < r.h; ++i) {
            if(ohm[i] != ohm[r.lam[i]]) return false;
        }
        return true;
    }
    weighted_ribbon_graph<T> dual() {
        std::vector<T> newohm(r.h);
        for(int i = 0; i < r.h; ++i) newohm[i] = T(1) / ohm[i];
        return weighted_ribbon_graph<T>(r.dual(), newohm);
    }
    weighted_graph<T> to_graph() {
        weighted_graph<T> g(r.v);
        for(int i = 0; i < r.h; ++i) {
            g[r.lam[i]].push_back(std::make_pair(ohm[i], r.lam[r.tau[i]]));
        }
        return g;
    }
    void set_weight(int i, T w) {
        ohm[i] = ohm[r.tau[i]] = w;
    }
};