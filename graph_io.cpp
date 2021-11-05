#pragma once

graph read_graph(std::istream &in) {
    int n, d, x;
    in >> n;
    graph g(n);
    for(int i = 0; i < n; ++i) {
        in >> d;
        for(int j = 0; j < d; ++j) {
            in >> x;
            g[i].push_back(x);
        }
    }
    return g;
}

void write_graph(std::ostream &out, graph &g) {
    out << g.size() << '\n';
    for(int i = 0; i < g.size(); ++i) {
        out << g[i].size() << ' ';
        for(int x : g[i]) out << x << ' ';
        out << '\n';
    }
}

ribbon_graph read_ribbon_graph(std::istream &in) {
    int v, h;
    in >> v >> h;
    h *= 2;
    std::vector<int> lam(h), tau(h), rho(h);
    for(int i = 0; i < h; ++i) in >> lam[i];
    for(int i = 0; i < h; ++i) in >> tau[i];
    for(int i = 0; i < h; ++i) in >> rho[i];
    return ribbon_graph(v, h, lam, tau, rho);
}

void write_ribbon_graph(std::ostream &out, ribbon_graph& r) {
    out << r.v << ' ' << r.edges() << '\n';
    for(int i = 0; i < r.h; ++i) out << r.lam[i] << ' ';
    out << '\n';
    for(int i = 0; i < r.h; ++i) out << r.tau[i] << ' ';
    out << '\n';
    for(int i = 0; i < r.h; ++i) out << r.rho[i] << ' ';
    out << '\n';
}

void pretty_print_perm(std::ostream &out, std::vector<int>& pi) {
    std::vector<int> done(pi.size(), 0), cur;
    for(int i = 0; i < pi.size(); ++i) {
        if(done[i]) continue;
        cur.push_back(i);
        done[i] = 1;
        while(!done[pi[cur.back()]]) {
            cur.push_back(pi[cur.back()]);
            done[cur.back()] = 1;
        }
        out << '(';
        for(int j = 0; j < cur.size(); ++j) {
            out << cur[j];
            if(j != cur.size() - 1) out << ',';
        }
        out << ')';
        cur.clear();
    }
}

void pretty_print_ribbon_graph(std::ostream &out, ribbon_graph& r) {
    out << r.v << ' ' << r.edges() << '\n';
    for(int i = 0; i < r.h; ++i) out << r.lam[i] << ' ';
    out << '\n';
    pretty_print_perm(out, r.tau);
    out << '\n';
    pretty_print_perm(out, r.rho);
    out << '\n';
}
