ribbon_graph squaretiles(int n) {
    int v = 2 * n + 2;
    int e = 2 * 8 * n;
    std::vector<int> lam(e), tau(e), rho(e);
    for(int i = 0; i < e; ++i) {
        tau[i] = i ^ 1;
    }
    for(int i = 0; i < n; ++i) {
        lam[4 * i] = lam[4 * i + 2] = 0;
        lam[4 * i + 1] = lam[4 * i + 3] = i + 1;
    }
    for(int i = 0; i < n; ++i) {
        lam[4 * (i + n)] = lam[4 * (i + n) + 2] = n + 1;
        lam[4 * (i + n) + 1] = lam[4 * (i + n) + 3] = i + n + 2;
    }
    for(int i = 0; i < n; ++i) {
        lam[4 * (i + 2 * n)] = lam[4 * (i + 2 * n) + 2] = i + 1;
        lam[4 * (i + 2 * n) + 1] = lam[4 * (i + 2 * n) + 3] = i + n + 2;
    }
    for(int i = 0; i < n; ++i) {
        lam[4 * (i + 3 * n)] = lam[4 * (i + 3 * n) + 2] = 0;
        lam[4 * (i + 3 * n) + 1] = lam[4 * (i + 3 * n) + 3] = n + 1;
    }
    for(int i = 0; i < n; ++i) {
        rho[4 * i + 1] = 4 * (i + 2 * n);
        rho[4 * (i + 2 * n)] = 4 * i + 3;
        rho[4 * i + 3] = 4 * (i + 2 * n) + 2;
        rho[4 * (i + 2 * n) + 2] = 4 * i + 1;
    }
    for(int i = 0; i < n; ++i) {
        rho[4 * (i + 2 * n) + 3] = 4 * (i + n) + 3;
        rho[4 * (i + n) + 3] = 4 * (i + 2 * n) + 1;
        rho[4 * (i + 2 * n) + 1] = 4 * (i + n) + 1;
        rho[4 * (i + n) + 1] = 4 * (i + 2 * n) + 3;
    }
    for(int i = 0; i < n; ++i) {
        rho[4 * (i + 3 * n)] = 4 * i;
        rho[4 * (i + 3 * n) + 2] = 4 * i + 2;
        rho[4 * i] = 4 * (i + 3 * n) - 2;
        rho[4 * i + 2] = 4 * (i + 3 * n);
    }
    rho[0] = 16 * n - 2;
    for(int i = 0; i < n; ++i) {
        rho[4 * (i + n)] = 4 * (i + 3 * n) + 1;
        rho[4 * (i + n) + 2] = 4 * (i + 3 * n) + 3;
        rho[4 * (i + 3 * n) + 1] = 4 * (i + n) + 2;
        rho[4 * (i + 3 * n) + 3] = 4 * (i + n) + 4;
    }
    rho[16 * n - 1] = 4 * n;
    return ribbon_graph(v, e, lam, tau, rho);
}

std::vector<int> verts2rho(std::vector<std::vector<int>>& verts, int e) {
    std::vector<int> rho(e);
    for(int i = 0; i < verts.size(); ++i) {
        for(int j = 0; j < verts[i].size() - 1; ++j) {
            rho[verts[i][j]] = verts[i][j + 1];
        }
        rho[verts[i][verts[i].size() - 1]] = verts[i][0];
    }
    return rho;
}

ribbon_graph brute_rho(int v, int e, std::vector<int> lam, std::vector<int> tau, int g) {
    std::vector<std::vector<int>> verts(v);
    for(int i = 0; i < e; ++i) {
        verts[lam[i]].push_back(i);
    }
    std::random_device rd;
    std::mt19937 rng(rd());
    while(true) {
        for(int i = 0; i < v; ++i) {
            std::shuffle(verts[i].begin(), verts[i].end(), rng);
        }
        std::vector<int> rho = verts2rho(verts, e);
        ribbon_graph r(v, e, lam, tau, rho);
        if(r.genus() == g && r.smooth()) {
            return r;
        }
    }
}

double curtime() { return static_cast<double>(clock()) / CLOCKS_PER_SEC; }
ribbon_graph anneal_rho(int v, int e, std::vector<int> lam, std::vector<int> tau, int g, double seconds) {
    std::vector<std::vector<int>> verts(v);
    for(int i = 0; i < e; ++i) {
        verts[lam[i]].push_back(i);
    }
    std::vector<int> to_shuffle;
    for(int i = 0; i < v; ++i) {
        if(verts[i].size() > 1) to_shuffle.push_back(i);
    }
    std::random_device rd;
    std::mt19937 rng(rd());
    for(int i = 0; i < v; ++i) {
        std::shuffle(verts[i].begin(), verts[i].end(), rng);
    }
    std::vector<int> rho = verts2rho(verts, e);
    ribbon_graph prv(v, e, lam, tau, rho);
    if(to_shuffle.empty()) return prv;
    int score = prv.genus(), iters = 0;
    double T0 = 100.0, T1 = 0.001, progress = 0, temp = T0, starttime = curtime();
    std::uniform_real_distribution<double> randfloat(0.0, 1.0);
    std::uniform_int_distribution<int> randvert(0, to_shuffle.size() - 1);
    while(true) {
        if(!(iters & ((1 << 4) - 1))) {
            progress = (curtime() - starttime) / seconds;
            temp = T0 * pow(T1 / T0, progress);
            if(progress > 1.0) break; 
        }
        int a = to_shuffle[randvert(rng)];
        std::uniform_int_distribution<int> randind(0, verts[a].size() - 2);
        int b = randind(rng);
        std::swap(verts[a][b], verts[a][b + 1]);
        rho = verts2rho(verts, e);
        ribbon_graph cur(v, e, lam, tau, rho);
        int newscore = cur.genus();
        if(newscore <= score || randfloat(rng) < exp((newscore - score) / temp)) {
            prv = cur;
            score = newscore;
            if(score == g) return prv;
        } else {
            std::swap(verts[a][b], verts[a][b + 1]);
        }
        iters++; 
    }
    return prv; 
}

std::pair<std::vector<int>, std::vector<int>> lamtau(graph& G) {
    int e = 0;
    for(int i = 0; i < G.size(); ++i) {
        e += G[i].size();
    }
    std::vector<int> tau(e), lam(e);
    for(int i = 0; i < e; ++i) {
        tau[i] = i ^ 1;
    }
    std::vector<std::pair<int,int>> edg;
    for(int i = 0; i < G.size(); ++i) {
        int loopcnt = 0;
        for(int j : G[i]) {
            if(i < j) edg.push_back(std::make_pair(i, j));
            else if(i == j) loopcnt++;
        }
        for(int j = 0; j < loopcnt / 2; ++j) edg.push_back(std::make_pair(i, i));
    }
    for(int i = 0; i < edg.size(); ++i) {
        lam[2 * i] = edg[i].first;
        lam[2 * i + 1] = edg[i].second;
    }
    return std::make_pair(lam, tau);
}

ribbon_graph bogo_embed(graph& G, int g) {
    auto lt = lamtau(G);
    return brute_rho(G.size(), lt.first.size(), lt.first, lt.second, g);
}

ribbon_graph annealed_embed(graph& G, int g, double seconds) {
    auto lt = lamtau(G);
    return anneal_rho(G.size(), lt.first.size(), lt.first, lt.second, g, seconds);
}
