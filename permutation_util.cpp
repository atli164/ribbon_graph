bool is_permutation(std::vector<int> pi) {
    sort(pi.begin(), pi.end());
    for(int i = 0; i < pi.size(); ++i)
        if(pi[i] != i) return false;
    return true;
}

int orbit_num(std::vector<int> pi) {
    unionfind uf(pi.size());
    for(int i = 0; i < pi.size(); ++i) {
        uf.unite(i, pi[i]);
    }
    return uf.c;
}

std::vector<int> composeperms(std::vector<int> alpha, std::vector<int> beta) {
    assert(alpha.size() == beta.size());
    std::vector<int> gamma(alpha.size());
    for(int i = 0; i < gamma.size(); ++i) {
        gamma[i] = alpha[beta[i]];
    }
    return gamma;
}

std::vector<std::vector<int>> factor_cycles(std::vector<int> pi) {
    std::vector<int> done(pi.size(), 0);
    std::vector<std::vector<int>> res;
    for(int i = 0; i < pi.size(); ++i) {
        if(done[i]) continue;
        std::vector<int> cyc;
        cyc.push_back(i);
        while(!done[cyc.back()]) {
            done[cyc.back()] = 1;
            cyc.push_back(pi[cyc.back()]);
        }
        res.push_back(cyc);
    }
    return res;
}