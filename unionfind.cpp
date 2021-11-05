#pragma once

struct unionfind {
    std::vector<int> p;
    int c;
    unionfind(int n) : p(n, -1), c(n) { }
    int find(int x) {
        return p[x] < 0 ? x : p[x] = find(p[x]);
    }
    bool united(int x, int y) {
        return find(x) == find(y);
    }
    void unite(int x, int y) {
        int a = find(x), b = find(y);
        if(a == b) return;
        if(p[a] > p[b]) std::swap(a, b);
        p[a] += p[b];
        p[b] = a;
        c--;
        return;
    }
    int size(int x) {
        return -p[find(x)];
    }
};