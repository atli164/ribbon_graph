#pragma once

void add_edge(graph& g, int i, int j) {
    g[i].push_back(j);
    g[j].push_back(i);
}

void delete_edge(graph& g, int x, int y) {
    if(x < 0 || x >= g.size()) return;
    if(y < 0 || y >= g.size()) return;
    for(int i = 0; i < g[x].size(); ++i) {
        if(g[x][i] == y) {
            std::swap(g[x][i], g[x].back());
            g[x].pop_back();
            break;
        }
    }
    for(int i = 0; i < g[y].size(); ++i) {
        if(g[y][i] == x) {
            std::swap(g[y][i], g[y].back());
            g[y].pop_back();
            break;
        }
    }
}

template<typename T>
void add_edge(weighted_graph<T>& g, int i, int j, T w) {
    g[i].push_back(std::make_pair(w, j));
    g[j].push_back(std::make_pair(w, i));
}

template<typename T>
void delete_edge(weighted_graph<T>& g, int x, int y) {
    if(x < 0 || x >= g.size()) return;
    if(y < 0 || y >= g.size()) return;
    for(int i = 0; i < g[x].size(); ++i) {
        if(g[x][i].second == y) {
            std::swap(g[x][i], g[x].back());
            g[x].pop_back();
            break;
        }
    }
    for(int i = 0; i < g[y].size(); ++i) {
        if(g[y][i].second == x) {
            std::swap(g[y][i], g[y].back());
            g[y].pop_back();
            break;
        }
    }
}

template<typename T>
weighted_graph<T> unit_weights(graph& g) {
    weighted_graph<T> h(g.size());
    for(int i = 0; i < g.size(); ++i) {
        for(int x : g[i]) {
            h[i].push_back(std::make_pair(T(1), x));
        }
    }
    return h;
}

template<typename T>
void change_weight(weighted_graph<T>& g, T w, int x, int y) {
    if(x < 0 || x >= g.size()) return;
    if(y < 0 || y >= g.size()) return;
    for(int i = 0; i < g[x].size(); ++i) {
        if(g[x][i].second == y) {
            g[x][i].first = w;
            break;
        }
    }
    for(int i = 0; i < g[y].size(); ++i) {
        if(g[y][i].second == x) {
            g[y][i].first = w;
            break;
        }
    }
}

template<typename T>
graph forget_weights(weighted_graph<T>& g) {
    graph h(g.size());
    for(int i = 0; i < g.size(); ++i) {
        for(auto p : g[i]) {
            h[i].push_back(p.second);
        }
    }
    return h;
}

std::ostream& operator <<(std::ostream& out, graph& g) {
    for(int i = 0; i < g.size(); ++i) {
        out << i << ": ";
        for(int x : g[i]) out << x << ' ';
        out << '\n';
    }
    return out;
}