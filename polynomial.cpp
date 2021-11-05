#pragma once

template<int n, typename T> struct polynom;

template <typename T> struct polynom<0, T> {
    T x;
    polynom(const T& y = T(0)) : x(y) { }
    bool is_zero() const { return x == T(0); }
    operator const T&() const { return x; }
    polynom& operator=(const T& y) { x = y; return *this; }
    T operator()() const { return x; }
    polynom operator+(const T& y) const { return polynom(x + y); }
    polynom operator-(const T& y) const { return polynom(x - y); }
    polynom operator*(const T& y) const { return polynom(x * y); }
    polynom operator/(const T& y) const { return polynom(x / y); }
    polynom operator-() const { return polynom(-x); }
    polynom& operator+=(const T& y) { x += y; return *this; }
    polynom& operator-=(const T& y) { x -= y; return *this; }
    polynom& operator*=(const T& y) { x *= y; return *this; }
    polynom& operator/=(const T& y) { x /= y; return *this; }
    bool operator==(const T& y) const { return x == y; }
    bool operator!=(const T& y) const { return x != y; }
    void print(std::ostream &s, int m = 0) const { s << x; }
    void swap(polynom &p) { std::swap(x, p.x); }
};

template<int n, typename T>
struct polynom {
    static polynom<n - 1, T> zero_val;
    std::vector<polynom<n - 1, T>> coeff;
    
    polynom() : coeff() { }
    polynom(const T& x) : coeff(1, polynom<n - 1, T>(x)) { }
    polynom(const polynom<n - 1, T>& p) : coeff(1, p) { }
    template<typename S>
    polynom(const polynom<n, S>& p) : coeff(p.degree() + 1, polynom<n - 1, T>()) {
        for(int i = 0; i < coeff.size(); ++i) coeff[i] = p[i];
        normalize();
    }
    template<typename S>
    polynom& operator=(const polynom<n, S>& p) {
        coeff.resize(p.degree() + 1);
        for(int i = 0; i < coeff.size(); ++i) coeff[i] = p[i];
        normalize();
        return *this;
    }
    void normalize() {
        int sz = coeff.size();
        while(sz && coeff[sz - 1].is_zero()) sz--;
        if(sz != coeff.size()) coeff.resize(sz);
    }
    void grow_to(int sz) {
        if(coeff.size() < sz) coeff.resize(sz);
    }
    int degree() const { return coeff.size() - 1; }
    const polynom<n - 1, T>& leading() const { return coeff.size() ? coeff.back() : zero_val; }
    bool is_zero() const { return coeff.size() == 0; }
    polynom<n - 1, T>& operator[](int i) {
        grow_to(i + 1);
        return coeff[i];
    }
    const polynom<n - 1, T>& operator[](int i) const {
        return i < coeff.size() ? coeff[i] : zero_val;
    }
    polynom operator*(const T& x) const {
        polynom p(*this);
        for(int i = 0; i < coeff.size(); ++i) p[i] *= x;
        return p;
    }
    polynom operator/(const T& x) const {
        polynom p(*this);
        for(int i = 0; i < coeff.size(); ++i) p[i] /= x;
        return p;
    }
    polynom& operator*=(const polynom& p) {
        return *this = *this * p;
    }
    polynom& operator*=(const T& x) {
        for(int i = 0; i < coeff.size(); ++i) coeff[i] *= x;
        return *this;
    }
    polynom& operator/=(const T& x) {
        for(int i = 0; i < coeff.size(); ++i) coeff[i] /= x;
        return *this;
    }
    friend polynom operator*(const T& x, const polynom& p) {
        polynom q(p);
        for(int i = 0; i < p.coeff.size(); ++i) q[i] *= x;
        return q;
    }
    polynom operator-() const {
        polynom p(*this);
        for(int i = 0; i < coeff.size(); ++i) p[i] = -p[i];
        return p;
    }
    polynom operator+(const polynom& p) const {
        polynom q(*this);
        q.grow_to(p.coeff.size());
        for(int i = 0; i < p.coeff.size(); ++i) q[i] += p[i];
        q.normalize();
        return q;
    }
    polynom operator-(const polynom& p) const {
        polynom q(*this);
        q.coeff.grow_to(p.coeff.size());
        for(int i = 0; i < p.coeff.size(); ++i) q[i] -= p[i];
        q.normalize();
        return q;
    }
    polynom& operator+=(const polynom& p) {
        grow_to(p.coeff.size());
        for(int i = 0; i < p.coeff.size(); ++i) coeff[i] += p[i];
        normalize();
        return *this;
    }
    polynom& operator-=(const polynom& p) {
        grow_to(p.coeff.size());
        for(int i = 0; i < p.coeff.size(); ++i) coeff[i] -= p[i];
        normalize();
        return *this;
    }
    polynom operator+(const polynom<n - 1, T>& p) const {
        polynom q(*this);
        q.grow_to(1);
        q[0] += p;
        q.normalize();
        return q;
    }
    polynom operator-(const polynom<n - 1, T>& p) const {
        polynom q(*this);
        q.grow_to(1);
        q[0] -= p;
        q.normalize();
        return q;
    }
    polynom operator+(const T& x) const {
        polynom p(*this);
        p.grow_to(1);
        p[0] += x;
        p.normalize();
        return p;
    }
    polynom operator-(const T& x) const {
        polynom p(*this);
        p.grow_to(1);
        p[0] -= x;
        p.normalize();
        return p;
    }
    polynom& operator+=(const polynom<n - 1, T>& p) {
        grow_to(1);
        coeff[0] += p;
        normalize();
        return *this;
    }
    polynom& operator-=(const polynom<n - 1, T>& p) {
        grow_to(1);
        coeff[0] -= p;
        normalize();
        return *this;
    }
    polynom& operator+=(const T& x) {
        grow_to(1);
        coeff[0] += x;
        normalize();
        return *this;
    }
    polynom& operator-=(const T& x) {
        grow_to(1);
        coeff[0] -= x;
        normalize();
        return *this;
    }
    bool operator==(const polynom& p) const {
        if(coeff.size() != p.coeff.size()) return false;
        for(int i = 0; i < coeff.size(); ++i) {
            if(coeff[i] != p.coeff[i]) return false;
        }
        return true;
    }
    bool operator!=(const polynom& p) const { return !(*this == p); }
    bool operator==(const T& x) const {
        if((x == T(0)) && coeff.empty()) return true;
        if(coeff.size() != 1) return false;
        return coeff[0] == x;
    }
    bool operator!=(const T& x) const { return !(*this == x); }
    void print(std::ostream &s, int m = n) const {
        if(is_zero()) s << T(0);
        else {
            int nonz = 0;
            for(int i = 0; i < coeff.size(); ++i) {
                if(!coeff[i].is_zero()) nonz++;
            }
            if(nonz > 1 && m != n) s << '(';
            bool fst = true;
            for(int i = 0; i < coeff.size(); ++i) {
                if(coeff[i].is_zero()) continue;
                if(fst) fst = false;
                else s << " + ";
                coeff[i].print(s, m);
                if(i > 0) {
                    s << '*';
                    if(m <= 26) s << (char)('A' + (m - n));
                    else s << "X_" << m - n;
                    if(i > 1) s << "^" << i;
                }
            }
            if(nonz > 1 && m != n) s << ')';
        }
    }
    friend std::ostream& operator<<(std::ostream& s, const polynom& p) { p.print(s); return s; }
    polynom operator*(const polynom & p) const {
        int d = std::max(p.degree() + degree(), -1);
        polynom r;
        r.grow_to(d + 1);
        if(!is_zero() && !p.is_zero())
            for(int i = 0; i < r.coeff.size(); ++i)
                for(int j = 0; j < coeff.size(); ++j)
                    if(i < j + p.coeff.size())
                        r[i] += coeff[j] * p[i - j];
        r.normalize();
        return r;
    }
    void swap(polynom & p) {
        coeff.swap(p.coeff);
    }
    polynom pow(int e) {
        if(e == 0) return polynom(T(1));
        polynom res(*this);
        for(int i = 0; i < e - 1; ++i) {
            res = res * (*this);
        }
        return res;
    }
    polynom<n - 1, T> operator()(T x) const {
        polynom<n - 1, T> res;
        T varmul = T(1);
        for(int i = 0; i < coeff.size(); ++i) {
            res += coeff[i] * varmul;
            varmul *= x;
        }
        return res; 
    }
};

template<int n, typename T>
polynom<n - 1, T> polynom<n, T>::zero_val;
