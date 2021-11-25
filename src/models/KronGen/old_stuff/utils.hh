
/// Calculates the required mean degree of a graph in order to generate a given
/// target clustering coefficient. The graph can either be complete or have
/// zero clustering.
/**
  * \param is_zero  Whether the graph in question has zero clustering or not
  * \param c_H      The clustering coefficient of the other Kronecker factor H
  * \param c        The target clustering coefficient
  * \param h        The mean degree of graph H
  * \param v        The degree distribution variance of H
*/
double get_mean_deg_c(const double c_H,
                      const double c,
                      const double h,
                      const double v)
{
    double g;

    // one Kronecker factor has clustering coefficient zero
    if (c_H >= c) {
        g = 3*pow(h, 2)*c_H-2*pow(h, 2)*c-3*h*c_H-h*c+6*h+3*v*c_H-2*c*v+c;
        g = pow(g, 2);
        g += -4*(-pow(h, 2)*c-2*h*c-c*v-c)*(pow(h, 2)*c_H-pow(h, 2)*c-h*c_H+h*c+c_H*v-c*v);
        g = sqrt(g);
        g += 3*pow(h, 2)*c_H -2*pow(h, 2)*c-3*h*c_H-h*c+6*h+3*c_H*v-2*c*v+c;
        g /= 2*c*(pow(h, 2)+2*h+v+1);
    }

    // one Kronecker factor is complete
    else {
        g = 8*pow(h, 2)*c_H+pow(h, 2)*c-9*pow(h, 2)-8*h*c_H+2*h*c+6*h+8*c_H*v-8*c*v+c-1;
        g *= c-1;
        g = sqrt(g);
        g += -2*pow(h, 2)*c_H+2*pow(h, 2)*c+2*h*c_H+h*c-3*h-2*c_H*v+2*c*v-c+1;
        g /= (2*(pow(h, 2)*c_H-pow(h, 2)*c-h*c_H-2*h*c+3*h+c_H*v-c*v-c+1));
    }

    // always select positive root
    return (fabs(g));

}


/// Return a list of possible vertex number factor pairs producing
// a desired product. Does not include (1, N) or (2, N/2)
vector_pt N_factors (const size_t N)
{
    vector<bool> candidates(round(static_cast<double>(N)/2.+2), true);
    candidates[2]=false;
    vector_pt res;
    for (int i = 2; i < round(static_cast<double>(N)/2+1); ++i) {
        if (not (N % i) && (candidates[i]==true)) {
            if (candidates[N/i] == true) {
                res.push_back({i, N/i});
                candidates[N/i] = false;
            }
        }
    }
    return res;
}

/// Return a list of possible vertex number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_N_factors (const size_t N)
{
    auto res = N_factors(N);
    size_t i = 1;
    while (res.empty()) {
        auto res2 = N_factors(N+i);
        res.insert(res.end(), res2.begin(), res2.end());
        if (N > i+2){
            res2=N_factors(N-i);
            res.insert(res.end(), res2.begin(), res2.end());
        }
        ++i;
    }
    return res;
}

/// Return a list of all possible mean degree factors producing a desired product
vector_pt mean_deg_factors (const size_t m)
{
    vector<bool> candidates(round(static_cast<double>(m)/2.+2), true);
    vector_pt res;
    candidates[1]=false;
    for (int i = 2; i < round(static_cast<double>(m)/2); ++i) {
        if (not ((m+1) % (i+1)) && (candidates[i]==true)) {
            auto k = ((m+1)/(i+1)-1);
            if (candidates[k] == true) {
                res.push_back({i, k});
                candidates[k] = false;
            }
        }
    }

    return res;
}

/// Return a list of possible mean degree number factor pairs producing a desired
/// product, or the closest possible factor pairs
vector_pt closest_mean_deg_factors (const size_t m)
{
    auto res = mean_deg_factors(m);
    size_t i = 1;
    while (res.empty()) {
        auto res2 = mean_deg_factors(m+i);
        res.insert(res.end(), res2.begin(), res2.end());
        if (m > i+3){
            res2=mean_deg_factors(m-i);
            res.insert(res.end(), res2.begin(), res2.end());
        }
        ++i;
    }

    return res;
}

/// Return two factors of a pair consisting of {num_vertices, mean_degree}
template<typename RNGType>
pair_pt get_factors_N_m (const size_t N,
                         const size_t m,
                         RNGType& rng)
{
    // Get the Kronecker factors of N, m and shuffle
    auto N_factors = Utils::closest_N_factors(N);
    auto m_factors = Utils::closest_mean_deg_factors(m);
    shuffle(N_factors.begin(), N_factors.end(), rng);
    shuffle(m_factors.begin(), m_factors.end(), rng);

    for (const auto& N_i : N_factors) {
        for (const auto& m_i : m_factors) {
            if ((N_i.first > m_i.first) and (N_i.second > m_i.second)) {
                return {{N_i.first, m_i.first}, {N_i.second, m_i.second}};
            }
            else if ((N_i.second > m_i.first) and (N_i.first > m_i.second)) {
                return {{N_i.second, m_i.first}, {N_i.first, m_i.second}};
            }
        }
    }

    return {{0, 0},{0, 0}};
}
