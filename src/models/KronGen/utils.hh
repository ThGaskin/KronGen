#ifndef UTOPIA_MODELS_KRONGEN_UTILS
#define UTOPIA_MODELS_KRONGEN_UTILS

namespace Utopia::Models::KronGen::Utils {

/// Calculate the clustering coefficient of a Kronecker product of two graphs G, H
/** This product is symmetric in G, H
 * \param c_G        The clustering coefficient of the first graph
 * \param c_H        The clustering coefficient of the second graph
 * \param g          The mean degree of the first graph
 * \param h          The mean degree of the second graph
 * \param var_g      The degree distribution variance of the first graph
 * \param var_h      The degree distribution variance of the second graph
 *
 * \return c         The clustering coefficient of the Kronecker product
 */
double Kronecker_clustering (const double c_G,
                             const double c_H,
                             const double g,
                             const double h,
                             const double v_g,
                             const double v_h)
{
    double c = c_G*c_H*(v_g+g*(g-1))*(v_h+h*(h-1));
    c += 3*g*c_H*(v_h+h*(h-1))+3*h*c_G*(v_g+g*(g-1));
    c += c_G*(v_g+g*(g-1))+c_H*(v_h+h*(h-1));
    c += 6*g*h;
    c /= (v_g*v_h+pow(1+g, 2)*v_h+pow(1+h, 2)*v_g+(g*h+g+h)*(g*h+g+h-1));

    return c;
}

// Calculate the mean_degree of a graph with a given clustering coefficient
// To do: test me
std::size_t get_mean_deg_c(const bool is_zero,
                           const double c_0, // c of Kronecker factor
                           const double c_1, // desired clustering coefficient
                           const double h, // mean degree of Kronecker factor
                           const double v) // variance of degree distribution
{
    // one Kronecker factor has clustering coefficient zero
    //solve (v(1+x)^2 + (xh+x+h)(xh+h+x-1))*q = p*(1+3x)(v+h^2-h)+6xh for x
    if (is_zero) {

      double g = 3*pow(h, 2)*c_0-2*pow(h, 2)*c_1-3*h*c_0-h*c_1+6*h+3*v*c_0-2*c_1*v+c_1;
      g = pow(g, 2);
      g += -4*(-pow(h, 2)*c_1-2*h*c_1-c_1*v-c_1)*(pow(h, 2)*c_0-pow(h, 2)*c_1-h*c_0+h*c_1+c_0*v-c_1*v);
      g = sqrt(g);
      g += 3*pow(h, 2)*c_0 -2*pow(h, 2)*c_1-3*h*c_0-h*c_1+6*h+3*c_0*v-2*c_1*v+c_1;
      g /= 2*c_1*(pow(h, 2)+2*h+v+1);

      std::cout<<"zero: "<<g<<std::endl;

      g = std::round(g);

      return g;

    }

    // one Kronecker factor is complete
    // solve (v(1+x)^2 + (xh+x+h)(xh+h+x-1))*q = p*(1+x)^2(v+h^2-h)+3hx^2+3xh+x^2-x for x
    else {

      double g = 8*pow(h, 2)*c_0+pow(h, 2)*c_1-9*pow(h, 2)-8*h*c_0+2*h*c_1+6*h+8*c_0*v-8*c_1*v+c_1-1;
      g *= c_1-1;
      g = sqrt(g);
      g += -2*pow(h, 2)*c_0+2*pow(h, 2)*c_1+2*h*c_0+h*c_1-3*h-2*c_0*v+2*c_1*v-c_1+1;
      g /= (2*(pow(h, 2)*c_0-pow(h, 2)*c_1-h*c_0-2*h*c_1+3*h+c_0*v-c_1*v-c_1+1));

      std::cout<<"complete: "<<g<<std::endl;

      // always select positive root
      return std::round(fabs(g));

    }
}

}
#endif
