#ifndef GPUFIT_BIEXP_CUH_INCLUDED
#define GPUFIT_BIEXP_CUH_INCLUDED

__device__ void calculate_biexp(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    REAL * user_info_float = (REAL*) user_info;
    REAL x = 0;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }
  
    value[point_index] = p[0] * exp(p[1] * x) + p[2] * exp(p[3] * x);

    // derivatives

    REAL * current_derivative = derivative + point_index;

    current_derivative[0 * n_points] = exp(p[1] * x);
    current_derivative[1 * n_points] = x * p[0] * exp(p[1] * x);
    current_derivative[2 * n_points] = exp(p[3] * x);
    current_derivative[3 * n_points] = x * p[2] * exp(p[3] * x);
}

#endif
