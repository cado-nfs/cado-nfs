#include "cado.h"
#include "utils.h"
#include "fb-slice-weight.hpp"

template<typename FB_ENTRY_TYPE>
double fb_slice_weight_estimator<FB_ENTRY_TYPE>::operator()(sl_t const & sl, it_t a, it_t b) const
{
    if (verbose_output_get(0, 4, 0) != NULL)
        return compare(sl, a, b);
    return max(sl, a, b);
}

class wurst : public _padded_pod<wurst> {
  double worst_err, worst_est, worst_val, mse;
  unsigned long nr;
public:
  void update(const double est, const double val) {
    const double err = est / val - 1.;
    mse += err * err;
    nr++;
    if (fabs(err) > fabs(worst_err)) {
      worst_err = err;
      worst_est = est;
      worst_val = val;
    }
    verbose_output_print(0, 4, "%0.3g (%.3g)", val, err);
  }
  void print() {
    if (nr > 0)
      verbose_output_print(0, 4, "%0.3g vs. %0.3g (rel err. %.3g, MSE: %.3g)",
                           worst_est, worst_val, worst_err, mse / nr);
  }
};

static wurst worst_max, worst_avg, worst_simpson, worst_mertens;

void print_slice_weight_estimator_stats()
{
    verbose_output_start_batch();
    verbose_output_print(0, 4, "# Worst weight errors: max = ");
    worst_max.print();
    verbose_output_print(0, 4, ", avg = ");
    worst_avg.print();
    verbose_output_print(0, 4, ", simpson = ");
    worst_simpson.print();
    verbose_output_print(0, 4, ", mertens = ");
    worst_mertens.print();
    verbose_output_print(0, 4, "\n");
    verbose_output_end_batch();
}

template <class FB_ENTRY_TYPE>
double
fb_slice_weight_estimator<FB_ENTRY_TYPE>::compare(sl_t const & sl, const it_t start, const it_t end) const
{
  const double _max = max(sl, start, end),
               avg = avg(sl, start, end),
               simpson = simpson(sl, start, end),
               mertens = mertens(sl, start, end),
               _sum = exact(sl, start, end);
  verbose_output_start_batch();
  verbose_output_print(0, 4, "# Slice [%zu, %zu] weight: max = ", start, end);
  worst_max.update(_max, _sum);
  verbose_output_print(0, 4, ", avg = ");
  worst_avg.update(avg, _sum);
  verbose_output_print(0, 4, ", simpson = ");
  worst_simpson.update(simpson, _sum);
  verbose_output_print(0, 4, ", mertens = ");
  worst_mertens.update(mertens, _sum);
  verbose_output_print(0, 4, ", exact = %0.3g\n", _sum);
  verbose_output_end_batch();
  return _sum;
}

