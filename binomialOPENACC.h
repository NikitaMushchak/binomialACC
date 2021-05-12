//
// Created by nikita.mushak on 4/30/2021.
//

#ifndef APPS_BINOMIALOPENACC_H
#define APPS_BINOMIALOPENACC_H

#include <cstddef>
#include "/home/nikita.mushak/apps/apps/shared/pricing/models/impl/implied_volatility/AmericanOptionsImpl.h"
#include <shared/dividend_container/model/DividendContainer.h>
#include "/home/nikita.mushak/apps/apps/shared/pricing/models/impl/fd_settings/solver_common.h"
using std::size_t;

class IYieldCurve;

int BinomialWithDiscreteDivs_ACC(const int call_put_flag, const int exercise_style,
                                 const int payment_style, struct ModelResult & result,
                                 const double S, const double X, IYieldCurve & r_curve,
                                 IYieldCurve & q_curve, const double T, const double sigma,
                                 const dividend_container::DividendContainer & dividends,
                                 const bool use_smoothing, const char *& err_str);
#endif //APPS_BINOMIALOPENACC_H
