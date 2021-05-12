//
// Created by nikita.mushak on 4/30/2021.
//

#include "binomialOPENACC.h"
//
// Created by nikita.mushak on 4/30/2021.
//
#include <shared/interest_rate/IYieldCurve.h>
//#include <shared/math/NewtonBrentSolver.hpp>
#include <shared/math/TbMinMax.h>
#include <shared/math/linear_interpolation.h>

#include "shared/math/math_constants.h"

#include "/home/nikita.mushak/apps/apps/shared/pricing/models/impl/BlackScholes.h"
#include "shared/math/probability_distributions.h"
#include "/home/nikita.mushak/apps/apps/shared/pricing/models/impl/fd_settings/solver_common.h"
#include "/home/nikita.mushak/apps/apps/shared/pricing/models/ModelResult.h"
#include "/home/nikita.mushak/apps/apps/shared/pricing/models/InstrumentBase.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <stdio.h>

using namespace pricing_model::option_properties;
using namespace tb_min_max;
using namespace dividend_container;

#ifndef DIVIDEND_TREATED_AS_ZERO
#define DIVIDEND_TREATED_AS_ZERO 1e-10
#endif

constexpr double trembling = 1.0 + 1e-6;


bool FindExerciseBoundary(const int call_put_flag, const double strike_price,
                          const int step, const double * const prices, const double * const values, const double pv,
                          double & exercise_boundary, ExerciseBoundaryState & exercise_boundary_state)
{
    int i;
    if (call_put_flag == CP_CALL) {
        for(i = step; i >= 0 && pv * values[i] < trembling * (prices[i] - strike_price); --i) {}
    } else {
        for(i = 0; i <= step && pv * values[i] < trembling * (strike_price - prices[i]); ++i) {}
    }

    if (call_put_flag == CP_CALL && i == step) {
        exercise_boundary_state = ExerciseBoundaryState::ABOVE;
        exercise_boundary = prices[step];
        return false;
    } else if (call_put_flag == CP_CALL && i == -1) {
        exercise_boundary_state = ExerciseBoundaryState::BELOW;
        exercise_boundary = prices[0];
        return false;
    }

    if (call_put_flag == CP_PUT && i == 0) {
        exercise_boundary_state = ExerciseBoundaryState::BELOW;
        exercise_boundary = prices[0];
        return false;
    } else if (call_put_flag == CP_PUT && i == step + 1) {
        exercise_boundary_state = ExerciseBoundaryState::ABOVE;
        exercise_boundary = prices[step];
        return false;
    }

    constexpr double tolerance = 1e-6;

    /*
     * expressing S from the equation
     * for Calls:
     * trembling * (S - X) = pv * (V_i + (V_{i+1} - V_i) * (S - S_i) / (S_{i+1} - S_i)) =>
     * S = (trembling * X + pv * V_{i-1} - pv * S_{i-1} * dvds) / (trembling - pv * dvds)
     *
     * for Puts:
     * trembling * (X - S) = pv * (V_i + (V_{i+1} - V_i) * (S - S_i) / (S_{i+1} - S_i)) =>
     * S = (trembling * X - pv * V_{i-1} + pv * S_{i-1} * dvds) / (trembling + pv * dvds)
     */
    if (call_put_flag == CP_CALL) {
        const double dvds = (values[i + 1] - values[i]) / (prices[i + 1] - prices[i]);
        const double spot_multiplier = (trembling - pv * dvds);

        if (std::abs(spot_multiplier) < tolerance) {
            exercise_boundary = values[i + 1];
            return true;
        }

        const double V_minus_S_dVdS = pv * (values[i] - dvds * prices[i]);
        exercise_boundary = (trembling * strike_price + V_minus_S_dVdS) / spot_multiplier;
    } else {
        const double dvds = (values[i] - values[i - 1]) / (prices[i] - prices[i - 1]);
        const double spot_multiplier = (trembling + pv * dvds);

        if (std::abs(spot_multiplier) < tolerance) {
            exercise_boundary = values[i - 1];
            return true;
        }

        const double V_minus_S_dVdS = pv * (values[i - 1] - dvds * prices[i - 1]);
        exercise_boundary = (trembling * strike_price - V_minus_S_dVdS) / spot_multiplier;
    }

    exercise_boundary_state = ExerciseBoundaryState::TOUCHED;

    return true;
}



int BinomialWithDiscreteDivs_ACC(const int call_put_flag, const int exercise_style,
                                 const int payment_style, struct ModelResult & result,
                                 const double S, const double X, IYieldCurve & r_curve,
                                 IYieldCurve & q_curve, const double T, const double sigma,
                                 const DividendContainer & dividends,
                                 const bool use_smoothing, const char *& err_str)
{

    constexpr int binomial_steps = 1024;
    constexpr int steps_to_prolong = 256;
    if ( (check_expired(T) == EXPIRED) || (check_expired(T) == AT_EXPIRATION) )
    {
        if (call_put_flag == CP_CALL)
        {
            result.fair_value = tbmax(S - X, 0.0);
            if (std::abs(S - X) < 1e-8) {
                result.delta = 0.5;
            }
            else {
                result.delta = ((S > X) ? 1.0 : 0.0);
            }
            result.gamma = 0.0;
            result.vega = 0.0;
#ifdef CONTINUOUS_TIME_DERIVATIVES
            result.theta = 0.0;
            result.charm = 0.0;
#endif
            result.rho = 0.0;
            result.vanna = 0.0;
            result.vomma = 0.0;
        }
        else if (call_put_flag == CP_PUT)
        {
            result.fair_value = tbmax(X - S, 0.0);
            if (std::abs(S - X) < 1e-8) {
                result.delta = -0.5;
            }
            else {
                result.delta = ( (S < X) ? -1.0 : 0.0);
            }
            result.gamma = 0.0;
            result.vega = 0.0;
#ifdef CONTINUOUS_TIME_DERIVATIVES
            result.theta = 0.0;
            result.charm = 0.0;
#endif
            result.rho = 0.0;
            result.vanna = 0.0;
            result.vomma = 0.0;
        }
        else {
            err_str = "Bad call-put type in binomial model";
            return EC_BAD_CallPutType;
        }
        return EC_SUCCESS;
    }

    if (S < 1e-8)
    {
        if (call_put_flag == CP_CALL)
        {
            result.fair_value = 0.0;
            result.delta = 0.0;
            result.gamma = 0.0;
            result.vega = 0.0;
#ifdef CONTINUOUS_TIME_DERIVATIVES
            result.theta = 0.0;
            result.charm = 0.0;
#endif
            result.rho = 0.0;
            result.vanna = 0.0;
            result.vomma = 0.0;
        }
        else if (call_put_flag == CP_PUT)
        {
            if (exercise_style == STYLE_EUROPEAN) {
                result.fair_value = X * exp(-r_curve.GetAverageValue(T) * T);
            }
            else {
                result.fair_value = X;
            }
            result.delta = -1.0;
            result.gamma = 0.0;
            result.vega = 0.0;
#ifdef CONTINUOUS_TIME_DERIVATIVES
            result.theta = 0.0;
            result.charm = 0.0;
#endif
            result.rho = 0.0;
            result.vanna = 0.0;
            result.vomma = 0.0;
        }
        else {
            err_str = "Bad call-put type in binomial model";
            return EC_BAD_CallPutType;
        }
        return EC_SUCCESS;
    }

    // check whether q is zero
    bool q_is_zero = q_curve.IsZero(T);

    bool all_divs_counted = false;

    if (dividends.Empty())
    {
//        std::cout<<"dividends.Empty() \n";
        all_divs_counted = true;

        if (q_is_zero && !r_curve.IsNegative(0, T) && call_put_flag == CP_CALL)
        {
            double r = r_curve.GetAverageValue(T);
            double q = q_curve.GetAverageValue(T);
            std::cout<<"Black Scholes!\n";
            return BlackScholes_Standard(CP_CALL, payment_style, &result, S, X,
                                         r, q, T, sigma, err_str);
        }

        // r_adj = q_adj = 0 at maturity due to zero rate_time
        if (q_is_zero && r_curve.IsZero(T) && call_put_flag == CP_PUT)
        {
            return BlackScholes_Standard(CP_PUT, payment_style, &result, S, X,
                                         0.0, 0.0, T, sigma, err_str);
        }
    }

    /*
     * The comparison aginst 401% is intended
     * taking into account calculation of the vega which is implemented via finite difference
     */
    if (sigma > 4.01) {
        err_str = "Adjusted volatility too big. Maximum allowed value is 400%";
        return EC_BAD_Volatility;
    }

    if (sigma < 1e-4) {
        err_str = "Adjusted volatility too small. Minimum allowed value is 0.01%";
        return EC_BAD_Volatility;
    }

    int i_curr_div = dividends.Count() - 1;

    int add_steps = 0;

    // 1) Find the step into which the first dividend falls
    // 2) Find discounted dividend
    //    D' = D * exp(-r * (t_div - t_first_div_step))
    // 3) Find lower boundary S'' = S * d^(first_div_step + 2)
    // 4) if S'' < D' then add_steps = steps_to_prolong and proceed as usual
    // 5) otherwise find M such that M = 2*M', 1<=M'<=steps_to_prolong/2, S''-D' at first_div_step > S'' * d^(2*M'):
    //    M' = min( ceil(0.5 * log_d(1 - D'/S'')), steps_to_prolong/2)
    // 6) Add M steps

    int i;
    bool only_percentage_divs = dividends.HasPercentageDividends();
    if (!dividends.Empty() && !dividends.DividendSignsEqual())
    {
        err_str = "Different signs for discrete and percentage dividends";
        return EC_FAILURE;
    }

    if (!dividends.Empty())
    {
        double dt_init = T / binomial_steps;
        double u_init = exp(sigma * sqrt(dt_init));
        double d_init = 1.0 / u_init;

        double first_div_time = dividends.GetTime(0);
        int first_div_step = static_cast<int>(std::floor(first_div_time / dt_init));
        double S_lower = S * std::pow(d_init, first_div_step + 2);
        double S_upper = S * std::pow(u_init, first_div_step + 2);

        bool negative_dividend = false;
        double D_disc = 0.0;

        D_disc = dividends.GetAbsolute(0) * exp(-r_curve.GetTimeIntegral(dt_init * first_div_step, first_div_time));
        if (D_disc < (-DIVIDEND_TREATED_AS_ZERO) ) {
            negative_dividend = true;
        }

        if (dividends.GetPercentage(0) < (-DIVIDEND_TREATED_AS_ZERO) )
        {
            negative_dividend = true;

            if (D_disc < (-DIVIDEND_TREATED_AS_ZERO)) {
                D_disc = tbmax(D_disc, dividends.GetPercentage(0) * S_upper);
            }
            else {
                D_disc = dividends.GetPercentage(0) * S_upper;
            }
        }
        else if (dividends.GetPercentage(0) > DIVIDEND_TREATED_AS_ZERO)
        {
            if (D_disc > DIVIDEND_TREATED_AS_ZERO) {
                D_disc = tbmin(D_disc, dividends.GetPercentage(0) * S_lower);
            }
            else {
                D_disc = dividends.GetPercentage(0) * S_lower;
            }
        }

        if (S_lower < D_disc) {
            add_steps = steps_to_prolong;
        } else if (negative_dividend) {
            int Mh_max = steps_to_prolong / 2;
            int Mh = tbmin(static_cast<int>(std::ceil(0.5 * (std::log(1.0 - D_disc/S_upper) / std::log(u_init)))), Mh_max);
            add_steps = 2*Mh;
        } else {
            int Mh_max = steps_to_prolong / 2;
            int Mh = tbmin(static_cast<int>(std::ceil(0.5 * (std::log(1.0 - D_disc/S_lower) / std::log(d_init)))), Mh_max);
            add_steps = 2*Mh;
        }
    }
    else {
        add_steps = 2;
    }

    add_steps = tbmax(add_steps, 2);

    double add_t = add_steps * T / binomial_steps;

    double t = T + add_t;
    int n_steps = binomial_steps + add_steps;

    bool move_term_condition = false;
//    if (callPutFlag == CALL_OPTION && q_is_zero) {
//        move_term_condition = true;
//    }

    if ( !dividends.Empty() && move_term_condition)
    {
        t = dividends.GetTime(dividends.Count()-1) + add_t;
    }

    double dt = t / n_steps;
    double u = exp(sigma * sqrt(dt));
    double uu = u * u;
    double d = 1.0 / u;

    double *prices = NULL;
    double option_vals[binomial_steps + steps_to_prolong + 1];
    double prices_even[binomial_steps + steps_to_prolong + 1];
    double prices_odd[binomial_steps + steps_to_prolong + 1];


    int step;

    int start_roll_down = 1;
    if (use_smoothing) {
        start_roll_down = 2;
    }

    double cur_div_absolute = 0.0;
    double cur_div_percentage = 0.0;
    double cur_div_mul_percentage = 1.0;
    // need to check whether there are dividends after start_roll_down step
    if (!dividends.Empty())
    {
        double roll_down_time_boundary = t - (start_roll_down - 1) * dt - 1e-3 * dt;
        for (;(i_curr_div >= 0) &&
              (dividends.GetTime(i_curr_div) + add_t > roll_down_time_boundary) &&
              (!all_divs_counted);)
        {
            cur_div_absolute += dividends.GetAbsolute(i_curr_div);
            cur_div_percentage += dividends.GetPercentage(i_curr_div);

            if (only_percentage_divs) {
                cur_div_mul_percentage *= 1 - dividends.GetPercentage(i_curr_div);
            }

            if (i_curr_div > 0) {
                i_curr_div--;
            }
            else {
                all_divs_counted = true;
            }
        }
    }

    double abs_cur_div_absolute = std::abs(cur_div_absolute);
    double abs_cur_div_percentage = std::abs(cur_div_percentage);

    double option_values_interpolated[binomial_steps + steps_to_prolong + 1];

    prices_even[0] = S * std::pow(d, n_steps);

    if (only_percentage_divs)
    {
        double percentage_mul_all_layers = 1.0;
        for (size_t idx = 0; idx < dividends.Count(); ++idx) {
            percentage_mul_all_layers *= 1 - dividends.GetPercentage(idx);
        }
        prices_even[0] *= percentage_mul_all_layers;
    }

    for(i = 1; i <= n_steps; ++i) {
        prices_even[i] = uu * prices_even[i-1];
    }

#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
    for(i = 0; i <= n_steps; ++i) {
        prices_odd[i] = prices_even[i] * u;
    }

    if (!use_smoothing) //if (!useSmoothing || !q_is_zero)
    {
        prices = prices_even;

        if (call_put_flag == CP_CALL)
        {
            if (exercise_style == STYLE_EUROPEAN)
            {
                //Final condition is right-shifted by cur_div_absolute if cur_div_absolute!=0.
                if (abs_cur_div_percentage <= DIVIDEND_TREATED_AS_ZERO) //cur_div_absolute != 0
                {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= n_steps; ++i) {
                        option_vals[i] = tbmax(prices[i] - (X + cur_div_absolute), 0.0 );
                    }
                }
                else if (abs_cur_div_absolute <= DIVIDEND_TREATED_AS_ZERO) //cur_div_percentage != 0
                {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= n_steps; ++i)
                    {
                        if (only_percentage_divs)
                        {
                            prices_odd[i] /= cur_div_mul_percentage;
                            prices_even[i] /= cur_div_mul_percentage;
                        }
                        option_vals[i] = tbmax(prices[i] - (X + cur_div_percentage * prices[i]), 0.0 );
                    }
                }
                else
                {
                    double d_div_p = cur_div_absolute / cur_div_percentage;
                    for (i = 0; prices[i] < d_div_p && i <= n_steps; ++i)
                    {
                        option_vals[i] = tbmax(prices[i] - (X + cur_div_percentage * prices[i]), 0.0 );
                    }

#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for (; i <= n_steps; ++i)
                    {
                        option_vals[i] = tbmax(prices[i] - (X + cur_div_absolute), 0.0 );
                    }
                }
            }
            else
            {
                //Intrinsic value exceeds right-shifted payoff functions for American calls if cur_div_absolute>0.
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for(i = 0; i <= n_steps; ++i) {
                    option_vals[i] = tbmax(prices[i] - X, 0.0 );
                }

                if (abs_cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || abs_cur_div_percentage > DIVIDEND_TREATED_AS_ZERO) {
                    if (i_curr_div == 0 && all_divs_counted
                        && (cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || cur_div_percentage > DIVIDEND_TREATED_AS_ZERO) )
                    {
                        result.exercise_boundary_div_reached = true;
                    }
                }
            }
        }
        else
        {
            //Final condition is right-shifted by cur_div_absolute if cur_div_absolute>0.
            //It works so for European and American put options, because the shifted curve exceeds the intrinsic value.
            if (abs_cur_div_percentage <= DIVIDEND_TREATED_AS_ZERO) //cur_div_absolute != 0
            {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for(i = 0; i <= n_steps; ++i) {
                    if (prices[i] <= cur_div_absolute) {
                        option_vals[i] = X;
                    }
                    else {
                        option_vals[i] = tbmax( X + cur_div_absolute - prices[i] , 0.0 );
                    }
                }
            }
            else if (abs_cur_div_absolute <= DIVIDEND_TREATED_AS_ZERO) //cur_div_percentage != 0
            {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for(i = 0; i <= n_steps; ++i)
                {
                    if (only_percentage_divs) {
                        prices_odd[i] /= cur_div_mul_percentage;
                        prices_even[i] /= cur_div_mul_percentage;
                    }
                    option_vals[i] = tbmax( X + cur_div_percentage * prices[i] - prices[i] , 0.0 );
                }
            }
            else
            {
                double d_div_p = cur_div_absolute / cur_div_percentage;
                for (i = 0; prices[i] < d_div_p && i <= n_steps; ++i)
                {
                    option_vals[i] = tbmax( X + cur_div_percentage * prices[i] - prices[i] , 0.0 );
                }

#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for (; i <= n_steps; ++i)
                {
                    option_vals[i] = tbmax( X + cur_div_absolute - prices[i] , 0.0 );
                }
            }

            if (abs_cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || abs_cur_div_percentage > DIVIDEND_TREATED_AS_ZERO) {
                if (i_curr_div == 0 && all_divs_counted
                    && (cur_div_absolute < -DIVIDEND_TREATED_AS_ZERO || cur_div_percentage < -DIVIDEND_TREATED_AS_ZERO) )
                {
                    result.exercise_boundary_div_reached = true;
                }
            }
        }
    }
    else
    {
        prices = prices_odd;
        int res = 0;

        double r_const = r_curve.GetAverageValue(T-dt, T);
        double q_const = q_curve.GetAverageValue(T-dt, T);
        double *input_prices = NULL;
        if (abs_cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || abs_cur_div_percentage > DIVIDEND_TREATED_AS_ZERO)
        {
            if (abs_cur_div_percentage <= DIVIDEND_TREATED_AS_ZERO) //cur_div_absolute != 0
            {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for (i = 0; i <= n_steps; ++i) {
                    option_values_interpolated[i] = tbmax(prices[i] - cur_div_absolute, 1e-100);
                }
            }
            else if (abs_cur_div_absolute <= DIVIDEND_TREATED_AS_ZERO) //cur_div_percentage != 0
            {
                if (only_percentage_divs)
                {
                    for (i = 0; i <= n_steps; ++i)
                    {
                        option_values_interpolated[i] = prices[i];
                        prices_odd[i] /= cur_div_mul_percentage;
                        prices_even[i] /= cur_div_mul_percentage;
                    }
                }
                else
                {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for (i = 0; i <= n_steps; ++i) {
                        option_values_interpolated[i] = tbmax(prices[i] - cur_div_percentage * prices[i], 1e-100);
                    }
                }
            }
            else {
                double d_div_p = cur_div_absolute / cur_div_percentage;
                for (i = 0; prices[i] < d_div_p && i <= n_steps; ++i) {
                    option_values_interpolated[i] = tbmax(prices[i] - cur_div_percentage * prices[i], 1e-100);
                }

#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for (; i <= n_steps; ++i) {
                    option_values_interpolated[i] = tbmax(prices[i] - cur_div_absolute, 1e-100);
                }
            }
            input_prices = option_values_interpolated;

            if (i_curr_div == 0 && all_divs_counted) {
                if (call_put_flag == CP_CALL
                    && (cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || cur_div_percentage > DIVIDEND_TREATED_AS_ZERO))
                {
                    result.exercise_boundary_div_reached = true;
                } else if (call_put_flag == CP_PUT
                           && (cur_div_absolute < -DIVIDEND_TREATED_AS_ZERO || cur_div_percentage < -DIVIDEND_TREATED_AS_ZERO)) {
                    result.exercise_boundary_div_reached = true;
                }
            }
        }
        else {
            input_prices = prices;
        }

        res = BlackScholesFairsOnPriceGrid(call_put_flag, payment_style, X, r_const, q_const, dt, sigma, n_steps, input_prices, option_vals, err_str);
        if (res != EC_SUCCESS) {
            return res;
        }

        if (exercise_style == AMERICAN_STYLE)
        {
            if (abs_cur_div_absolute <= DIVIDEND_TREATED_AS_ZERO)
            {
                double payoff;
                if (call_put_flag == CP_CALL)
                {
                    for (i = n_steps - 1; i >= 0 ; --i)
                    {
                        payoff = prices[i] - X;
                        if (payoff > option_vals[i]) {
                            option_vals[i] = payoff;
                        }
                        else {
                            break;
                        }
                    }
                }
                else
                {
                    for (i = 0; i < n_steps; ++i)
                    {
                        payoff = X - prices[i];
                        if (payoff > option_vals[i]) {
                            option_vals[i] = payoff;
                        }
                        else {
                            break;
                        }
                    }
                }
            }
            else
            {
                if (call_put_flag == CP_CALL)
                {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for (i = 0; i < n_steps; ++i) {
                        option_vals[i] = tbmax(prices[i] - X, option_vals[i]);
                    }
                }
                else
                {
                    for (i = 0; i < n_steps; ++i)
                    {
                        double X_ = X;

                        if (r_const < 0 && payment_style == PAYMENT_STYLE_EQUITY) {
                            X_ *= std::exp(-r_const * dt);
                        }

                        if (prices[i] < cur_div_absolute) {
                            option_vals[i] = X_;
                        }
                        else {
                            option_vals[i] = tbmax(X - prices[i], option_vals[i]);
                        }
                    }
                }
            }

        }
    }

    double t_prev_step = t - dt * start_roll_down - dt * 1e-3;
    int steps_done = use_smoothing ? 1 : 0;
    double t_left_node_r = 1e100;
    double t_left_node_q = 1e100;
    double r_dt = 0;
    double q_dt = 0;
    double growth_factor = 1.0;
    double p_up = 1.0;
    double p_down = 0;
    double pv = 1.0;
    double p_up_pv = 1.0;
    double p_down_pv = 0.0;

    double r_dt_next = 0;
    double q_dt_next = 0;
    bool r_prev_is_correct = false;
    bool q_prev_is_correct = false;
    double t_left_node_r_next = 1e100;
    double t_left_node_q_next = 1e100;
    double is_set_r_next = false;
    double is_set_q_next = false;
    double t_next_step = 0;
    double t_cur_step = 0;
    for(step = n_steps - start_roll_down; step >= add_steps; --step, t_prev_step -= dt, ++steps_done){
        double t_cur = dt * step;
        t_next_step = t_cur - add_t;
        t_cur_step = t_next_step + dt;

        if (t_next_step < t_left_node_r)
        {
            if (r_prev_is_correct)
            {
                if (t_left_node_r < t_cur_step)
                {
                    r_curve.GetTimeLeftNode(t_left_node_r, t_left_node_r_next);

                    if (t_next_step - dt < t_left_node_r_next)
                    {
                        r_prev_is_correct = false;
                        r_dt = r_curve.GetTimeIntegral(t_next_step, t_cur_step);
                    }
                    else
                    {
                        r_dt_next = r_curve.GetTimeIntegral(t_next_step - dt, t_next_step);
                        r_dt += (t_next_step - t_left_node_r) * (r_dt - r_dt_next) / dt;
                        is_set_r_next = true;
                    }
                }
                else
                {
                    if (is_set_r_next)
                    {
                        r_dt = r_dt_next;
                        is_set_r_next = false;
                    }
                    else
                    {
                        r_dt = r_curve.GetTimeIntegral(t_next_step, t_cur_step);
                        r_prev_is_correct = false;
                    }
                }
            }
            else
            {
                r_dt = r_curve.GetTimeIntegral(t_next_step, t_cur_step);
                r_curve.GetTimeLeftNode(t_cur_step, t_left_node_r_next);
                if (!(t_next_step < t_left_node_r_next && t_left_node_r_next < t_cur_step))
                {
                    r_prev_is_correct = true;
                }
            }

            growth_factor = exp(r_dt - q_dt);
            p_up = (growth_factor - d) / (u - d);
            p_down = 1.0 - p_up;

            if (payment_style == PAYMENT_STYLE_EQUITY) {
                pv = exp(-r_dt);
            }

            p_up_pv = p_up * pv;
            p_down_pv = p_down * pv;

            if (!q_is_zero && (growth_factor <= d)) {
                err_str = "Volatility too small for binomial model (q > r)";
                return EC_BAD_bin_Volatility_too_small;
            }

            if (p_up > 0.9999999999) {
                err_str = "Volatility too small for binomial model. Probability of up-move >= 1.";
                return EC_BAD_bin_Volatility_too_small;
            }

            if (t_cur_step < t_left_node_r) {
                if (t_left_node_r_next < t_left_node_r)
                {
                    t_left_node_r = t_left_node_r_next;
                }
                else
                {
                    r_curve.GetTimeLeftNode(t_cur_step, t_left_node_r);
                }
            }
        }

        if (t_next_step < t_left_node_q) {

            if (q_prev_is_correct)
            {
                if (t_left_node_q < t_cur_step)
                {
                    q_curve.GetTimeLeftNode(t_left_node_q, t_left_node_q_next);
                    if (t_next_step - dt < t_left_node_q_next)
                    {
                        q_prev_is_correct = false;
                        q_dt = q_curve.GetTimeIntegral(t_next_step, t_cur_step);
                    }
                    else
                    {
                        q_dt_next = q_curve.GetTimeIntegral(t_next_step - dt, t_next_step);
                        q_dt += (t_next_step - t_left_node_q) * (q_dt - q_dt_next) / dt;
                        is_set_q_next = true;
                    }
                }
                else
                {
                    if (is_set_q_next)
                    {
                        q_dt = q_dt_next;
                        is_set_q_next = false;
                    }
                    else
                    {
                        q_dt = q_curve.GetTimeIntegral(t_next_step, t_cur_step);
                        q_prev_is_correct = false;
                    }
                }
            }
            else
            {
                q_dt = q_curve.GetTimeIntegral(t_next_step, t_cur_step);
                q_curve.GetTimeLeftNode(t_cur_step, t_left_node_q_next);
                if (!(t_next_step < t_left_node_q_next && t_left_node_q_next < t_cur_step))
                {
                    q_prev_is_correct = true;
                }
            }

            growth_factor = exp(r_dt - q_dt);
            p_up = (growth_factor - d) / (u - d);
            p_down = 1.0 - p_up;
            p_up_pv = p_up * pv;
            p_down_pv = p_down * pv;

            if (!q_is_zero && (growth_factor <= d)) {
                err_str = "Volatility too small for binomial model (q > r)";
                return EC_BAD_bin_Volatility_too_small;
            }

            if (p_up > 0.9999999999) {
                err_str = "Volatility too small for binomial model. Probability of up-move >= 1.";
                return EC_BAD_bin_Volatility_too_small;
            }

            if (t_cur_step < t_left_node_q) {
                if (t_left_node_q_next < t_left_node_q)
                {
                    t_left_node_q = t_left_node_q_next;
                }
                else
                {
                    q_curve.GetTimeLeftNode(t_cur_step, t_left_node_q);
                }
            }
        }

        cur_div_absolute = 0.0;
        cur_div_percentage = 0.0;
        cur_div_mul_percentage = 1.0;
        for (;(!all_divs_counted) && (i_curr_div >= 0) &&
              (dividends.GetTime(i_curr_div)+add_t > t_prev_step);)
        {
            cur_div_absolute += dividends.GetAbsolute(i_curr_div) *
                                exp(-r_dt * (dividends.GetTime(i_curr_div) - t_next_step) / dt);

            cur_div_percentage += dividends.GetPercentage(i_curr_div);

            if (only_percentage_divs)
            {
                cur_div_mul_percentage *= 1 - dividends.GetPercentage(i_curr_div);
            }

            if (i_curr_div > 0) {
                i_curr_div--;
            }
            else {
                all_divs_counted = true;
            }
        }

        abs_cur_div_absolute = std::abs(cur_div_absolute);
        abs_cur_div_percentage = std::abs(cur_div_percentage);

        if (abs_cur_div_absolute > DIVIDEND_TREATED_AS_ZERO || abs_cur_div_percentage > DIVIDEND_TREATED_AS_ZERO)
        {
            int k = 1;
            double D = 0.0;

            if (only_percentage_divs)
            {
                double *prices_even_head = &prices_even[(steps_done + 1) / 2];
                double *prices_odd_head  = &prices_odd[steps_done / 2];
                double inv_cur_div_mul_percentage = 1.0 / cur_div_mul_percentage;
                for(i = 0; i <= step; ++i)
                {
                    prices_even_head[i] *= inv_cur_div_mul_percentage;
                    prices_odd_head[i]  *= inv_cur_div_mul_percentage;
                    option_values_interpolated[i] = p_up * option_vals[i + 1] + p_down * option_vals[i];
                }
            }
            else
            {
                for(i = 0; i <= step; ++i)
                {
                    if (abs_cur_div_percentage <= DIVIDEND_TREATED_AS_ZERO) { //cur_div_absolute != 0
                        D = cur_div_absolute;
                    }
                    else if (abs_cur_div_absolute <= DIVIDEND_TREATED_AS_ZERO) { //cur_div_percentage != 0
                        D = cur_div_percentage * prices[i];
                    }
                    else {
                        if (cur_div_absolute > DIVIDEND_TREATED_AS_ZERO && cur_div_percentage > DIVIDEND_TREATED_AS_ZERO) {
                            D = tbmin(cur_div_percentage * prices[i], cur_div_absolute);
                        }
                        else if (cur_div_absolute < (-DIVIDEND_TREATED_AS_ZERO) && cur_div_percentage < (-DIVIDEND_TREATED_AS_ZERO)) {
                            D = tbmax(cur_div_percentage * prices[i], cur_div_absolute);
                        }
                    }

                    double S_disc = tbmax(prices[i] * u - D, 0.0);
                    while ((S_disc >= prices[k]) && (k <= step)) {
                        k++;
                    }

                    option_values_interpolated[i] = (option_vals[k] - option_vals[k - 1]) /
                                                    (prices[k] - prices[k - 1]) *
                                                    (S_disc - prices[k - 1]) +
                                                    option_vals[k - 1];

                    option_values_interpolated[i] = tbmax(option_values_interpolated[i], 0.0);
                }
            }

            // done with interpolation, now we can move down in the tree
            if (steps_done % 2 == 1) {
                prices = &prices_even[(steps_done + 1) / 2];
            }
            else {
                prices = &prices_odd[steps_done / 2];
            }

            if (exercise_style == STYLE_AMERICAN)
            {
                // exercise boundary analysis on the last step
                if(step == add_steps && result.get_exercise_boundary)
                {
                    FindExerciseBoundary(call_put_flag, X, step, &prices[0], &option_values_interpolated[0], pv, result.exercise_boundary, result.exercise_boundary_state);
                }

                // Analyze "exercise boundary div reached" before the application of intrinsic value correction
                // since dividend is clinged to the next time layer
                if (all_divs_counted) { // Consider the first dividend
                    if(result.get_exercise_boundary_div)
                    {
                        FindExerciseBoundary(call_put_flag, X, step, &prices[0], &option_values_interpolated[0], pv, result.exercise_boundary_div, result.exercise_boundary_div_state);
                    }

                    double cur_r_minus_q = r_curve.GetTimeIntegral(0, t_next_step) - q_curve.GetTimeIntegral(0, t_next_step);
                    /*
                     * Since u = exp(sigma * sqrt(dt)),
                     * S_n = S*u^n = S*d^{-n}, where n is the relative index of price step from the initial spot price S,
                     * A distance between steps within temporal layer is S*u*u - S = S(u^2 - 1)
                     * we calculate forward_S = S*exp((r-q)*t_div), equate them S_n = S_forward to find the node where forward_S lays
                     * S*u^(n*2) = S*exp((r-q)*dt*step) => sigma * sqty(dt) * n * 2 = (r-q)*dt*n =>
                     * n = 0.5 * (r-q)*t_div / (sigma * sqty(dt))
                     *
                     * step / 2 is the index of S if step % 2 == 0 or the index of left node that is closer to S if step % 2 == 1.
                     * Since we do not have S in odd price grid we add 1 additional step
                     */
                    size_t step_odd = step % 2;
                    int index_left = step / 2 + static_cast<int>((step_odd + cur_r_minus_q / (sigma * sqrt(dt))) * 0.5);
                    int index_right = index_left + 1;

                    if (index_left >= 0 && index_right <= step) {
                        double forward_S = S * exp(cur_r_minus_q);
                        double payoff = ((call_put_flag == CP_CALL)
                                         ? forward_S - X : X - forward_S);

                        double option_val_at_div = pv * LInterpTwoPoints(prices[index_left], prices[index_right], forward_S,
                                                                         option_values_interpolated[index_left],
                                                                         option_values_interpolated[index_right]);

                        result.exercise_boundary_div_reached = option_val_at_div < trembling * payoff;

                        result.forward_S = forward_S;
                        result.fair_at_next_dividend = option_val_at_div;
                        result.delta_at_next_dividend =  (option_values_interpolated[index_right] - option_values_interpolated[index_left]) / (prices[index_right] - prices[index_left]);
                    }
                }

                if (call_put_flag == CP_CALL) {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= step; ++i) {
                        option_vals[i] = pv * option_values_interpolated[i];
                        option_vals[i] = tbmax(option_vals[i], prices[i] - X);
                    }
                } else {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= step; ++i) {
                        option_vals[i] = pv * option_values_interpolated[i];
                        option_vals[i] = tbmax(option_vals[i], X - prices[i]);
                    }
                }
            }
            else
            {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                for (i = 0; i <= step; ++i) {
                    option_vals[i] = pv * option_values_interpolated[i];
                }
            }
        }
        else
        {
            //option_values_interpolated is temporary array here
            //operations are splitted for auto-vectorization
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
            for (i = 0; i <= step; ++i) {
                option_values_interpolated[i] =  p_up_pv * option_vals[i + 1];
            }

#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
            for (i = 0; i <= step; ++i) {
                option_vals[i] = p_down_pv * option_vals[i] + option_values_interpolated[i];
            }

            if (steps_done % 2 == 1) {
                prices = &prices_even[(steps_done + 1) / 2];
            }
            else {
                prices = &prices_odd[steps_done / 2];
            }

            if(exercise_style == STYLE_AMERICAN)
            {
                if(step == add_steps && result.get_exercise_boundary)
                {
                    FindExerciseBoundary(call_put_flag, X, step, &prices[0], &option_vals[0], /* pv */1.0, result.exercise_boundary, result.exercise_boundary_state);
                }

                //double payoff;
                if (call_put_flag == CP_CALL)
                {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= step; ++i) {
                        option_vals[i] = tbmax(option_vals[i], prices[i] - X);
                    }

                    /*
                    // This variant leads to performance degradation due to backward loop direction
                    for(i = step; i >= 0 && (payoff = (prices[i] - X)) > 0 ; --i) {
                        option_vals[i] = tbmax(option_vals[i], payoff);
                    }
                    */
                } else {
#ifdef __clang__
#pragma clang loop vectorize_width(4) //will work with clang 3.5 or higher
#endif
                    for(i = 0; i <= step; ++i) {
                        /*
                        payoff = X - prices[i];
                        if (payoff <= 0) {
                            break;
                        }
                        */
                        option_vals[i] = tbmax(option_vals[i], X - prices[i]);
                    }
                }
            }
        }
    }

    size_t index_center = add_steps / 2;
    size_t index_left = index_center - 1;
    size_t index_right = index_center + 1;

    result.fair_value = option_vals[index_center];

    double intrinsic_value = (call_put_flag == CP_CALL
                              ? tbmax(S - X, 0.0)
                              : tbmax(X - S, 0.0));

    if ( std::abs(result.fair_value) < trembling * intrinsic_value)
    {
        result.delta = (call_put_flag == CP_CALL ? 1.0 : -1.0);
        result.gamma = 0.0;
    }
    else
    {
        // Calculating Delta and Gamma with the account of non-equidistant grid
        double dy_r = option_vals[index_right] - option_vals[index_center];
        double dy_l = option_vals[index_center] - option_vals[index_left];
        double dx_r = prices[index_right] - prices[index_center];
        double dx_l = prices[index_center] - prices[index_left];
        double dx_sum = dx_l + dx_r;

        result.delta = dy_r * dx_l / dx_r / dx_sum;
        result.delta += dy_l * dx_r / dx_l / dx_sum;
        result.gamma = (dy_r / dx_r - dy_l / dx_l) * 2.0 / dx_sum;
    }

#ifdef CONTINUOUS_TIME_DERIVATIVES
    result.theta = (result.theta - option_vals[price_index]) * 0.5 / ( t / steps);
#endif

    return EC_SUCCESS;
}
