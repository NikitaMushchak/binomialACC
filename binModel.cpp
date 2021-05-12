//
// Created by nikita.mushak on 4/28/2021.
//
//#include <strategy/type/Double.h>

//#include "shared/pricing/models/impl/Binomial.h"
#include "shared/dividend_container/model/DividendContainer.h"
#include "shared/pricing/models/ModelResult.h"
#include "shared/interest_rate/YieldCurve.h"
#include "shared/interest_rate/YieldCurveWrapper.h"

#include <chrono>
#include <sstream>
#include <string>
#include <array>
#include <string_view>

//#include "binomialCut.h"
#include "binomialOPENACC.h"

int main(){
    auto callput       = CP_CALL;
    auto exercise =  STYLE_EUROPEAN;//STYLE_AMERICAN;
    auto payment  = PAYMENT_STYLE_EQUITY;
    auto smoothing      = false;

    const int numOfIterations = 1;

    const int numberOfOptions = 1;
    // Generate Options
//    std::vector<CallOption> optionsVector;
//    GenerateOptions(optionsVector, numberOfOptions);
    size_t iter_test;
    iter_test = 0;
    std::vector<double> timeOfIteration;
    for (int iter = 0; iter < numOfIterations; ++iter){
        timeOfIteration.clear();
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < numberOfOptions; ++i) {
//
//            struct ModelResult result;
//            result.get_fair_value = 1;
//            result.get_delta = 1;
//            result.get_gamma = 1;
//            result.get_vega = 1;
//            result.get_rho = 1;
//            result.get_vanna = 1;
//            result.get_vomma = 1;
//
//            double underlyingPrice = 110 + i * 0.25, strikePrice = 100 - i * 025, rT = 0.43013698630137 + i * 0.000002,
//                    sT = 0.43013698630137 + i * 0.0000020002,
//                    eT = 0.43013698630137 + i * 0.0000020001,
//                    sigma = 0.2 + i * 0.000000002;
//
//            YieldCurve r(tbricks::instrument_parameters::ZCYC());
//            YieldCurveWrapper r_curve(0.01, r, tbricks::Double()),
//                    q_curve(0.01, r, tbricks::Double());
//
//            PrecisionSetting precision = PrecisionSetting::VERY_PRECISE;
//            const char *err = nullptr;
//
//            dividend_container::DividendContainer divs;
//
//            BinomialWithDiscreteDivs(callput, exercise, payment, result, underlyingPrice, strikePrice,
//                                     r_curve, q_curve, rT, sT, eT, sigma, divs, precision, err, smoothing);
//            struct ModelResult result1;
//            result1.get_fair_value = 1;
//            result1.get_delta = 1;
//            result1.get_gamma = 1;
//            result1.get_vega = 1;
//            result1.get_rho = 1;
//            result1.get_vanna = 1;
//            result1.get_vomma = 1;
//
//            BinomialWithDiscreteDivs_Cut(callput, exercise, payment, result1,
//                                                underlyingPrice, strikePrice,  r_curve,
//                                              q_curve, rT, sigma, divs, smoothing, err);

            struct ModelResult result2;
            result2.get_fair_value = 1;
            result2.get_delta = 1;
            result2.get_gamma = 1;
            result2.get_vega = 1;
            result2.get_rho = 1;
            result2.get_vanna = 1;
            result2.get_vomma = 1;
            BinomialWithDiscreteDivs_ACC(callput, exercise, payment, result2,
                                         underlyingPrice, strikePrice,  r_curve,
                                         q_curve, rT, sigma, divs, smoothing, err);
//            state.SetLabel(label);
            iter_test++;

        }
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds =
                std::chrono::duration_cast < std::chrono::duration < double >> (
                        end - start);
        timeOfIteration.push_back(elapsed_seconds.count());

//        state.SetIterationTime(elapsed_seconds.count());
    }
    double totalTime = 0.;
    double averageTime = 0.;

    for(auto iterTime : timeOfIteration)
        totalTime +=iterTime;


    averageTime=totalTime/numOfIterations;

//    printf("result = %2.10f \n",result.fair_value);
    std::cout<<"number of iters = "<<iter_test<<'\n';
    std::cout<<"host aver of #  "<<iter_test<<" with num of iteration : "<<numOfIterations
             <<" time = "<< averageTime<<" sec"<<'\n';
    return 0;
}

