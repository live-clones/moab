#define TARGET_TEST_GROUP "ShapeSizeOrientTargetMetricTests"
#include "TargetMetricTest.hpp"

using namespace MBMesquite;

#include "AWShapeSizeOrientNB1.hpp"
#include "TShapeSizeOrientB1.hpp"
#include "TShapeSizeOrientB2.hpp"
#include "TShapeSizeOrientNB1.hpp"

//                            NAME     !SHAPE !SIZE !ORIENT BARRIER
TEST_METRIC_WITH_HESS( AWShapeSizeOrientNB1, false, false, false, false, 0.0 );
TEST_METRIC_WITH_HESS( TShapeSizeOrientB1, false, false, false, true, 0.0 );
TEST_METRIC_WITH_HESS( TShapeSizeOrientB2, false, false, false, true, 0.0 );
TEST_METRIC_WITH_HESS( TShapeSizeOrientNB1, false, false, false, false, 0.0 );
