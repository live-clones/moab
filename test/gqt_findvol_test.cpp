#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"

#include "moab/GeomQueryTool.hpp"
#include "moab/GeomTopoTool.hpp"

using namespace moab;

Core* MBI;
GeomTopoTool* GTT;
GeomQueryTool* GQT;

Tag id_tag = 0;

const std::string input_file = TestDir + "/find_vol_test_geom.h5m";

void find_volume_tests();

int main() {

  MBI = new Core();
  ErrorCode rval = MBI->load_file(input_file.c_str());
  MB_CHK_SET_ERR(rval, "Failed to load test file");

  GTT = new GeomTopoTool(MBI);
  GQT = new GeomQueryTool(GTT);

  // initialize the rest of the GQT
  GQT->initialize();

  // build OBBs with one-vol tree
  GTT->construct_obb_trees(true);


  int result = 0;
  result += RUN_TEST(find_volume_tests);

  return result;
}

ErrorCode id_lookup(EntityHandle eh, int& id) {
  ErrorCode rval;
  if (!id_tag) {
    // can be replaced with MBI->globalId_tag() soon
    rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, id_tag);
    MB_CHK_ERR(rval);
  }

  rval = MBI->tag_get_data(id_tag, &eh, 1, (void*)&id);
  MB_CHK_SET_ERR(rval, "Failed to lookup volume id");

  return MB_SUCCESS;
}

// Test case structure
// It is possible in some cases
// to have 2 correct answers (overlaps)
// A result of zero means no volume should be found

struct FindVolTestResult {
  double pnt[3];
  double dir[3];
  int resultA;
  int resultB;
};

void find_volume_tests() {

  ErrorCode rval;

  const struct FindVolTestResult tests[] = {
    { {-0.1, 0.0, 0.0}, { 0.0, 0.0, 0.0,}, 2, -1},
    { { 0.7, 0.0, 0.0}, { 0.0, 0.0, 0.0,}, 3, -1},
    { { 3.0, 0.0, 0.0}, { 0.0, 0.0, 0.0,}, 1, -1},
    { {-5.0, 0.0, 0.0}, { 1.0, 0.0, 0.0,}, 4, -1} };

  int num_tests = sizeof(tests)/sizeof(FindVolTestResult);

  EntityHandle volume_found;
  int vol_id;
  for (int i = 0; i < num_tests; i++) {

    const FindVolTestResult& test = tests[i];

    const double* direction = NULL;
    if (test.dir[0] != 0.0 || test.dir[1] != 0.0 || test.dir[2] != 0.0) {
      direction = test.dir;
    }

    rval = GQT->find_volume(test.pnt,
                            volume_found,
                            direction);
    // if not found, we will check later
    if (rval != MB_ENTITY_NOT_FOUND) {
      MB_CHK_SET_ERR_CONT(rval, "Failed in find_volume");
    }

    rval = id_lookup(volume_found, vol_id);
    MB_CHK_SET_ERR_CONT(rval, "Failed in id lookup");

    // if the expected result is valid (non-negative),
    // make sure there's a match. Otherwise, set true.
    bool a = test.resultA >= 0 ? vol_id == test.resultA : true;
    bool b = test.resultB >= 0 ? vol_id == test.resultB : true;
    std::cout << "Volume found id: " << vol_id << "\n";
    // make sure both checks passed
    CHECK(a && b);

    // reset result and id for safety
    volume_found = 0;
    vol_id = -1;
  }
}
