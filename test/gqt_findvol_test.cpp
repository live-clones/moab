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

  int result = 0;

  // with no global OBB tree (will defer to find_volume_slow)
  result += RUN_TEST(find_volume_tests);

  // build OBBs with one-vol tree
  GTT->construct_obb_trees(true);

  // using the global OBB tree
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

// Test Geometry \\

// The test geometry consists of 3 cubes, one of which
// overlaps another. Volumes 1 and 2 have edges with
// length 1. Volume 3 has an edge length of 0.5.
// Volume 1 is centered at (3, 0, 0). Volume 2 is
// centered on the origin. Volume 3 is centered on (0.5, 0, 0).
// Volume 4 is the implicit complement.

// XY slice of geometry:

// ###########################                   #############################
// #                         #                   #                           #
// #                         #                   #                           #
// #                ####################         #                           #
// #                #        #         #         #                           #
// #      Vol 2     #        #  Vol 3  #         #           Vol 1           #
// #                #        #         #         #                           #
// #                #        #         #         #                           #
// #                ####################         #                           #
// #                         #                   #                           #
// #                         #     Vol 4 (IC)    #                           #
// ###########################                   #############################

void find_volume_tests() {

  ErrorCode rval;

  const struct FindVolTestResult tests[] = {
    // one point unambiguously placed in each volume
    // and the implicit complement
    { {-0.1, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 2, -1},   // 1
    { { 0.6, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 3, -1},   // 2
    { { 3.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 1, -1},   // 3
    { {-5.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 4,  0},   // 4
    // Point on the negative side of the geometry
    { {-5.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, 0, -1},   // 5
    { {-5.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 0, -1},   // 6
    // Point on the positive side of the geometry
    { {10.0, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 1,  0},   // 7
    { {10.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 0, -1},   // 8
    { {10.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, 0, -1},   // 9
    // Point between the volumes
    { { 1.5, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 0,  4},   // 10
    { { 1.5, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 4, -1},   // 11
    { { 1.5, 0.0, 0.0}, { 1.0, 0.0, 0.0}, 4, -1},   // 12
    { { 1.5, 0.0, 0.0}, { 0.0, 1.0, 0.0}, 0, -1},   // 13
    // Point in the overlap of vols 2 & 3
    { { 0.4, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 2,  3},   // 14
    { { 0.4, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 2,  3},   // 15
    { { 0.4, 0.0, 0.0}, { 1.0, 0.0, 0.0}, 2,  3},   // 16
    // Point in Vol 3 w/ different directions applied
    { { 0.6, 0.0, 0.0}, { 0.0, 0.0, 0.0}, 3, -1},   // 17
    { { 0.6, 0.0, 0.0}, {-1.0, 0.0, 0.0}, 3, -1},   // 18
    { { 0.6, 0.0, 0.0}, { 1.0, 0.0, 0.0}, 3, -1} }; // 19

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

    std::cout << "Test " << i << ". Volume found id: " << vol_id << "\n";
    // make sure at least one of these checks passed
    CHECK(vol_id == test.resultA || vol_id == test.resultB);

    // reset result and id for safety
    volume_found = 0;
    vol_id = -1;
  }
}
