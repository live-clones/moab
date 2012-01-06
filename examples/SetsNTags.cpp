#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include <iostream>

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }

    // get the material set tag handle
  moab::Tag mtag;
  moab::ErrorCode rval;
  const char *tag_nms[] = {"MATERIAL_SET", "DIRICHLET_SET", 
                           "NEUMANN_SET"};
  moab::Range sets, set_ents;

    // instantiate & load a file
  moab::Interface *mb = new moab::Core();
  rval = mb->load_file(argv[1]);

    // loop over set types
  for (int i = 0; i < 3; i++) {
    rval = mb->tag_get_handle(tag_nms[i], 1, moab::MB_TYPE_INTEGER, mtag);

      // get all the sets of that type in the mesh
    sets.clear();
    rval = mb->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &mtag,
                                            NULL, 1, sets);

      // iterate over each set, getting entities
    moab::Range::iterator set_it;
    for (set_it = sets.begin(); set_it != sets.end(); set_it++)  {
      moab::EntityHandle this_set = *set_it;

        // get the id for this set
      int set_id;
      rval = mb->tag_get_data(mtag, &this_set, 1, &set_id);

        // get the entities in the set, recursively
      rval = mb->get_entities_by_handle(this_set, set_ents, true);

      std::cout << tag_nms[i] << " " << set_id << " has " 
                << set_ents.size() << " entities." << std::endl;
      set_ents.clear();
    }
  }
}