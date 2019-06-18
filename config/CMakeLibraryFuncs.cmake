# Fix up a library, create appropriate defines for shared, static
#
# Proj:    The name of the project
# LibName: The name of the library with the capitalization it will have
macro(SciFixupLibrary Proj LibName DeclName)

  string(TOLOWER "${LibName}" libname)
  string(TOUPPER "${LibName}" LIBNAME)
  string(TOUPPER "${DeclName}" DECLNAME)

# Set the version
  if (BUILD_SHARED_LIBS)
    if (VERSION_MAJOR)
      set_target_properties(${LibName} PROPERTIES
          SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}
          VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
      )
    endif ()
    message(STATUS "LibName = ${LibName}.")
    target_compile_definitions(${LibName}
      PRIVATE ${DeclName}_EXPORTS
      # ${Proj}_DLL # Added at the top
    )
  else ()
    target_compile_definitions(${LibName} PRIVATE ${DECLNAME}_STATIC_DEFINE)
  endif ()

endmacro()

