## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================

############################################################################
#####
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")

# add Scotch library?
SET ( USE_SCOTCH "" CACHE STRING "Use SCOTCH tool for renumbering (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_SCOTCH PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_SCOTCH OR USE_SCOTCH STREQUAL "" )
  # Variable is not provided by user
  FIND_PACKAGE(SCOTCH QUIET)

ELSE ()
  IF ( USE_SCOTCH )
    # User wants to use scotch
    FIND_PACKAGE(SCOTCH)
    IF ( NOT SCOTCH_FOUND )
      MESSAGE ( FATAL_ERROR "Scotch library not found:"
        "Using scotch reduce the execution time of mmg3d "
        "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
        "If you have already installed Scotch and want to use it, "
        "please set the CMake variable or environment variable SCOTCH_DIR "
        "to your scotch directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

If ( SCOTCH_FOUND )
  add_definitions(-DUSE_SCOTCH)

  MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARIES}")
  SET( LIBRARIES ${SCOTCH_LIBRARIES} ${LIBRARIES})
ENDIF()


############################################################################
#####
#####         LinearElasticity
#####
############################################################################
# add LinearElasticity library?
SET(ELAS_DIR "" CACHE PATH "Installation directory for Elas")

SET ( USE_ELAS "" CACHE STRING "Use the Elas library for lagrangian motion option (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_ELAS PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_ELAS OR USE_ELAS STREQUAL ""  )
  INCLUDE(FindElas)

ELSE()
  IF ( USE_ELAS )
    # User wants to use elas
    INCLUDE(FindElas)
    IF ( NOT ELAS_FOUND )
      MESSAGE ( FATAL_ERROR "Elas is a library to solve the linear elasticity "
        "problem (see https://github.com/ISCDtoolbox/LinearElasticity to"
        " download it). "
        "This library is needed to use the lagrangian motion option. "
        "If you have already installed Elas and want to use it, "
        "please set the CMake variable or environment variable ELAS_DIR "
        "to your Elas directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

############################################################################
#####
#####         VTK (to parse (p)vtp/(p)vtu files )
#####
############################################################################
# add the VTK library ?
SET ( USE_VTK "" CACHE STRING "Use VTK I/O (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_VTK PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_VTK OR USE_VTK STREQUAL "" OR USE_VTK )
  # USE_VTK is not false, ie that it is true or empty (or contains a fake value)

  # Handle vtk components name change between v8.2 and v9
  # Before v9
  FIND_PACKAGE(VTK QUIET)
  IF ( VTK_FOUND )
    message (STATUS "VTK_VERSION: ${VTK_VERSION}")
    IF (VTK_VERSION VERSION_LESS "9.0.0")
      find_package(VTK  COMPONENTS
        vtkCommonCore
        vtkCommonDataModel
        vtkIOLegacy
        vtkIOParallel
        vtkIOParallelXML
        vtkIOXML
        QUIET)
    ELSE()
      # After v9
      FIND_PACKAGE(VTK COMPONENTS
        CommonCore
        CommonDataModel
        IOLegacy
        IOParallel
        IOParallelXML
        IOXML
        QUIET)
    ENDIF()

  ELSEIF ( USE_VTK )
    # USE_VTK is not empty so user explicitely ask for VTK: raise an error
    MESSAGE(FATAL_ERROR "VTK library not found.")
  ENDIF()
ENDIF()

IF ( VTK_FOUND )
  add_definitions(-DUSE_VTK)

  MESSAGE ( STATUS "Compilation with VTK: add vtk, vtp and vtu I/O." )

  IF( "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" LESS 8.90 )
    INCLUDE ( ${VTK_USE_FILE} )
  ENDIF()

  SET( LIBRARIES ${VTK_LIBRARIES} ${LIBRARIES} )
ENDIF ( )

############################################################################
#####
#####         StarPU
#####
############################################################################
# Find STARPU library?
SET(STARPU_DIR "" CACHE PATH "Installation directory for StarPU")

SET ( USE_STARPU "" CACHE STRING "Use the StarPU scheduler for shared memory parallelization (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_STARPU PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_STARPU OR USE_STARPU STREQUAL ""  )
  # Variable is not provided by user
  FIND_PACKAGE(STARPU QUIET)

ELSE ()

  IF ( USE_STARPU )
    # User wants to use starPU
    FIND_PACKAGE(STARPU)
    IF ( NOT STARPU_FOUND )
      MESSAGE ( FATAL_ERROR "StarPU library not found:"
      "If you have already installed StarPU and want to use it, "
      "please set the CMake variable or environment variable STARPU_DIR "
      "to your StarPU directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

If ( STARPU_FOUND )
  add_definitions(-DUSE_STARPU)

  MESSAGE(STATUS
    "Compilation with starPU: ${STARPU_LIBRARY_DIRS}/${STARPU_LIBRARIES}")
  SET( LIBRARIES ${STARPU_LIBRARIES} ${LIBRARIES})

  # FindSTARPU returns library names without paths: add paths
  link_directories(AFTER ${STARPU_LIBRARY_DIRS})

ENDIF()


############################################################################
#####
#####         FXT
#####
############################################################################

## warning Useless I think: FindSTARPU already search for FXT

# Find FXT library for StarPU execution trace?
SET(FXT_DIR "" CACHE PATH "Installation directory for FXT")

SET ( USE_FXT "" CACHE STRING "Use the FXT trace tool for StarPU (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_FXT PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_FXT OR USE_FXT STREQUAL ""  )
  # Variable is not provided by user
  FIND_PACKAGE(FXT QUIET)

ELSE ()

  IF ( USE_FXT )
    # User wants to use FXT
    FIND_PACKAGE(FXT)
    IF ( NOT FXT_FOUND )
      MESSAGE ( FATAL_ERROR "FXT library not found:"
      "If you have already installed FXT and want to use it, "
      "please set the CMake variable or environment variable FXT_DIR "
      "to your FXT directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

If ( FXT_FOUND )
  add_definitions(-DUSE_FXT)

  MESSAGE(STATUS
    "Compilation with FXT: ${FXT_LIBRARIES}")
  SET( LIBRARIES ${FXT_LIBRARIES} ${LIBRARIES})
ENDIF()


############################################################################
#####
#####        Metis
#####
############################################################################

# Find METIS library?
SET(METIS_DIR "" CACHE PATH "Installation directory for METIS")
SET ( USE_METIS "" CACHE STRING "Use the Metis graph partitionner (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_METIS PROPERTY STRINGS "ON" "OFF" " ")

IF ( NOT DEFINED USE_METIS OR USE_METIS STREQUAL ""  )
  # Variable is not provided by user
  FIND_PACKAGE(METIS QUIET)

ELSE ()
  IF ( USE_METIS )
    # User wants to use METIS
    FIND_PACKAGE(METIS)
    IF ( NOT METIS_FOUND )
      MESSAGE ( FATAL_ERROR "METIS library not found:"
      "If you have already installed METIS and want to use it, "
      "please set the CMake variable or environment variable METIS_DIR "
      "to your METIS directory.")
    ENDIF ( )
  ENDIF ( )

ENDIF ( )

If ( METIS_FOUND )
  add_definitions(-DUSE_METIS)

  MESSAGE(STATUS
    "Compilation with METIS: ${METIS_LIBRARIES}")
  SET( LIBRARIES ${METIS_LIBRARIES} ${LIBRARIES})
ENDIF()
