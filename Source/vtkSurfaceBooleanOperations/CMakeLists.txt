PROJECT(SurfaceBooleanOperations)


CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


OPTION(BUILD_PARAVIEW_PLUGIN "Build paraview plugin" OFF)


IF(BUILD_PARAVIEW_PLUGIN)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ELSE(BUILD_PARAVIEW_PLUGIN)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
ENDIF(BUILD_PARAVIEW_PLUGIN)


FIND_PATH(GTS_LIBRARY_DIR 
      NAMES libgts.a libgts.so 
      PATHS /usr/lib /sw/lib NO_DEFAULT_PATH)
FIND_PATH(GLIB_LIBRARY_DIR 
      NAMES libglib-2.0.a libglib-2.0.so libglib-2.0.dylib
      PATHS /usr/lib /sw/lib NO_DEFAULT_PATH)

FIND_PATH(GTS_INCLUDE_DIR 
      NAMES gts.h 
      PATHS /usr/include /sw/include
      DOC "Directory which contains gts.h" NO_DEFAULT_PATH)
FIND_PATH(GLIB_INCLUDE_DIR 
      NAMES glib.h 
      PATHS /usr/include/glib-2.0 /sw/include/glib-2.0
      DOC "Directory which contains glib.h" NO_DEFAULT_PATH)
FIND_PATH(GLIBCONFIG_INCLUDE_DIR 
      NAMES glibconfig.h 
      PATHS /usr/lib/glib-2.0/include /sw/lib/glib-2.0/include 
      DOC "Directory which contains glibconfig.h" NO_DEFAULT_PATH)



SET(GTS_LIBS gts glib-2.0)
LINK_DIRECTORIES(${GTS_LIBRARY_DIR} ${GLIB_LIBRARY_DIR})
INCLUDE_DIRECTORIES(
  ${GTS_INCLUDE_DIR}
  ${GLIB_INCLUDE_DIR}
  ${GLIBCONFIG_INCLUDE_DIR}
)



ADD_EXECUTABLE(example-usage example-usage.cxx vtkSurfaceBooleanOperations.cxx)
TARGET_LINK_LIBRARIES(example-usage vtkGraphics vtkIO vtkCommon ${GTS_LIBS})

IF(BUILD_PARAVIEW_PLUGIN)

  ADD_PARAVIEW_PLUGIN(vtkSurfaceBooleanOperations "1.0" 

SERVER_MANAGER_XML vtkSurfaceBooleanOperations_SM.xml 

SERVER_MANAGER_SOURCES vtkSurfaceBooleanOperations.cxx

GUI_RESOURCE_FILES vtkSurfaceBooleanOperations_GUI.xml)

  TARGET_LINK_LIBRARIES(vtkSurfaceBooleanOperations ${GTS_LIBS})

ENDIF(BUILD_PARAVIEW_PLUGIN)



