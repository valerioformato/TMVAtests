# - Locate XGBoost library
# Defines:
#
#  XGBoost_FOUND
#  XGBoost_INCLUDE_DIR
#  XGBoost_INCLUDE_DIRS (not cached)
#  XGBoost_LIBRARY
#  XGBoost_LIBRARIES (not cached)
#  XGBoost_LIBRARY_DIRS (not cached)
#  XGBoost_EXTRA_INCLUDE_DIR


find_path(XGBoost_INCLUDE_DIR xgboost/c_api.h
          HINTS $ENV{XGBoost_ROOT_DIR}/include ${XGBoost_ROOT_DIR}/include
)

find_path(XGBoost_EXTRA_INCLUDE_DIR rabit/c_api.h
          HINTS ${XGBoost_INCLUDE_DIR}/../rabit/include
          $ENV{XGBoost_ROOT_DIR}/rabit/include ${XGBoost_ROOT_DIR}/rabit/include
)

find_path(XGBoost_DMLC_INCLUDE_DIR dmlc/data.h
          HINTS ${XGBoost_INCLUDE_DIR}/../dmlc-core/include
          $ENV{XGBoost_ROOT_DIR}/dmlc-core/include ${XGBoost_ROOT_DIR}/dmlc-core/include
)

find_library(XGBoost_LIBRARY NAMES xgboost
          HINTS $ENV{XGBoost_ROOT_DIR}/lib ${XGBoost_ROOT_DIR}/lib
)

find_library(XGBoost_STATIC_LIBRARY NAMES libxgboost.a
          HINTS $ENV{XGBoost_ROOT_DIR}/lib ${XGBoost_ROOT_DIR}/lib
)

find_library(XGBoost_EXTRA_STATIC_LIBRARY NAMES librabit.a
          HINTS $ENV{XGBoost_ROOT_DIR}/rabit/lib ${XGBoost_ROOT_DIR}/rabit/lib
)

find_library(XGBoost_DMLC_STATIC_LIBRARY NAMES libdmlc.a
          HINTS $ENV{XGBoost_ROOT_DIR}/dmlc-core ${XGBoost_ROOT_DIR}/dmlc-core
)

# handle the QUIETLY and REQUIRED arguments and set XGBoost_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XGBoost DEFAULT_MSG XGBoost_INCLUDE_DIR XGBoost_LIBRARY  XGBoost_EXTRA_INCLUDE_DIR)

mark_as_advanced(XGBoost_FOUND XGBoost_INCLUDE_DIR XGBoost_LIBRARY XGBoost_EXTRA_INCLUDE_DIR)

set(XGBoost_INCLUDE_DIRS ${XGBoost_INCLUDE_DIR} ${XGBoost_EXTRA_INCLUDE_DIR} ${XGBoost_DMLC_INCLUDE_DIR})
set(XGBoost_LIBRARIES ${XGBoost_LIBRARY})
set(XGBoost_STATIC_LIBRARIES ${XGBoost_STATIC_LIBRARY} ${XGBoost_EXTRA_STATIC_LIBRARY} ${XGBoost_DMLC_STATIC_LIBRARY})
get_filename_component(XGBoost_LIBRARY_DIRS ${XGBoost_LIBRARY} PATH)
