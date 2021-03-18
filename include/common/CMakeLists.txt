
set(COMMON_SOURCES common.hpp)

add_library(common OBJECT ${COMMON_SOURCES})
target_link_libraries(common)