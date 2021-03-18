
FetchContent_GetProperties(r-index)

add_executable(split_fa split_fa.cpp)
target_link_libraries(split_fa klib z)
