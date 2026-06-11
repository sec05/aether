find_program(CLANG_FORMAT_EXE NAMES clang-format)
if(NOT CLANG_FORMAT_EXE)
  message(FATAL_ERROR "clang-format not found")
endif()

file(GLOB_RECURSE AETHER_SOURCE_FILES
  "${CMAKE_SOURCE_DIR}/src/*.cpp"
  "${CMAKE_SOURCE_DIR}/include/*.hpp"
  "${CMAKE_SOURCE_DIR}/tests/*.cpp"
)

if(AETHER_SOURCE_FILES)
  execute_process(
    COMMAND ${CLANG_FORMAT_EXE} --dry-run --Werror ${AETHER_SOURCE_FILES}
    RESULT_VARIABLE FORMAT_CHECK_RESULT
  )

  if(NOT FORMAT_CHECK_RESULT EQUAL 0)
    message(FATAL_ERROR "clang-format check failed")
  endif()
endif()
