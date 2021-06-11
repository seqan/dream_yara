cmake_minimum_required (VERSION 3.8)

if (MODEL MATCHES ".*MemCheck.*")
    set (_VALGRIND_FLAG --valgrind)
else ()
    set (_VALGRIND_FLAG)
endif ()

find_package (PythonInterp 3)

if (PYTHONINTERP_FOUND)
  add_test (NAME app_test_dream_yara
            COMMAND ${PYTHON_EXECUTABLE}
                    ${CMAKE_SOURCE_DIR}/test/cli/run_tests.py
                    ${_VALGRIND_FLAG}
                    ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR})
endif ()
