cmake_minimum_required(VERSION 3.10)
project(openga VERSION 1.0.3 LANGUAGES CXX)

include(CMakePackageConfigHelpers)

find_package(Threads)

# Options
option(Examples "Build the project examples." OFF)

set(OPENGA_TARGET_NAME                ${PROJECT_NAME})
set(OPENGA_CONFIG_INSTALL_DIR         "lib/cmake/${PROJECT_NAME}" CACHE INTERNAL "")
set(OPENGA_INCLUDE_BUILD_DIR          "${PROJECT_SOURCE_DIR}/include/")
set(OPENGA_INCLUDE_INSTALL_DIR        "include")
set(OPENGA_TARGETS_EXPORT_NAME        "${PROJECT_NAME}Targets")
set(OPENGA_CMAKE_CONFIG_TEMPLATE      "cmake/config.cmake.in")
set(OPENGA_CMAKE_CONFIG_DIR           "${CMAKE_CURRENT_BINARY_DIR}")
set(OPENGA_CMAKE_VERSION_CONFIG_FILE  "${OPENGA_CMAKE_CONFIG_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(OPENGA_CMAKE_PROJECT_CONFIG_FILE  "${OPENGA_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Config.cmake")
set(OPENGA_CMAKE_PROJECT_TARGETS_FILE "${OPENGA_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Targets.cmake")

add_library(openga INTERFACE)
target_include_directories(openga INTERFACE 
	$<BUILD_INTERFACE:${NLOHMANN_JSON_INCLUDE_BUILD_DIR}>
	$<INSTALL_INTERFACE:include>
)
target_link_libraries(openga INTERFACE Threads::Threads)

write_basic_package_version_file(
    ${OPENGA_CMAKE_VERSION_CONFIG_FILE} COMPATIBILITY SameMajorVersion
)

configure_file(
    ${OPENGA_CMAKE_CONFIG_TEMPLATE}
    ${OPENGA_CMAKE_PROJECT_CONFIG_FILE}
    @ONLY
)

install(
	DIRECTORY ${OPENGA_INCLUDE_BUILD_DIR}
	DESTINATION ${OPENGA_INCLUDE_INSTALL_DIR}
)

install(
	FILES ${OPENGA_CMAKE_PROJECT_CONFIG_FILE} ${OPENGA_CMAKE_VERSION_CONFIG_FILE}
	DESTINATION ${OPENGA_CONFIG_INSTALL_DIR}
)

export(
	TARGETS ${OPENGA_TARGET_NAME}
	NAMESPACE ${PROJECT_NAME}::
	FILE ${OPENGA_CMAKE_PROJECT_TARGETS_FILE}
)

install(
	TARGETS ${OPENGA_TARGET_NAME}
	EXPORT ${OPENGA_TARGETS_EXPORT_NAME}
	INCLUDES DESTINATION ${OPENGA_INCLUDE_INSTALL_DIR}
)

install(
	EXPORT ${OPENGA_TARGETS_EXPORT_NAME}
	NAMESPACE ${PROJECT_NAME}::
	DESTINATION ${OPENGA_CONFIG_INSTALL_DIR}
)

# build configuration for all the examples
if(Examples)
	set(EXAMPLES_DIR "examples")

	add_executable(so-1 "${EXAMPLES_DIR}/so-1/example_so1.cpp")
	target_include_directories(so-1 PRIVATE ${OPENGA_INCLUDE_BUILD_DIR})
	target_link_libraries(so-1 PRIVATE openga)
	
	add_executable(mo-1 "${EXAMPLES_DIR}/mo-1/example_mo1.cpp")
	target_include_directories(mo-1 PRIVATE ${OPENGA_INCLUDE_BUILD_DIR})
	target_link_libraries(mo-1 PRIVATE openga)
	
	add_executable(mo-dtlz2 "${EXAMPLES_DIR}/mo-dtlz2/mo-dtlz2.cpp")
	target_include_directories(mo-dtlz2 PRIVATE ${OPENGA_INCLUDE_BUILD_DIR})
	target_link_libraries(mo-dtlz2 PRIVATE openga)
	
	add_executable(so-rastrigin "${EXAMPLES_DIR}/so-rastrigin/so-rastrigin.cpp")
	target_include_directories(so-rastrigin PRIVATE ${OPENGA_INCLUDE_BUILD_DIR})
	target_link_libraries(so-rastrigin PRIVATE openga)
endif()
