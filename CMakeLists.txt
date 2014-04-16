cmake_minimum_required(VERSION 2.8.9)

project(QuantitativeIndicesExt)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/4.3/Extensions/QuantitativeIndicesExt")
set(EXTENSION_CATEGORY "Quantification")
set(EXTENSION_CONTRIBUTORS "Ethan Ulrich, Markus Van Tol, Reinhard Beichel (Univ. of Iowa)")
set(EXTENSION_DESCRIPTION "Slicer extension used to calulate quantitative features from an input volume and label map.")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/QuantitativeIndicesExt/QuantitativeIndicesExt.png")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-QuantitativeIndicesExt-screenshot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
SlicerMacroExtractRepositoryInfo(VAR_PREFIX QuantitativeIndicesExt)
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/CMake/vtkQuantitativeIndicesExtVersionConfigure.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/vtkQuantitativeIndicesExtVersionConfigure.h
    )

#-----------------------------------------------------------------------------
add_subdirectory(QuantitativeIndicesCLI)
add_subdirectory(QuantitativeIndicesTool)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
