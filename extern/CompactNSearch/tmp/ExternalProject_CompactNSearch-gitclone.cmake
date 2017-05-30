if("1.0.0" STREQUAL "")
  message(FATAL_ERROR "Tag for git checkout should not be empty.")
endif()

set(run 0)

if("D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitinfo.txt" IS_NEWER_THAN "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitclone-lastrun.txt")
  set(run 1)
endif()

if(NOT run)
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch'")
endif()

# try the clone 3 times incase there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "d:/Program Files/Git/cmd/git.exe" clone "https://github.com/InteractiveComputerGraphics/CompactNSearch.git" "ExternalProject_CompactNSearch"
    WORKING_DIRECTORY "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/InteractiveComputerGraphics/CompactNSearch.git'")
endif()

execute_process(
  COMMAND "d:/Program Files/Git/cmd/git.exe" checkout 1.0.0
  WORKING_DIRECTORY "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: '1.0.0'")
endif()

execute_process(
  COMMAND "d:/Program Files/Git/cmd/git.exe" submodule init
  WORKING_DIRECTORY "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to init submodules in: 'D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch'")
endif()

execute_process(
  COMMAND "d:/Program Files/Git/cmd/git.exe" submodule update --recursive 
  WORKING_DIRECTORY "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitinfo.txt"
    "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitclone-lastrun.txt"
  WORKING_DIRECTORY "D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'D:/Projects/MySPlishSPlasH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-stamp/ExternalProject_CompactNSearch-gitclone-lastrun.txt'")
endif()

