# SPDX-License-Identifier: Unlicense

target_sources( ${PROJECT_NAME}
    PRIVATE
        RegisterExtension.cpp
        WiggleKitServer.h
        WiggleKitServer.cpp
)

target_include_directories( ${PROJECT_NAME}
    PRIVATE
        "src"
)
