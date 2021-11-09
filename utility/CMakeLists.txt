set(UTILITY_SRC
        command_parser.cpp
        command_parser.h
        graph_operations.cpp
        graph_operations.h)

add_library(utility SHARED
        ${UTILITY_SRC})