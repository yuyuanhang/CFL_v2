set(MATCHING_SRC
        matching_command.cpp
        matching_command.h
        filter_vertices.cpp
        filter_vertices.h
        generate_filtering_plan.cpp
        generate_filtering_plan.h
        build_table.cpp
        build_table.h
        generate_query_plan.cpp
        generate_query_plan.h
        evaluate_query.cpp
        evaluate_query.h)

add_library(matching SHARED
        ${MATCHING_SRC})