# ESP-IDF component: ARBITER orientation filter.
# Place this directory under <project>/components/arbiter (or add it to
# EXTRA_COMPONENT_DIRS) and add "arbiter" to the REQUIRES of your main
# component.
idf_component_register(
    SRCS "arbiter.c"
    INCLUDE_DIRS "include"
)
