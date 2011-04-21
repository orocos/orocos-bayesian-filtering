ADD_SUBDIRECTORY( matrix )
ADD_SUBDIRECTORY( rng )

### Generate config.h : HAS NO ACCESS TO SUBDIRECTORY VARIABLES !!!
CONFIGURE_FILE( ${PROJ_SOURCE_DIR}/src/wrappers/config.h.in
    ${PROJ_BINARY_DIR}/src/wrappers/config.h)

