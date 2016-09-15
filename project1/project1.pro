TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    special.cpp \
    error.cpp \
    LU_decomp.cpp

HEADERS += \
    LU_decomp.h

LIBS += -larmadillo -llapack -lblas
