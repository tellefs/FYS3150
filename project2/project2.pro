TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


SOURCES += \
    main.cpp

HEADERS +=

LIBS += -llapack -lblas -larmadillo
