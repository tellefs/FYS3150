TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt


SOURCES += \
    main.cpp

HEADERS +=

LIBS += -llapack -lblas -larmadillo

macx: LIBS += -L$$PWD/../../../../../../../usr/local/Cellar/armadillo/7.400.2/lib/ -larmadillo.7.40.2

INCLUDEPATH += $$PWD/../../../../../../../usr/local/Cellar/armadillo/7.400.2/include
DEPENDPATH += $$PWD/../../../../../../../usr/local/Cellar/armadillo/7.400.2/include
