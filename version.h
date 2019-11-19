#ifndef VERSION_H
#define VERSION_H
namespace AutoVersion{
    //Date Version Types
    static const char DATE[] = "16";
    static const char MONTH[] = "Nov";
    static const char YEAR[] = "2019";
    static const char UBUNTU_VERSION_STYLE[] = "0.6.0";
    //Software Status
    static const char STATUS[] = "Alpha";
    static const char STATUS_SHORT[] = "a";
    //Standard Version Type
    static const long MAJOR = 7;
    static const long MINOR = 1;
    static const long BUILD = 0;
    static const long REVISION = 1;
    //Miscellaneous Version Types
    static const long BUILDS_COUNT = 2;
    #define RC_FILEVERSION 1,0,0,0
    #define RC_FILEVERSION_STRING "1, 0, 0, 0\0"
    static const char FULLVERSION_STRING[] = "0.7.1.1";

    //These values are to keep track of your versioning state, don't modify them.
    static const long BUILD_HISTORY = 2;

}
#endif //VERSION_H

