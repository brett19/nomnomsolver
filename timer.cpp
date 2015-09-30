#define WIN32_LEAN_AND_MEAN
#include <windows.h>

unsigned int hrTime() {
    return GetTickCount();
}

void sleepWait(unsigned int msecs) {
    Sleep(msecs);
}
