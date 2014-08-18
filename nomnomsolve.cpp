// nomnomsolve.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <unordered_map>
#include <list>
#include <map>
#include <deque>
#include <queue>
#include <array>
#include <assert.h>
#include <boost/pool/object_pool.hpp>

typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;

#define TRACKPARENTS 1 // Whether to keep track of parents and print full game solutions (costs memory)
#define ESTIMBEST 400 // The estimated best move count
#define MAXCOST 256 // The highest possible cost returned from the heuristic
#define TESTNUMBER 9 // The test to run
#define MAXPLRS 4 // The most players possible in a single game
#define MAXSTACK 6 // Used for optimizing around maximum stacked items
#define CULLBMAP 1 // Use Seen Map
#define THREADBATCH 10000

typedef uint8_t GRID[4][6];
#define GV(itm, num) (itm|(num<<5))
#define GI(id,y,x) (id[y][x]&0x1F)
#define GN(id,y,x) (id[y][x]>>5)
#define GE(id,y,x) (id[y][x]==0)
#define GS(id,y,x,itm,num) id[y][x]=GV(itm,num)
#define GR(id,y,x) id[y][x]=0

template<typename T>
class cppalloc {
public:
    T * malloc() { return (T*)new T; }
    void free(T * p) { delete p; }
};

template<typename T>
class boostalloc {
    boost::object_pool<T> ba;
public:
    T * malloc() { return ba.malloc(); }
    void free(T * p) { ba.free(p); }
};


cppalloc<class AStarOpen> nodePool;
cppalloc<GRID> gridPool;

#include "tests.h"

uint32_t lastStatus = 0;
volatile unsigned long powCount = 0;
volatile unsigned long poofCount = 0;
volatile unsigned long earlyCount = 0;
volatile unsigned long testCount = 0;
volatile unsigned long solutionCount = 0;
volatile unsigned long pendCount = 0;

class GridCmp {
public:
    bool operator()(const GRID* a, const GRID* b) {
        return memcmp(a, b, sizeof(GRID)) < 0;
    }
};

class AStarOpen {
public:
    GRID grd;

    volatile unsigned long refcount;
    unsigned char mindist;
    unsigned char moves;

#if TRACKPARENTS
    AStarOpen *parent;
#endif // TRACKPARENTS

    void addref() {
#pragma omp atomic
        refcount++;
    }
    void decref() {
#pragma omp atomic
        refcount--;
        if (refcount == 0) {
#if TRACKPARENTS
            if (parent) {
                parent->decref();
                parent = nullptr;
            }
#endif
            nodePool.free(this);
#pragma omp atomic
            poofCount++;
        }
    }

};

typedef std::map<GRID*, uint8_t, GridCmp> GRIDBMAP;
typedef std::deque<AStarOpen*> OPENLISTGRP;
typedef std::vector<AStarOpen*> OPENLIST;
GRIDBMAP seenBMap;
OPENLISTGRP openListX[MAXCOST];
volatile unsigned long lowOpenCost = MAXCOST;
volatile unsigned long fastestFound = ESTIMBEST;

struct ThreadData {
    OPENLIST curOpenList;
    OPENLIST newOpenList;
    GRIDBMAP seenList;
} *tdata;
#pragma omp threadprivate(tdata)

const int pxs[5][4][2][2] = {
    {
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,1}, {0,2}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,1}, {0,2}},
        {{5,1}, {5,2}},
        {{0,0}, {0,0}},
        {{0,0}, {0,0}}
    }, {
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}},
        {{0,0},{0,0}}
    }, {
        {{0,0}, {0,1}},
        {{0,2}, {0,3}},
        {{5,0}, {5,1}},
        {{5,2}, {5,3}}
    }
};


uint32_t startTime;
int xtime() {
    return GetTickCount() - startTime;
}

void copyGrid(GRID& dest, const GRID& src) {
    memcpy(&dest, &src, sizeof(GRID));
}

bool checkWin(GRID& in) {
    for (int iq = 0; iq < TTLPLAYRS; ++iq) {
        const int px1 = pxs[TTLPLAYRS][iq][0][0];
        const int py1 = pxs[TTLPLAYRS][iq][0][1];
        const int px2 = pxs[TTLPLAYRS][iq][1][0];
        const int py2 = pxs[TTLPLAYRS][iq][1][1];
        const uint8_t i1i = i[iq][0][0];
        const uint8_t i1n = i[iq][0][1];
        const uint8_t i2i = i[iq][1][0];
        const uint8_t i2n = i[iq][1][1];

        if (GI(in,py1,px1) != i1i && GI(in,py1,px1) != i2i) return false;
        if (GI(in,py2,px2) != i2i && GI(in,py2,px2) != i1i) return false;

        if (GN(in,py1,px1) != i1n && GN(in,py1,px1) != i2n) return false;
        if (GN(in,py2,px2) != i2n && GN(in,py2,px2) != i1n) return false;
    }
    return true;
}

void moveItem(GRID& in, int fx, int fy, int tx, int ty) {
    GS(in,ty,tx, GI(in,fy,fx), GN(in,fy,fx)+GN(in,ty,tx));
    GR(in,fy,fx);
}

void printBoard(GRID& in, AStarOpen *parent);

struct FoundItem {
    unsigned char x;
    unsigned char y;
    char num;
};

int findClosestItem(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num);

struct CacheTables {
    unsigned char t2[MAXSTACK] [2][4] [6][4][MAXSTACK] [6][4][MAXSTACK];
    unsigned char t3[MAXSTACK] [2][4] [6][4][MAXSTACK] [6][4][MAXSTACK] [6][4][MAXSTACK];
} caches;

void buildTable() {
    FILE *rfh = fopen("table.cache", "rb");
    if (rfh != NULL) {
        fread(&caches, sizeof(CacheTables), 1, rfh);
        fclose(rfh);
        return;
    }

    int xlkp[] = { 0, 5 };
    FoundItem itms[3];

    for (int tn = 0; tn < MAXSTACK; ++tn) {
        for (int _sx = 0; _sx < 2; ++_sx) {
            int sx = xlkp[_sx];
            for (int sy = 0; sy < 4; ++sy) {

                for (int ix1 = 0; ix1 < 6; ++ix1) {
                    itms[0].x = ix1;
                    for (int iy1 = 0; iy1 < 4; ++iy1) {
                        itms[0].y = iy1;
                        for (int in1 = 0; in1 < MAXSTACK; ++in1) {
                            itms[0].num = in1;
                            for (int ix2 = 0; ix2 < 6; ++ix2) {
                                itms[1].x = ix2;
                                for (int iy2 = 0; iy2 < 4; ++iy2) {
                                    itms[1].y = iy2;
                                    for (int in2 = 0; in2 < MAXSTACK; ++in2) {
                                        itms[1].num = in2;

                                        caches.t2[tn][_sx][sy][ix1][iy1][in1][ix2][iy2][in2] = 
                                            findClosestItem(itms, 2, 0, 0, -1, sx, sy, tn);

                                        for (int ix3 = 0; ix3 < 6; ++ix3) {
                                            itms[2].x = ix2;
                                            for (int iy3 = 0; iy3 < 4; ++iy3) {
                                                itms[2].y = iy2;
                                                for (int in3 = 0; in3 < MAXSTACK; ++in3) {
                                                    itms[2].num = in2;

                                                    caches.t3[tn][_sx][sy][ix1][iy1][in1][ix2][iy2][in2][ix3][iy3][in3] = 
                                                        findClosestItem(itms, 3, 0, 0, -1, sx, sy, tn);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

            }
        }
    }

    FILE *wfh = fopen("table.cache", "wb");
    if (wfh != NULL) {
       fwrite(&caches, sizeof(CacheTables), 1, wfh);
        fclose(wfh);
    }
}

int findClosestItem(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num) {
    int bestDist = -1;
    unsigned int key = 1;
    const FoundItem *mlist = list;
    for (; mlist->num > 0; mlist++, key <<= 1) {
        if (ignore&key) continue;
        const FoundItem& itm = *mlist;

        int nnum = num - itm.num;
        if (nnum < 0) continue;

        int ndist = abs(sx-itm.x) + abs(sy-itm.y);
        if (nnum == 0) return dist + ndist;
        if (bDist > -1 && dist + ndist >= bDist) {
            continue;
        }

        int xdist = findClosestItem(list, numList, ignore|key, dist+ndist, bestDist, itm.x, itm.y, nnum);
        if (xdist >= 0 && (bestDist == -1 || xdist < bestDist)) {
            bestDist = xdist;
        }
    }
    return bestDist;
}

int findClosestItemFast(const FoundItem *list, int numList, unsigned int ignore, int dist, int bDist, int sx, int sy, int num) {
#if 1
    const int xlkp[] = {0, -1, -1, -1, -1, 5};
    if (numList == 1) {
        if (list->num == num) {
            return abs(sx-list->x) + abs(sy-list->y);
        } else {
            return -1;
        }
    } else if (numList == 2) {
        return caches.t2[num]
            [xlkp[sx]][sy]
            [list[0].x][list[0].y][list[0].num]
            [list[1].x][list[1].y][list[1].num];
    } else if (numList == 3) {
        return caches.t3[num]
            [xlkp[sx]][sy]
            [list[0].x][list[0].y][list[0].num]
            [list[1].x][list[1].y][list[1].num]
            [list[2].x][list[2].y][list[2].num];
    } else {
        return findClosestItem(list, numList, ignore, dist, bDist, sx, sy, num);
    }
#else
    return findClosestItem(list, numList, ignore, dist, bDist, sx, sy, num);
#endif
}

int findClosestPlrItem(const GRID& in, int iq) {
    const int sx1 = pxs[TTLPLAYRS][iq][0][0];
    const int sy1 = pxs[TTLPLAYRS][iq][0][1];
    const int itm1 = i[iq][0][0];
    const int num1 = i[iq][0][1];
    const int sx2 = pxs[TTLPLAYRS][iq][1][0];
    const int sy2 = pxs[TTLPLAYRS][iq][1][1];
    const int itm2 = i[iq][1][0];
    const int num2 = i[iq][1][1];

    int itemCnt1 = 0;
    FoundItem items1[MAXSTACK];
    int itemCnt2 = 0;
    FoundItem items2[MAXSTACK];
    
    memset(items1, 0, sizeof(items1));
    memset(items2, 0, sizeof(items2));

    for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 6; ++ix) {
            if (GI(in,iy,ix) == itm1) {
                items1[itemCnt1].x = ix;
                items1[itemCnt1].y = iy;
                items1[itemCnt1].num = GN(in,iy,ix);
                itemCnt1++;
            }
            if (GI(in,iy,ix) == itm2) {
                items2[itemCnt2].x = ix;
                items2[itemCnt2].y = iy;
                items2[itemCnt2].num = GN(in,iy,ix);
                itemCnt2++;
            }
        }
    }

    int cost11 = findClosestItem(items1, itemCnt1, 0, 0, -1, sx1, sy1, num1);
    int cost12 = findClosestItem(items1, itemCnt1, 0, 0, -1, sx2, sy2, num1);
    if (cost11 == -1 && cost12 == -1) return -1;

    int cost21 = findClosestItem(items2, itemCnt2, 0, 0, -1, sx2, sy2, num2);
    int cost22 = findClosestItem(items2, itemCnt2, 0, 0, -1, sx1, sy1, num2);
    if (cost21 == -1 && cost22 == -1) return -1;

    int ttlCost = 0;
    if (cost11 == -1 || cost12 < cost11) {
        ttlCost += cost12;
    } else {
        ttlCost += cost11;
    }
    if (cost21 == -1 || cost22 < cost21) {
        ttlCost += cost22;
    } else {
        ttlCost += cost21;
    }
    return ttlCost;
}

int calcDist(const GRID& in) {
    int ttlCost = 0;
    for (int iq = 0; iq < TTLPLAYRS; ++iq) {
        int cost = findClosestPlrItem(in, iq);
        if (cost == -1) {
            return -1;
        }
        ttlCost += cost;
    }
    assert(ttlCost < MAXCOST);
    return ttlCost;
}

void printBoard(GRID& in, AStarOpen *parent) {
    for (int iy = 0; iy < 4; ++iy) {
        printf("| ");
        for (int ix = 0; ix < 6; ++ix) {
            printf("%02dx%02d ", GN(in,iy,ix), GI(in,iy,ix));
        }

        if (parent) {
            printf ("| | ");
            GRID& pg = parent->grd;
            for (int ipx = 0; ipx < 6; ipx++) {
                if (GI(in,iy,ipx) != GI(pg,iy,ipx) || GN(in,iy,ipx) != GN(pg,iy,ipx)) {
                    if (GE(in,iy,ipx)) {
                        printf("----- ");
                    } else {
                        printf("%02dx%02d ", GN(in,iy,ipx), GI(in,iy,ipx));
                    }
                } else {
                    printf("      ");
                }
            }
        }

        printf("|\n");
    }
}

#if TRACKPARENTS
void recursePrintBoard(AStarOpen *itm) {
    if (itm->parent) {
        recursePrintBoard(itm->parent);
    }

    printf("\\_____________________________________/\n");
    printBoard(itm->grd, itm->parent);
    printf("\\_____________________________________/ \\_____________________________________/\n");
}
#else
void recursePrintBoard(AStarOpen *itm) {
    printf("\\_____________________________________/\n");
    printBoard(itm->grd, nullptr);
    printf("\\_____________________________________/\n");
}
#endif

#if CULLBMAP && _OPENMP
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
    GRIDBMAP::iterator seenItm = seenBMap.find((GRID*)&grd);
    if (seenItm != seenBMap.end()) {
        if (seenItm->second <= moves) {
            return true;
        }
        seenItm->second = moves;
        return false;
    }
    seenItm = tdata->seenList.find((GRID*)&grd);
    if (seenItm != tdata->seenList.end()) {
        if (seenItm->second <= moves) {
            return true;
        }
        seenItm->second = moves;
        return false;
    }

    GRID *ngrid = gridPool.malloc();
    memcpy(ngrid, &grd, sizeof(GRID));
    tdata->seenList.insert(std::pair<GRID*,unsigned int>(ngrid, moves));
    return false;
}
#elif CULLBMAP
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
    GRIDBMAP::iterator seenItm = seenBMap.find((GRID*)&grd);
    if (seenItm != seenBMap.end()) {
        //assert(memcmp(seenItm->first, &itm.grd, sizeof(GRID)) == 0);
        if (seenItm->second <= moves) {
            return true;
        }
        seenItm->second = moves;
    } else {
        GRID *ngrid = gridPool.malloc();
        memcpy(ngrid, &grd, sizeof(GRID));
        seenBMap.insert(std::pair<GRID*,unsigned int>(ngrid, moves));
    }
    return false;
}
#else
__forceinline bool wasGridSeen(const GRID& grd, int moves) {
    return false;
}
#endif

void addToOpen(const GRID& grd, unsigned int moves, AStarOpen *parent) {
    int minCost = calcDist(grd);
    if (minCost == -1) {
        earlyCount++;
        return;
    }

    assert(minCost >= 0 && minCost < MAXCOST);
    assert(moves + minCost < MAXCOST);

    if (moves + minCost >= fastestFound) {
        earlyCount++;
        return;
    }

    if (wasGridSeen(grd, moves)) {
        earlyCount++;
        return;
    }

    AStarOpen *itm = nodePool.malloc();

#pragma omp atomic
    powCount++;

    copyGrid(itm->grd, grd);
    itm->mindist = minCost;
    itm->moves = moves;
    itm->refcount = 1;

#if TRACKPARENTS
    if (parent) {
        parent->addref();
        itm->parent = parent;
    } else {
        itm->parent = nullptr;
    }
#endif

#if _OPENMP
    tdata->newOpenList.push_back(itm);
#else
    openListX[minCost].push_back(itm);
    if (minCost < lowOpenCost) {
        lowOpenCost = minCost;
    }
#endif
}

int addMoves(AStarOpen *owner) {
    GRID & in = owner->grd;
    int moves = owner->moves;

    testCount++;
#if !_OPENMP

    static unsigned int lasttick = 0;
    if (testCount % 10000 == 0) {
        unsigned int thistick = GetTickCount();

        int kmps = -1;
        if (thistick - lasttick > 0) {
            kmps = 10000 / (thistick - lasttick);
            lasttick = thistick;
        }
        printf("[% 7d] tested %dk (%d best) [% 5dkmps %dk,%dk,%dk,%dk]\n", xtime(), testCount/1000, solutionCount?fastestFound:-1, kmps, powCount/1000, poofCount/1000, (powCount-poofCount)/1000, earlyCount/1000);
    }
#endif

    if (owner->moves + owner->mindist >= fastestFound) {
        return 3;
    }

    if (checkWin(in)) {
        recursePrintBoard(owner);
        printf("[% 7d] Found winner using %d moves!\n", xtime(), moves);

        fastestFound = moves;
        solutionCount++;

#if 0
        for (int z = 255; z >= fastestFound; --z) {
            for(OPENLISTGRP::iterator j = openListX[z].begin(); j != openListX[z].end(); ++j) {
                (*j)->decref();
            }
            openListX[z].clear();
        }
#endif

        return 1000;
    }

    GRID tgrd;

    for (int iy = 0; iy < 4; ++iy) {
        for (int ix = 0; ix < 6; ++ix) {
            if (GE(in,iy,ix)) {
                // Cell is empty, can't do anything
                continue;
            }

            // Up
            if (iy >= 1) {
                if (GE(in,iy-1,ix) || GI(in,iy-1,ix) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(tgrd, in);
                    moveItem(tgrd, ix, iy, ix, iy-1);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Down
            if (iy <= 2) {
                if (GE(in,iy+1,ix) || GI(in,iy+1,ix) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(tgrd, in);
                    moveItem(tgrd, ix, iy, ix, iy+1);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Left
            if (ix >= 1) {
                if (GE(in,iy,ix-1) || GI(in,iy,ix-1) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(tgrd, in);
                    moveItem(tgrd, ix, iy, ix-1, iy);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }

            // Right
            if (ix <= 4) {
                if (GE(in,iy,ix+1) || GI(in,iy,ix+1) == GI(in,iy,ix)) {
                    // Move
                    copyGrid(tgrd, in);
                    moveItem(tgrd, ix, iy, ix+1, iy);
                    addToOpen(tgrd, owner->moves+1, owner);
                }
            }
        }
    }

    return 0;
}

void doSearch() {
#if _OPENMP
#pragma omp parallel
    while(true) {
        uint32_t batchTotal = 0;

#pragma omp critical
        {
            // Merge everything back to global state, we do this to pull
            //   anything added by the master thread in the beginning into
            //   the global state
            for (OPENLIST::iterator i = tdata->newOpenList.begin(); i != tdata->newOpenList.end(); ++i) {
                openListX[(*i)->mindist].push_back(*i);
                if ((*i)->mindist < lowOpenCost) {
                    lowOpenCost = (*i)->mindist;
                }
            }
            tdata->newOpenList.clear();

            seenBMap.insert(tdata->seenList.begin(), tdata->seenList.end());
            tdata->seenList.clear();

            // Create a list of things for this thread
            batchTotal = 0;
            while (lowOpenCost < 256 && batchTotal < THREADBATCH) {
                if (openListX[lowOpenCost].empty()) {
                    lowOpenCost++;
                    continue;
                }

                AStarOpen* itm = openListX[lowOpenCost].front();
                openListX[lowOpenCost].pop_front();

                tdata->curOpenList.push_back(itm);
                batchTotal++;
            }
        }

#pragma omp atomic
        pendCount += batchTotal;

        uint32_t curTime = GetTickCount();
        if (curTime - lastStatus > 1000) {
            lastStatus = curTime;
            printf("[% 7d] tested %dk (%d best) [%dk %dk %dk %dk %dk]\n", 
                xtime(),
                testCount/1000,
                solutionCount?fastestFound:-1,
                powCount/1000,
                poofCount/1000,
                (powCount-poofCount)/1000,
                earlyCount/1000,
                pendCount/1000);
        }

        if (batchTotal == 0) {
#if 0
            printf("T%d ran out of items\n", omp_get_thread_num());
            break;
#else
            if (pendCount == 0) {
                break;
            }
            Sleep(10);
            continue;
#endif
        }

        //printf("[% 7d] T%d executing %d items\n", xtime(), omp_get_thread_num(), tdata->curOpenList.size());

        while (!tdata->curOpenList.empty()) {
            AStarOpen *itm = tdata->curOpenList.back();
            tdata->curOpenList.pop_back();

            addMoves(itm);
            itm->decref();
        }

        //printf("[% 7d] T%d finished executing\n", xtime(), omp_get_thread_num());

#pragma omp atomic
        pendCount -= batchTotal;
    }

#else
    while (lowOpenCost < 256) {
        if (openListX[lowOpenCost].empty()) {
            lowOpenCost++;
            continue;
        }

        AStarOpen* itm = openListX[lowOpenCost].front();
        openListX[lowOpenCost].pop_front();

        addMoves(itm);
        itm->decref();
    }
#endif
}

int _tmain(int argc, _TCHAR* argv[])
{
    omp_set_num_threads(omp_get_num_procs()-1);

    startTime = GetTickCount();

    printf("ESTIMBEST: %d\n", ESTIMBEST);
    printf("MAXCOST: %d\n", MAXCOST);
    printf("TESTNUMBER: %d\n", TESTNUMBER);
    printf("MAXPLRS: %d\n", MAXPLRS);
    printf("MAXSTACK: %d\n", MAXSTACK);
    printf("CULLBMAP: %d\n", CULLBMAP);
    printf("_OPENMP: %d\n", _OPENMP);
    printf("TTLPLAYRS: %d\n", TTLPLAYRS);

    printf("sizeof(GRID): %d\n", sizeof(GRID));
    printf("sizeof(AStarOpen): %d\n", sizeof(AStarOpen));
    printf("sizeof(openListX): %d\n", sizeof(openListX));
    printf("\n");

    printf("checkWin(winTest): %d\n", checkWin(winTest));
    printf("checkWin(start): %d\n", checkWin(start));
    printf("\n");

    buildTable();
    printf("[% 7d] Cache Map Ready\n", xtime());

#if _OPENMP
#pragma omp parallel
    tdata = new ThreadData();
#endif

    addToOpen(start, 0, nullptr);
    doSearch();
    printf("[% 7d] Completed! [%d %d %d %d]\n", xtime(), powCount, poofCount, powCount-poofCount, earlyCount);

#if _OPENMP
#pragma omp parallel
    delete tdata;
#endif

    system("PAUSE");
    return 0;
}

