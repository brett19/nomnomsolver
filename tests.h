#ifdef TESTNUMBER

/*
GRID startH = {
    {GV(7,1), GV(8,1), GV(9,1), GV(7,1), GV(8,1), GV(9,1)},
    {GV(7,1), GV(1,1), GV(2,1), GV(3,1), GV(0,0), GV(0,0)},
    {GV(7,1), GV(4,1), GV(5,1), GV(6,1), GV(0,0), GV(0,0)},
    {GV(7,1), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)}
};
GRID startT = {
    {GV(0,0), GV(2,1), GV(3,1), GV(0,0), GV(4,1), GV(5,1)},
    {GV(6,1), GV(7,1), GV(8,1), GV(9,1), GV(10,1), GV(11,1)},
    {GV(12,1), GV(0,0), GV(2,1), GV(3,1), GV(0,0), GV(0,0)},
    {GV(0,0), GV(6,1), GV(7,1), GV(8,1), GV(9,1), GV(10,1)}
};
*/

#if TESTNUMBER == 1 // Simple - Top Covered, Simple Rearrange [18]
GRID start = {
    {GV(7,1), GV(7,1), GV(7,1), GV(7,1), GV(7,1), GV(7,1)},
    {GV(0,0), GV(1,1), GV(2,1), GV(3,1), GV(0,0), GV(0,0)},
    {GV(0,0), GV(4,1), GV(5,1), GV(6,1), GV(0,0), GV(0,0)},
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)}
};
char i[][2][2] = {
    {{6,1}, {3,1}},
    {{4,1}, {2,1}}
};
#elif TESTNUMBER == 2 // Medium - Surrounded [25]
GRID start = {
    {GV(1,1), GV(1,1), GV(1,1), GV(1,1), GV(1,1), GV(1,1)},
    {GV(1,1), GV(2,1), GV(3,1), GV(4,1), GV(5,1), GV(1,1)},
    {GV(1,1), GV(5,1), GV(4,1), GV(3,1), GV(2,1), GV(1,1)},
    {GV(1,1), GV(1,1), GV(1,1), GV(1,1), GV(1,1), GV(1,1)},
};
char i[2][2][2] = {
    {{4,2}, {5,2}},
    {{2,2}, {3,2}}
};
#elif TESTNUMBER == 3 // Stupid Hard - One Space
GRID start = {
    {GV(1,1), GV(2,1), GV(3,1), GV(4,1), GV(5,1), GV(2,1)},
    {GV(0,0), GV(6,1), GV(2,1), GV(7,1), GV(8,1), GV(8,1)},
    {GV(9,1), GV(4,1), GV(10,1), GV(1,1), GV(4,1), GV(6,1)},
    {GV(2,1), GV(7,1), GV(3,1), GV(7,1), GV(2,1), GV(5,1)}
};
char i[][2][2] = {
    {{2,5}, {5,2}},
    {{4,3}, {3,2}}
};
#elif TESTNUMBER == 4
GRID start = {
    {GV(1,1), GV(2,1), GV(3,1), GV(4,1), GV(5,1), GV(6,1)},
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)},
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)},
    {GV(7,1), GV(7,1), GV(8,1), GV(10,1), GV(11,1), GV(12,1)}
};
char i[][2][2] = {
    {{1,1}, {2,1}},
    {{3,1}, {4,1}},
    {{5,1}, {6,1}},
    {{7,2}, {8,1}}
};
#elif TESTNUMBER == 5
GRID start = {
    {GV(1,1), GV(2,1), GV(3,1), GV(4,1), GV(5,1), GV(2,1)},
    {GV(0,0), GV(6,1), GV(2,1), GV(7,1), GV(8,1), GV(8,1)},
    {GV(9,1), GV(4,1), GV(10,1), GV(1,1), GV(4,1), GV(6,1)},
    {GV(2,1), GV(7,1), GV(3,1), GV(7,1), GV(2,1), GV(5,1)}
};
char i[][2][2] = {
    {{2,5}, {5,2}},
    {{4,3}, {3,2}},
    {{7,3}, {8,2}},
    {{10,1}, {1,2}}
};
#elif TESTNUMBER == 6
GRID start = {
    {GV(0,0), GV(1,1), GV(2,1), GV(5,1), GV(1,1), GV(0,0)},
    {GV(1,1), GV(1,1), GV(3,1), GV(4,1), GV(1,1), GV(1,1)},
    {GV(1,1), GV(1,1), GV(4,1), GV(3,1), GV(1,1), GV(1,1)},
    {GV(0,0), GV(1,1), GV(2,1), GV(5,1), GV(1,1), GV(0,0)},
};
char i[][2][2] = {
    {{3,2}, {5,2}},
    {{4,2}, {2,2}}
};
#elif TESTNUMBER == 7 // ? - Tom Random [30]
GRID start = {
    {GV(11,1), GV(7,1), GV(7,1), GV(0,0), GV(0,0), GV(3,1)},
    {GV(1,1), GV(10,1), GV(9,1), GV(7,1), GV(6,1), GV(9,1)},
    {GV(11,1), GV(1,1), GV(5,1), GV(0,0), GV(0,0), GV(12,1)},
    {GV(4,1), GV(7,1), GV(3,1), GV(12,1), GV(4,1), GV(11,1)}
};

char i[][2][2] = {
    {{3,1}, {9,1}},
    {{9,1}, {11,2}}
};
#elif TESTNUMBER == 8 // Fast - Tom, 1 Player [22]
GRID start = {
    {GV(0,0), GV(1,1), GV(2,1), GV(1,1), GV(3,1), GV(0,0)},
    {GV(0,0), GV(1,1), GV(3,1), GV(1,1), GV(2,1), GV(0,0)},
    {GV(0,0), GV(1,1), GV(2,1), GV(1,1), GV(3,1), GV(0,0)},
    {GV(0,0), GV(1,1), GV(3,1), GV(1,1), GV(2,1), GV(0,0)}
};

char i[][2][2] = {
    {{3,4}, {2,4}}
};
#elif TESTNUMBER == 9 // ? - Tom Random 
GRID start = {
    {GV(0,0), GV(0,0), GV(11,1), GV(9,1), GV(9,1), GV(10,1)},
    {GV(10,1), GV(11,1), GV(12,1), GV(4,1), GV(10,1), GV(2,1)},
    {GV(12,1), GV(10,1), GV(0,0), GV(2,1), GV(12,1), GV(0,0)},
    {GV(4,1), GV(11,1), GV(2,1), GV(12,1), GV(3,1), GV(3,1)}
};

char i[][2][2] = {
    {{3,1}, {9,1}},
    {{2,2}, {11,2}},
    {{12,2}, {11,1}},
    {{10,3}, {12,2}}
};
#elif TESTNUMBER == 10 // ? - Real Level 5
GRID start = {
    {GV(5,1), GV(5,1), GV(0,0), GV(0,0), GV(5,1), GV(5,1)},
    {GV(6,1), GV(5,1), GV(1,1), GV(4,1), GV(5,1), GV(6,1)},
    {GV(6,1), GV(5,1), GV(3,1), GV(2,1), GV(5,1), GV(6,1)},
    {GV(5,1), GV(5,1), GV(0,0), GV(0,0), GV(5,1), GV(5,1)}
};

char i[][2][2] = {
    {{1,1}, {2,1}},
    {{3,1}, {4,1}}
};
#endif

const int ttlPlyrs = sizeof(i) / sizeof(char[2][2]);

GRID winTest = {
    {GV(1,1), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(5,1)},
    {GV(2,1), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(6,1)},
    {GV(3,1), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(7,2)},
    {GV(4,1), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(8,1)}
};

GRID distTest = {
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)},
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)},
    {GV(0,0), GV(2,3), GV(0,0), GV(2,2), GV(0,0), GV(2,2)},
    {GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0), GV(0,0)}
};

#endif
