#pragma once
/*
Clasic haar 2D wavelet

Type    is all signed type;
SizeWH  is number pow 2 and >= 2: 2 4 8 16....

Data, DataTemp, DataWavelet is Type array [SizeWH * SizeWH]

DataWavelet = [[A][B][C][D]], [[B B][C C][D D]], [[B B B B][C C C C][D D D D]], ...

Example using
float DataIn[512 * 512];
float DataOut[512 * 512];
float DataTemp[512 * 512];
float DataInInv[512 * 512];

// image to wavelet data
Wavelet2DHaar<float>(512, DataIn, DataTemp, DataOut);

// wavelet data to image
Wavelet2DHaarInverse<float>(512, DataOut, DataTemp, DataInInv);

// wavelet data to debug image
Wavelet2DHaarDebug<float>(512, DataOut, DataTemp, 127);

*/

inline bool isPow2(unsigned int Val) {
    return ((Val & (~(Val - 1))) == Val);
}

template <typename Type>
inline void Wavelet2DHaar(unsigned int SizeWH, Type* Data, Type* DataTemp, Type* DataWavelet) {
    assert(isPow2(SizeWH));
    assert(SizeWH > 1);
    assert(Data);
    assert(DataTemp);
    assert(DataWavelet);

    Type* DataIn = Data;
    Type* DataOut = DataTemp;

    unsigned int TempOffset = 0;
    unsigned int SizeWH_ = SizeWH;
    unsigned int SizeWH_2 = SizeWH / 2;

    while (SizeWH_ >= 2) {
        unsigned int SizeBlock = SizeWH_2 * SizeWH_2;

        if (SizeWH_ == 2) DataOut = DataWavelet;

        for (unsigned int y = 0, InIDX = 0, OutIDX = 0; y < SizeWH_2; y++) {
            for (unsigned int x = 0; x < SizeWH_2; x++) {
                unsigned int InIDX_2 = InIDX;
                Type A = DataIn[InIDX_2];

                InIDX_2 += 1;
                Type B = DataIn[InIDX_2];

                InIDX_2 += SizeWH_;
                Type D = DataIn[InIDX_2];

                InIDX_2 -= 1;
                Type C = DataIn[InIDX_2];

                Type LL = (A + B + C + D) / (Type)4;
                Type HL = (A - B + C - D) / (Type)4;
                Type LH = (A + B - C - D) / (Type)4;
                Type HH = (A - B - C + D) / (Type)4;

                unsigned int OutIDX_2 = OutIDX;
                DataOut[OutIDX_2] = LL;

                OutIDX_2 += SizeBlock;
                DataWavelet[OutIDX_2] = HL;

                OutIDX_2 += SizeBlock;
                DataWavelet[OutIDX_2] = LH;

                OutIDX_2 += SizeBlock;
                DataWavelet[OutIDX_2] = HH;

                OutIDX += 1;

                InIDX += 2;
            }

            InIDX += SizeWH_;
        }

        DataIn = DataOut;

        TempOffset += SizeBlock;
        DataOut = &DataTemp[TempOffset];

        SizeWH_ = SizeWH_2;
        SizeWH_2 /= 2;
    }
}

template <typename Type>
inline void Wavelet2DHaarInverse(unsigned int SizeWH, Type* DataWavelet, Type* DataTemp, Type* Data) {
    assert(isPow2(SizeWH));
    assert(SizeWH > 1);
    assert(DataWavelet);
    assert(DataTemp);
    assert(Data);

    Type* DataIn = DataWavelet;
    Type* DataOut = DataTemp;

    unsigned int SizeWH_2 = 1;
    unsigned int SizeWH_ = 2;

    while (SizeWH_ <= SizeWH) {
        unsigned int SizeBlock = SizeWH_2 * SizeWH_2;

        for (unsigned int y = 0, InIDX = 0, OutIDX = 0; y < SizeWH_2; y++) {
            for (unsigned int x = 0; x < SizeWH_2; x++) {
                unsigned int OutIDX_2 = OutIDX;
                Type LL = DataIn[OutIDX_2];

                OutIDX_2 += SizeBlock;
                Type HL = DataWavelet[OutIDX_2];

                OutIDX_2 += SizeBlock;
                Type LH = DataWavelet[OutIDX_2];

                OutIDX_2 += SizeBlock;
                Type HH = DataWavelet[OutIDX_2];

                Type A = (LL + HL + LH + HH);
                Type B = (LL - HL + LH - HH);
                Type C = (LL + HL - LH - HH);
                Type D = (LL - HL - LH + HH);

                unsigned int InIDX_2 = InIDX;
                DataOut[InIDX_2] = A;

                InIDX_2 += 1;
                DataOut[InIDX_2] = B;

                InIDX_2 += SizeWH_;
                DataOut[InIDX_2] = D;

                InIDX_2 -= 1;
                DataOut[InIDX_2] = C;

                OutIDX += 1;

                InIDX += 2;
            }

            InIDX += SizeWH_;
        }

        DataIn = DataOut;

        if (DataOut == DataTemp)    DataOut = Data;
        else                        DataOut = DataTemp;

        SizeWH_2 = SizeWH_;
        SizeWH_ *= 2;
    }

    if (Data != DataIn) memcpy(Data, DataIn, sizeof(Type) * SizeWH * SizeWH);
}

// slow function
template <typename Type>
inline void Wavelet2DHaarDebug(unsigned int SizeWH, Type* DataWavelet, Type* Data, Type AddZero = 0) {
    assert(isPow2(SizeWH));
    assert(SizeWH > 1);
    assert(DataWavelet);
    assert(Data);

    unsigned int InIDX = 0;
    Data[0] = DataWavelet[InIDX++];

    for (unsigned int CurSizeWH = 1; CurSizeWH < SizeWH; CurSizeWH *= 2) {
        for (unsigned int y = 0; y < CurSizeWH; y++) {
            for (unsigned int x = 0; x < CurSizeWH; x++) {
                unsigned int OutIDX = SizeWH * (0 + y) + (CurSizeWH + x);
                Data[OutIDX] = DataWavelet[InIDX++] + AddZero;
            }
        }

        for (unsigned int y = 0; y < CurSizeWH; y++) {
            for (unsigned int x = 0; x < CurSizeWH; x++) {
                unsigned int OutIDX = SizeWH * (CurSizeWH + y) + (0 + x);
                Data[OutIDX] = DataWavelet[InIDX++] + AddZero;
            }
        }

        for (unsigned int y = 0; y < CurSizeWH; y++) {
            for (unsigned int x = 0; x < CurSizeWH; x++) {
                unsigned int OutIDX = SizeWH * (CurSizeWH + y) + (CurSizeWH + x);
                Data[OutIDX] = DataWavelet[InIDX++] + AddZero;
            }
        }
    }
}
