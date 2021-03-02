#include <cassert>
#include "BinHeap.h"

template <class ElementValue>
BinHeapMin<ElementValue>::BinHeapMin(int maxSize):
    capacity(maxSize),
    numElems(0),
    indexArray(0),
    elements(new ElementValue[maxSize])
{}

template <class ElementValue>
void BinHeapMin<ElementValue>::bubbleUp(int i) {
    assert(0 <= i && i < numElems);

    if (i <= 0)
        return;
    int parent = (i - 1) / 2;
    while (i > 0 && elements[i] < elements[parent]) {
        swap(i, parent);
        i = parent;
        parent = (i - 1) / 2;
    }
}

template <class ElementValue>
void BinHeapMin<ElementValue>::bubbleDown(int i) {
    int s0 = 2*i + 1;
    int s1 = s0 + 1;
    int s = s0;
    if (s1 < numElems && elements[s1] < elements[s])
        s = s1;
    while (s < numElems && elements[s] < elements[i]) {
        swap(i, s);
        i = s;
        s0 = 2*i + 1;
        s1 = s0 + 1;
        s = s0;
        if (s1 < numElems && elements[s1] < elements[s])
            s = s1;
    }
}

template <class ElementValue>
void BinHeapMin<ElementValue>::removeRoot() {
    assert(numElems > 0);

    elements[0] = elements[numElems - 1];
    if (indexArray != 0) {
        int ext = indexArray[numElems - 1];
        indexArray[0] = ext;
        if (heapIndex != 0)
            heapIndex[ext] = 0;
    }
    --numElems;
    bubbleDown(0);
}

template <class ElementValue>
void BinHeapMin<ElementValue>::orderHeap() { // Called for arbitrary contents
    int k = numElems / 2;
    while (k >= 0) {
        bubbleDown(k);
        --k;
    }
}

template <class ElementValue>
BinHeapMin<ElementValue>::~BinHeapMin() {
    delete[] elements;
}

// Provide explicit instantiations for
// necessary BinHeapMin classes

//... template class BinHeapMin<int>;
template class BinHeapMin<char>;
template class BinHeapMin<double>;
