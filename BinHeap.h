#ifndef BIN_HEAP_H
#define BIN_HEAP_H

// Binary heap with minimal element in top

template <class ElementValue> class BinHeapMin {
public:
    int capacity;
    int numElems;
    int* indexArray; // Index of heap element in some external array
    int* heapIndex;  // Index of element of external array in heap
    //... double* elements;
    ElementValue* elements; // Array of heap elements

    BinHeapMin(int maxSize);
    void bubbleUp(int i);   // Called after a change (decrement) of element
    void bubbleDown(int i);
    void removeRoot();
    void orderHeap();       // Called for arbitrary contents
    ElementValue root() const { return elements[0]; }
    int rootIndex() const {
        return indexArray[0];
    }

    ~BinHeapMin();

private:
    void swap(int i, int j) {
        //... double t = elements[i];
        ElementValue t = elements[i];
        elements[i] = elements[j];
        elements[j] = t;
        if (indexArray != 0) {
            int ext0 = indexArray[i];
            int ext1 = indexArray[j];
            indexArray[i] = ext1;
            indexArray[j] = ext0;
            if (heapIndex != 0) {
                heapIndex[ext0] = j;
                heapIndex[ext1] = i;
            }
        }
    }
};

#endif
