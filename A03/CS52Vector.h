// Michael Helfand
// 1653490
// Assignment #: 3
// Submission Date: 11/10/2024
//
// Description: A program that uses the struct vector and helper functions to act as the C++ library std::vector
#include <cstdlib>  //C++ version of stdlib.h
#include <cstring>  //C++ version of string.h
#include <string>   //C++ string library
#include <iostream>
#include <stdexcept>
#include <algorithm>

struct Vector {
    int size = 0;
    int capacity = 0;
    int* data = nullptr;
};

// Construct a Vector with a specified size and initialize with `initVal`.
void construct_vector(Vector& v, int size = 0, int initVal = 0) {
    v.size = size;
    v.capacity = size > 0 ? size : 1; // Ensure at least capacity of 1
    v.data = new int[v.capacity];
    for (int i = 0; i < v.size; ++i) {
        v.data[i] = initVal;
    }
}

// Destroy the vector and free allocated memory.
void destroy_vector(Vector& v) {
    delete[] v.data;
    v.data = nullptr;
    v.size = 0;
    v.capacity = 0;
}

// Returns a pointer to the first element in the vector.
int* data(const Vector& v) {
    return v.data;
}

// Returns the number of elements in the vector.
int size(const Vector& v) {
    return v.size;
}

// Returns a reference to the element at index `i`, or throws an exception if out of range.
int& at(Vector& v, int i) {
    if (i < 0 || i >= v.size) {
        throw std::out_of_range("Index out of range");
        
    }
    return v.data[i];
}

// Returns a reference to the last element in the vector, or throws an exception if empty.
int& back(Vector& v) {
    if (v.size == 0) {
        throw std::out_of_range("Vector is empty");
    }
    return v.data[v.size - 1];
}

// Returns the capacity (allocated storage) of the vector.
int capacity(const Vector& v) {
    return v.capacity;
}

// Erases the elements of the vector but does not change capacity.
void clear(Vector& v) {
    v.size = 0;
}

// Checks if the vector is empty.
bool empty(const Vector& v) {
    return v.size == 0;
}

// Returns a reference to the first element in the vector, or throws an exception if empty.
int& front(Vector& v) {
    if (v.size == 0) {
        throw std::out_of_range("Vector is empty");
    }
    return v.data[0];
}

// Adds an element to the end of the vector, resizing if necessary.
void push_back(Vector& v, int element) {
    if (v.size == v.capacity) {
        // Increase capacity if needed
        v.capacity = v.capacity * 2;
        int* data_new = new int[v.capacity];
                for (int i = 0; i < v.size; ++i) {
                    data_new[i] = v.data[i];
                }
                delete[] v.data;
                v.data = data_new;
    }
    v.data[v.size++] = element;
}

// Searches for a key in the vector and returns its index, or -1 if not found.
int find(Vector& v, int key) {
    for (int i = 0; i < v.size; ++i) {
        if (v.data[i] == key) {
            return i;
        }
    }
    return -1; // Not found
}

// Sorts the vector in ascending order.
void sort(Vector& v) {
    std::sort(v.data, v.data + v.size);
}

// Multiplies two vectors element-wise and returns a new vector with the result.
// Assumes vectors are of the same size, otherwise throws an exception.
Vector multiply(Vector& a, Vector& b) {
    if (a.size != b.size) {
        throw std::invalid_argument("Vectors must be of the same size for multiplication.");
    }
    Vector newV;
    construct_vector(newV, a.size);
    for (int i = 0; i < a.size; ++i) {
        newV.data[i] = a.data[i] * b.data[i];
    }
    return newV;
}


