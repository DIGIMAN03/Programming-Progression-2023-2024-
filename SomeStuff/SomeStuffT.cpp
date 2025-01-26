// MICHAEL HELFAND
// ID: 1653490
// Assignment 5
// 12/07/2024
// This program is a template class SomeStuff taht allwos us to create and manipulate data of any given type such as add new bits of data or clear it or find specified pieces of data all using the class SomeStuff

#ifndef SOMESTUFFT_HPP
#define SOMESTUFFT_HPP

#include <iostream>
#include <stdexcept> // For std::out_of_range
#include <string>    // For std::string if needed

namespace CS52 {

template <typename T>
class SomeStuff {
public:
    // Constructors
    SomeStuff();                              // Default constructor
    SomeStuff(int size, T init_val = T());    // Overloaded constructor
    SomeStuff(const SomeStuff<T>& other);    // Copy constructor

    // Destructor
    ~SomeStuff();

    // Public member functions
    T& at(int i) const;                       // Access element with bounds checking
    void clear();                             // Clear the container
    T* data() const;                          // Return raw data pointer
    bool empty() const;                       // Check if container is empty
    void pop_back();                          // Remove the last element
    void push_back(const T& element);         // Add a new element
    int size() const;                         // Return current size
    int capacity() const;                     // Return current capacity

    // Overloaded operators
    T& operator[](int index);                 // Access element by index
    SomeStuff<T>& operator=(const SomeStuff<T>& other); // Assignment operator

private:
    int _size = 0;
    int _capacity = 0;
    T* _data = nullptr;

    // Helper function to resize the internal storage
    void resize(int new_capacity);
};

// Implementation

// Default constructor
template <typename T>
SomeStuff<T>::SomeStuff() : _size(0), _capacity(0), _data(nullptr) {}

// Overloaded constructor
template <typename T>
SomeStuff<T>::SomeStuff(int size, T init_val)
    : _size(size), _capacity(size), _data(new T[size]) {
    for (int i = 0; i < size; ++i) {
        _data[i] = init_val;
    }
}

// Copy constructor
template <typename T>
SomeStuff<T>::SomeStuff(const SomeStuff<T>& other)
    : _size(other._size), _capacity(other._capacity), _data(new T[other._capacity]) {
    for (int i = 0; i < _size; ++i) {
        _data[i] = other._data[i];
    }
}

// Destructor
template <typename T>
SomeStuff<T>::~SomeStuff() {
    delete[] _data;
}

// Access element with bounds checking
template <typename T>
T& SomeStuff<T>::at(int i) const {
    if (i < 0 || i >= _size) {
        throw std::out_of_range("out-of-bounds");
    }
    return _data[i];
}

// Clear the container
template <typename T>
void SomeStuff<T>::clear() {
    delete[] _data;
    _data = nullptr;
    _size = 0;
    _capacity = 0;
}

// Return raw data pointer
template <typename T>
T* SomeStuff<T>::data() const {
    return _data;
}

// Check if container is empty
template <typename T>
bool SomeStuff<T>::empty() const {
    return _size == 0;
}

// Remove the last element
template <typename T>
void SomeStuff<T>::pop_back() {
    if (_size > 0) {
        --_size;
    }
}

// Add a new element
template <typename T>
void SomeStuff<T>::push_back(const T& element) {
    if (_size == _capacity) {
        resize((_capacity == 0) ? 1 : _capacity * 2);
    }
    _data[_size++] = element;
}

// Return current size
template <typename T>
int SomeStuff<T>::size() const {
    return _size;
}

// Return current capacity
template <typename T>
int SomeStuff<T>::capacity() const {
    return _capacity;
}

// Access element by index
template <typename T>
T& SomeStuff<T>::operator[](int index) {
    return _data[index];
}

// Assignment operator
template <typename T>
SomeStuff<T>& SomeStuff<T>::operator=(const SomeStuff<T>& other) {
    if (this != &other) {
        delete[] _data;
        _size = other._size;
        _capacity = other._capacity;
        _data = new T[_capacity];
        for (int i = 0; i < _size; ++i) {
            _data[i] = other._data[i];
        }
    }
    return *this;
}

// Resize the internal storage
template <typename T>
void SomeStuff<T>::resize(int new_capacity) {
    T* new_data = new T[new_capacity];
    for (int i = 0; i < _size; ++i) {
        new_data[i] = _data[i];
    }
    delete[] _data;
    _data = new_data;
    _capacity = new_capacity;
}

} // namespace CS52

#endif // SOMESTUFFT_HPP
