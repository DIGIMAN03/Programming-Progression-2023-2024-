// Michael Helfand
// 1653490
// 12/4/24.
// A program that uses the class some stuff and the interface something to preform operations on sets of data

#include "SomeStuff.h"
#include <string>
#include <iostream>
#include <stdexcept>

namespace CS52 {

// Default constructor
SomeStuff::SomeStuff() : _size(0), _capacity(0), _data(nullptr) {}

// Overloaded constructor
SomeStuff::SomeStuff(int size, int int_val)
    : _size(size), _capacity(size), _data(new int[size]) {
    for (int i = 0; i < _size; ++i) {
        _data[i] = int_val;
    }
}

// Copy constructor
SomeStuff::SomeStuff(const SomeStuff &other)
    : _size(other._size), _capacity(other._capacity), _data(new int[other._capacity]) {
    for (int i = 0; i < _size; ++i) {
        _data[i] = other._data[i];
    }
}

// Destructor
SomeStuff::~SomeStuff() {
    delete[] _data;
}

// Access element at index
int& SomeStuff::at(int i) const {
    if(i<0 || i>=_size)
        throw std::string("out-of-bounds");
    
    return _data[i];
}

// Clear all elements
void SomeStuff::clear() {
    _size = 0;
}

// Get raw pointer to data
int* SomeStuff::data() const {
    return _data;
}

// Check if empty
bool SomeStuff::empty() const {
    return _size == 0;
}

// Remove last element
void SomeStuff::pop_back() {
    if (_size > 0) --_size;
}

// Add element to the end
void SomeStuff::push_back(int element) {
    if (_size == _capacity) {
        int new_capacity = (_capacity == 0) ? 1 : _capacity * 2;
        int* new_data = new int[new_capacity];

        for (int i = 0; i < _size; ++i) {
            new_data[i] = _data[i];
        }

        delete[] _data;
        _data = new_data;
        _capacity = new_capacity;
    }
    _data[_size++] = element;
}

// Override size
int SomeStuff::size() const {
    return _size;
}

// Override capacity
int SomeStuff::capacity() const {
    return _capacity;
}

// Override type
std::string SomeStuff::type() const {
    return "CS52::SomeStuff";
}

// Overload operator[]
int& SomeStuff::operator[](int index) {
    return at(index);
}

// Overload assignment operator
SomeStuff& SomeStuff::operator=(const SomeStuff &other) {
    if (this != &other) {
        delete[] _data;
        _size = other._size;
        _capacity = other._capacity;
        _data = new int[_capacity];
        for (int i = 0; i < _size; ++i) {
            _data[i] = other._data[i];
        }
    }
    return *this;
}

// Overload addition operator
SomeStuff SomeStuff::operator+(const SomeStuff &rhs) {
    int new_size = std::min(_size, rhs._size);
    SomeStuff result(new_size, 0);
    for (int i = 0; i < new_size; ++i) {
        result._data[i] = this->_data[i] + rhs._data[i];
    }
    return result;
}

// Overload stream insertion operator
std::ostream& operator<<(std::ostream &os, SomeStuff &obj) {
    for (int i = 0; i < obj.size(); ++i) {
        os << obj[i] << " ";
    }
    return os;
}

} // namespace CS52


