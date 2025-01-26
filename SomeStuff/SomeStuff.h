// Michael Helfand
// 1653490
// 12/4/24.
// A program that uses the class some stuff and the interface something to preform operations on sets of data



#ifndef SOMESTUFF_H
#define SOMESTUFF_H

#include "Something.h" // Include the base class

#include <iostream>
#include <string>
#include <stdexcept> // for exceptions
namespace CS52 {

class SomeStuff : public Something {
public:
    // Constructors
    SomeStuff();                          // Default constructor
    SomeStuff(int size, int int_val);     // Overloaded constructor
    SomeStuff(const SomeStuff &other);    // Copy constructor

    // Destructor
    ~SomeStuff();

    // Public member functions
    int& at(int i) const;
    void clear();
    int* data() const;
    bool empty() const;
    void pop_back();
    void push_back(int element);

    // Overrides
    int size() const override;
    int capacity() const override;
    std::string type() const override;

    // Overloaded operators
    int& operator[](int index);
    SomeStuff& operator=(const SomeStuff &other);
    SomeStuff operator+(const SomeStuff &rhs);
    friend std::ostream& operator<<(std::ostream &os, SomeStuff &obj);

private:
    int _size = 0;
    int _capacity = 0;
    int* _data = nullptr;

};

} // namespace CS52

#endif

