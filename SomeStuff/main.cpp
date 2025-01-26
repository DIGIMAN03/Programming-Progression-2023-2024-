
// Michael Helfand
// 1653490
// 12/4/24.
// A program that uses the class some stuff and the interface something to preform operations
#include "SomeStuffT.hpp"
#include <string>
#include <iostream>
#include <stdexcept>
/*
void info(CS52::SomeStuff &c) {
    std::cout << "Size is: " << c.size() << "\n";
    std::cout << "Capacity is: " << c.capacity() << "\n";
    std::cout << "Contents: ";
    for (int i = 0; i < c.size(); i++) {
        std::cout << c[i] << " ";
    }
    std::cout << "\n";
} // info
void info(CS52::SomeStuff &c) {
    std::cout << "Size is: " << c.size() << "\n";
    std::cout << "Capacity is: " << c.capacity() << "\n";
    std::cout << "Contents: ";
    for (int i = 0; i < c.size(); i++) {
        std::cout << c[i] << " ";
    }
    std::cout << "\n";
} // info
 */

template <class T>
void info(CS52::SomeStuff<T>& somestuff) {
    std::cout << "Size is: " << somestuff.size() << "\n";
    std::cout << "Capacity is: " << somestuff.capacity() << "\n";
    std::cout << "Contents: ";
    for (int i = 0; i < somestuff.size(); i++) {
        std::cout<< somestuff.at(i) << " ";
    }
    std::cout << "\n";
}

int main(int argc, const char * argv[]) {
    // Test with double
    CS52::SomeStuff<double> some_doubles(5, -1.0);
    std::cout << "CS52::SomeStuff<double> some_doubles(5, -1.0);\n";
    info(some_doubles);

    try {
        // Throw an exception for out-of-range access
        some_doubles.at(-1);
    } catch (const std::exception& e) {
        std::cout << "\nException handling:\n" << e.what() << "\n";
    }

    // Test with bool
    CS52::SomeStuff<bool> some_bools(5, true);
    std::cout << "\nCS52::SomeStuff<bool> some_bools(5, true);\n";
    some_bools.push_back(false);
    info(some_bools);

    // Test with std::string
    std::cout << "\nTesting std::string as an element:\n";
    CS52::SomeStuff<std::string> some_strings;
    some_strings.push_back("this is great");
    some_strings.push_back("I love stuff!!!");
    some_strings.push_back("Amazing stuff!");
    for ( int i = 0; i < some_strings.size(); i ++) {
    std :: cout << some_strings[ i ] << "\n" ;}

    // Test with a 2D matrix of int
    std::cout << "\nTesting a 2D int matrix of stuff:\n";
    CS52::SomeStuff<CS52::SomeStuff<int>> matrix_of_stuff;
    int rows = 5, cols = 3;
    for (int i = 0; i < rows; i++) {
        CS52::SomeStuff<int> some_ints(cols, i);
        matrix_of_stuff.push_back(some_ints);
    }

    for (int i = 0; i < matrix_of_stuff.size(); i++) {
        CS52::SomeStuff<int> some_ints = matrix_of_stuff.at(i);
        for (int j = 0; j < some_ints.size(); j++) {
            std::cout << some_ints.at(j) << " ";
        }
        std::cout << "\n";
    }
    CS52 :: SomeStuff < std::pair<std::string, int>> pair_of_somestuff (3 ,std::pair<std::string, int> {"year", 2024});
  
    std::cout <<"Size is "<< pair_of_somestuff.size();
    std :: cout << " Contents : \n" ;
    for ( int i = 0; i < pair_of_somestuff.size(); i ++){
        std::pair<std::string,int> p = pair_of_somestuff[i];
        std :: cout << p.first << " " << p.second + i << "\n" ;
    }
    std :: cout << "\n" ;
    std :: cout << "au revoir :)\n" ;
    char stop ; std :: cin >> stop ;
    // Default constructor, push_back, at methods
    /*
       std::cout << "// Default constructor, push_back, at\nCS52::SomeStuff a;\n";
       CS52::SomeStuff a;
       std::cout << "\na.push_back(10); a.push_back(88)\n";
       a.push_back(10);
       a.push_back(88);
       std::cout << "a.at(0) = 99;\n";
       a.at(0) = 99;
       info(a);

       // b
       std::cout << "\n// Overloaded constructor, [] op, at, empty, clear, exception handling: at\nCS52::SomeStuff b(2, 5);\n";
       CS52::SomeStuff b(2, 5);
       std::cout << "b is " << b << "\n";

       std::cout << "\n// Add more elements to b\n";
       std::cout << "b.push_back(10); b.push_back(2);\n";
       std::cout << "b.push_back(99); b.push_back(-5);\n";
       b.push_back(10);
       b.push_back(2);
       b.push_back(99);
       b.push_back(-5);

       std::cout << "\n// Array index [] and at()\n";
       std::cout << "b[0] = 25; b[1] = 1;\n";
       b[0] = 25;
       b[1] = 1;
       std::cout << "b.at(0); b.at(1);\n";
       std::cout << b.at(0) << " " << b.at(1) << "\n";

       std::cout << "\n\n// empty method, size, and capacity\n";
       std::cout << "b.empty() " << (b.empty() ? "True" : "False") << "\n";
       std::cout << "b.clear()\n";
       b.clear();
       info(b);
       std::cout << "b.empty() " << (b.empty() ? "True" : "False") << "\n";

       std::cout << "\n\n// Exception handling:\n";
       try {
           std::cout << "\nb.at(9) = ";
           std::cout << b.at(9);
       } catch (std::string msg) {
           std::cerr << "\n" << msg << std::endl;
       }

       // c
       std::cout << "\n// Copy constructor, copy assignment, pop_back,\n// capacity, size\nCS52::SomeStuff c(b);\n";
       std::cout << "\n// Add more elements to b\n";
       std::cout << "b.push_back(11); b.push_back(7);\n";
       std::cout << "b.push_back(3); b.push_back(23);\n";
       b.push_back(11);
       b.push_back(7);
       b.push_back(3);
       b.push_back(23);
       CS52::SomeStuff c(b);
       info(c);

       // d
       std::cout << "\n// Copy assignment\n";
       CS52::SomeStuff d;
       d = c;
       std::cout << "CS52::SomeStuff d = c;\n";
       std::cout << "d is " << d << "\n";

       std::cout << "\n// Size vs capacity()\n";
       std::cout << "d.size() is " << d.size() << "\n";
       std::cout << "d.capacity() is " << d.capacity() << "\n";

       std::cout << "\n// pop_back()\n";
       std::cout << "d.pop_back(); d.pop_back();\n";
       d.pop_back();
       d.pop_back();

       std::cout << "\n// Size vs capacity\n";
       std::cout << "d.size() is " << d.size() << "\n";
       std::cout << "d.capacity() is " << d.capacity() << "\n";

       // Test overloaded + operator
       CS52::SomeStuff m(5, 3);
       CS52::SomeStuff n(3, 0);
       n[0] = 55;
       n[1] = 100;
       n[2] = 500;

       CS52::SomeStuff o = m + n;

       std::cout << "Test overloaded + operator\n";
       std::cout << "m = [3, 3, 3, 3, 3], n = [55, 100, 500]\n";
       std::cout << "o = m + n = [" << o << "]\n";

       const int SIZE = 2;
       CS52::Something *something[SIZE];
       something[0] = new CS52::SomeStuff(2, 1);
       something[1] = new CS52::SomeStuff(5, 3);

       // Print out the type and length of each Something element
       std::cout << "\nTest polymorphism\n";
       for (int i = 0; i < SIZE; i++) {
           std::cout << "something[" << i << "] is a " << something[i]->type()
                     << " : " << *(dynamic_cast<CS52::SomeStuff *>(something[i])) << "\n";
       }
       for (int i = 0; i < 2; i++) {
           delete something[i];
       }

       char stop;
       std::cin >> stop;*/
       return 0;
}

