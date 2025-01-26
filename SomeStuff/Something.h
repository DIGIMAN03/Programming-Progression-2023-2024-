//
//  Something.h
//  SomeStuff
//
//  Created by Michael Helfand on 12/4/24.
//

//Something class Declaration
#ifndef SOMETHING_H
#define SOMETHING_H

#include <string>

namespace CS52 {

class Something
{
public :
virtual int size () const = 0;
virtual int capacity () const = 0;
virtual std :: string type () const = 0;
virtual ~ Something () {}
};}
#endif
 
